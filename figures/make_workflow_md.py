"""Write WORKFLOW.md and the standalone figure SVGs from figures/workflow.html.

    .venv/bin/python figures/make_workflow_md.py

figures/workflow.html is the source. GitHub cannot show the page's inline
SVGs in Markdown, so each one is written as a standalone file (the page's
light-theme colours and fonts filled in) and WORKFLOW.md, at the repository
root, links them. The text, tables and captions are converted as they are.
Run it again after any edit to the page.
"""

import html
import os
import re
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PAGE = os.path.join(HERE, "workflow.html")
OUT = os.path.join(ROOT, "WORKFLOW.md")
# The page's figures, in order: (file name, alt text, write it here). The map
# is written by make_coverage_map.py, with its land-outline attribution.
FIGURES = [("fig1_reference_build", "How the reference is built", True),
           ("fig2_coverage_map", "Coverage map of the reference", False),
           ("fig3_prediction_pipeline", "How one curve becomes a prediction", True),
           ("fig4_validation_scale", "Accuracy at each validation level", True)]
FONTS = {"sans": "IBM Plex Sans, Helvetica, Arial, sans-serif",
         "serif": "Source Serif 4, Georgia, serif",
         "mono": "IBM Plex Mono, Menlo, monospace"}
WIDTH = 78


def palette(page):
    """{css variable: value} from the page's light-theme :root block."""
    root = re.search(r":root\s*{([^}]*)}", page).group(1)
    return dict(re.findall(r"--([\w-]+):\s*([^;]+);", root))


def standalone(svg, colours):
    """An inline page SVG as a self-contained file."""
    for name, value in colours.items():
        if name not in FONTS:
            svg = svg.replace(f"var(--{name})", value.strip())
    for name, stack in FONTS.items():
        svg = svg.replace(f"var(--{name})", stack)
    svg = svg.replace("currentColor", colours["ink"].strip())
    _, _, w, h = re.search(r'viewBox="([^"]+)"', svg).group(1).split()
    svg = svg.replace("<svg ", f'<svg xmlns="http://www.w3.org/2000/svg" '
                               f'width="{float(w):.0f}" height="{float(h):.0f}" ', 1)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            "<!-- py_HB_texture, written from figures/workflow.html by "
            "figures/make_workflow_md.py -->\n" + svg + "\n")


def inline(s):
    """Inline HTML to Markdown."""
    s = re.sub(r'<span class="m">(.*?)</span>', r"`\1`", s, flags=re.S)
    s = re.sub(r"<code>(.*?)</code>", r"`\1`", s, flags=re.S)
    s = re.sub(r"<em>(.*?)</em>", r"*\1*", s, flags=re.S)
    s = re.sub(r"<strong>(.*?)</strong>", r"**\1**", s, flags=re.S)
    s = re.sub(r"<(?!/?sub\b)[^>]+>", "", s)          # keep <sub>, drop the rest
    s = html.unescape(s).replace(" ", " ")
    return re.sub(r"\s+", " ", s).strip()


def para(s):
    lines = textwrap.wrap(inline(s), WIDTH, break_long_words=False,
                          break_on_hyphens=False)
    # a wrapped line must not start like a list item, heading or quote
    for i in range(1, len(lines)):
        while re.match(r"(\d+\.|[-+*]) |#|>", lines[i]):
            head, _, rest = lines[i].partition(" ")
            lines[i - 1] += " " + head
            lines[i] = rest
    return "\n".join(lines)


def table(t):
    rows = re.findall(r"<tr>(.*?)</tr>", t, flags=re.S)
    cells = [re.findall(r"<t([hd])([^>]*)>(.*?)</t[hd]>", r, flags=re.S) for r in rows]
    head = cells[0]
    out = ["| " + " | ".join(inline(c[2]) for c in head) + " |",
           "|" + "|".join("---:" if 'class="num"' in c[1] else "---" for c in head) + "|"]
    out += ["| " + " | ".join(inline(c[2]) for c in r) + " |" for r in cells[1:]]
    return "\n".join(out)


def main():
    page = open(PAGE, encoding="utf-8").read()
    colours = palette(page)
    svgs = re.findall(r"<svg\b.*?</svg>", page, flags=re.S)
    assert len(svgs) == len(FIGURES), f"{len(svgs)} figures on the page"
    for (name, _, write), svg in zip(FIGURES, svgs):
        if write:
            open(os.path.join(HERE, name + ".svg"), "w", encoding="utf-8").write(
                standalone(svg, colours))

    body = page[page.index('<div class="wrap">'):]
    body = re.sub(r"<svg\b.*?</svg>", "", body, flags=re.S)
    token = re.compile(
        r'<p class="eyebrow">.*?</p>|<h1>(?P<h1>.*?)</h1>'
        r'|<span class="stage-n">(?P<n>\d+)</span>\s*<h2>(?P<h2>.*?)</h2>'
        r"|<h3[^>]*>(?P<h3>.*?)</h3>|<figure>.*?<figcaption>(?P<cap>.*?)</figcaption>\s*</figure>"
        r'|<table>(?P<tab>.*?)</table>|<div class="note">(?P<note>.*?)</div>'
        r"|<footer>(?P<foot>.*?)</footer>|<p[^>]*>(?P<p>.*?)</p>", re.S)
    out = ["<!-- Written from figures/workflow.html by figures/make_workflow_md.py;"
           " edit the page, not this file. -->"]
    fig = 0
    for m in token.finditer(body):
        g = m.groupdict()
        if g["h1"]:
            out.append("# " + inline(g["h1"]))
        elif g["h2"]:
            out.append(f"## {int(g['n'])}. {inline(g['h2'])}")
        elif g["h3"]:
            out.append("### " + inline(g["h3"]))
        elif g["cap"] is not None:
            name, alt, _ = FIGURES[fig]
            fig += 1
            out.append(f"![{alt}](figures/{name}.svg)")
            out.append(para(f"<em>Figure {fig}.</em> " + g["cap"]))
        elif g["tab"]:
            out.append(table(g["tab"]))
        elif g["note"]:
            ps = re.findall(r"<p>(.*?)</p>", g["note"], flags=re.S)
            out.append("\n>\n".join("\n".join("> " + l for l in para(p).split("\n"))
                                    for p in ps))
        elif g["foot"]:
            out.append("---")
            out.append(para(g["foot"]))
        elif g["p"]:
            out.append(para(g["p"]))
    assert fig == len(FIGURES)
    open(OUT, "w", encoding="utf-8").write("\n\n".join(out) + "\n")
    print(f"wrote {os.path.relpath(OUT, ROOT)} and "
          f"{sum(f[2] for f in FIGURES)} figure SVGs")


if __name__ == "__main__":
    main()
