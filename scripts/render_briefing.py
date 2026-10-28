"""Render docs/briefing/BRIEFING.md to docs/briefing/BRIEFING.pdf (markdown -> HTML -> xhtml2pdf).

Requires `pip install markdown xhtml2pdf`. Fonts: DejaVu Sans from matplotlib (covers ≥, −, ≈, µ, α).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import markdown
import matplotlib
from reportlab.lib.fonts import addMapping
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from xhtml2pdf import pisa
from xhtml2pdf.default import DEFAULT_FONT

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "docs/briefing/BRIEFING.md"
OUT = ROOT / "docs/briefing/BRIEFING.pdf"
FONTS = Path(matplotlib.get_data_path()) / "fonts/ttf"

CSS = f"""
@page {{ size: a4; margin: 2cm 2cm 2cm 2cm;
        @frame footer {{ -pdf-frame-content: footer; bottom: 0.9cm; margin-left: 2cm; margin-right: 2cm; height: 0.6cm; }} }}
body {{ font-family: DejaVu; font-size: 9.6pt; line-height: 1.38; color: #1b1b1a; }}
h1 {{ font-size: 18pt; color: #0b0b0b; margin-bottom: 2pt; line-height: 1.2; }}
h2 {{ font-size: 13pt; color: #1d5fae; margin-top: 14pt; margin-bottom: 4pt; -pdf-keep-with-next: true; }}
h3 {{ font-size: 11pt; color: #0b0b0b; margin-top: 9pt; margin-bottom: 3pt; -pdf-keep-with-next: true; }}
p {{ margin: 0 0 6pt 0; text-align: justify; }}
li {{ margin-bottom: 2pt; }}
em {{ color: #3d3c39; }}
code {{ font-family: DejaVuMono; font-size: 8.3pt; color: #3d3c39; }}
img {{ margin-top: 6pt; }}
table {{ border-top: 0.8pt solid #52514e; border-bottom: 0.8pt solid #52514e; margin: 4pt 0 8pt 0; }}
th {{ font-size: 8.4pt; text-align: left; padding: 3pt 4pt 2pt 4pt; border-bottom: 0.6pt solid #52514e; background-color: #f2f1ed; }}
td {{ font-size: 8.4pt; padding: 3pt 4pt 2pt 4pt; border-bottom: 0.3pt solid #d6d5d0; vertical-align: top; }}
blockquote {{ background-color: #f4f3ef; border-left: 2.5pt solid #9b9a95; padding: 5pt 8pt 2pt 8pt; margin: 6pt 0 8pt 0; }}
hr {{ color: #d6d5d0; }}
#footer {{ font-size: 7.5pt; color: #6d6c68; text-align: center; }}
"""


def register_fonts():
    """Register DejaVu (regular/bold/italic) with reportlab and expose it to xhtml2pdf's CSS font lookup."""
    files = {"DejaVu": "DejaVuSans.ttf", "DejaVu-Bold": "DejaVuSans-Bold.ttf", "DejaVu-Oblique": "DejaVuSans-Oblique.ttf",
             "DejaVu-BoldOblique": "DejaVuSans-BoldOblique.ttf", "DejaVuMono": "DejaVuSansMono.ttf"}
    for name, f in files.items():
        pdfmetrics.registerFont(TTFont(name, str(FONTS / f)))
    addMapping("DejaVu", 0, 0, "DejaVu"); addMapping("DejaVu", 1, 0, "DejaVu-Bold")
    addMapping("DejaVu", 0, 1, "DejaVu-Oblique"); addMapping("DejaVu", 1, 1, "DejaVu-BoldOblique")
    for v in (0, 1):
        for w in (0, 1):
            addMapping("DejaVuMono", v, w, "DejaVuMono")
    DEFAULT_FONT["dejavu"] = "DejaVu"
    DEFAULT_FONT["dejavumono"] = "DejaVuMono"


LIST = re.compile(r"^(\s*)([-*]|\d+\.)\s")


def normalise_markdown(md: str) -> str:
    """GitHub-style -> python-markdown: blank line before a list that follows a paragraph line, and nested
    indentation (2-3 spaces per level) widened to 4 spaces per level."""
    out, prev = [], ""
    for line in md.splitlines():
        m = LIST.match(line)
        lead = len(line) - len(line.lstrip(" "))
        if line.strip() and lead and not line.lstrip().startswith("|"):
            line = " " * (4 * ((lead + 1) // 3 if lead < 6 else 2)) + line.lstrip(" ")
        if m and not m.group(1) and prev.strip() and not LIST.match(prev) and not prev.startswith(" "):
            out.append("")
        out.append(line)
        prev = line
    return "\n".join(out)


def link_callback(uri, rel):
    p = (SRC.parent / uri).resolve()
    return str(p) if p.exists() else uri


def main() -> int:
    register_fonts()
    md = normalise_markdown(SRC.read_text())
    body = markdown.markdown(md, extensions=["tables", "sane_lists"])
    # figures: full text width
    body = body.replace('<img alt="Figure', '<img width="480" alt="Figure')
    html = (f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>"
            f"<div id='footer'>Agentic kinase activity inference — research briefing · page <pdf:pagenumber/></div>"
            f"{body}</body></html>")
    with open(OUT, "wb") as fh:
        res = pisa.CreatePDF(html, dest=fh, link_callback=link_callback, encoding="utf-8")
    print(("FAILED " if res.err else "wrote ") + str(OUT))
    return 1 if res.err else 0


if __name__ == "__main__":
    sys.exit(main())
