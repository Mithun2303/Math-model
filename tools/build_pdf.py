"""Markdown -> HTML (python-markdown + KaTeX, server-side) -> A4 PDF via headless Chromium (Brave). Usage: build_pdf.py [md] [pdf]"""
import html, json, re, subprocess, sys, tempfile
from pathlib import Path
import markdown

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
BRAVE = "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser"
md_path = Path(sys.argv[1]) if len(sys.argv) > 1 else RES / "COMPLETE_RESULTS_REPORT.md"
pdf_path = Path(sys.argv[2]) if len(sys.argv) > 2 else RES / "COMPLETE_RESULTS_REPORT.pdf"
html_path = RES / "COMPLETE_RESULTS_REPORT.html"

text = md_path.read_text()

# ---- protect fenced code, then math
codes = []
def keep_code(m):
    codes.append(m.group(0)); return f"CODEPH{len(codes)-1}END"
text = re.sub(r"^```.*?^```", keep_code, text, flags=re.S | re.M)

maths = []
def keep_disp(m):
    maths.append((m.group(1).strip(), True)); return f"\n\nMATHPH{len(maths)-1}END\n\n"
text = re.sub(r"\$\$(.+?)\$\$", keep_disp, text, flags=re.S)
def keep_inl(m):
    maths.append((m.group(1), False)); return f"MATHPH{len(maths)-1}END"
text = re.sub(r"(?<![\\$\w])\$(?!\s)([^$\n]+?)(?<!\s)\$(?![\w$])", keep_inl, text)

# ---- render math with KaTeX
with tempfile.TemporaryDirectory() as td:
    jin, jout = Path(td, "in.json"), Path(td, "out.json")
    jin.write_text(json.dumps([{"tex": t, "display": d} for t, d in maths]))
    subprocess.run(["node", str(ROOT / "tools" / "katex_render.js"), str(jin), str(jout)], check=True)
    rendered = json.loads(jout.read_text())
bad = [(i, r) for i, r in enumerate(rendered) if not r["ok"]]
for i, r in bad:
    print("KATEX ERROR:", r["err"], "| TeX:", r["tex"][:120])
print(f"math expressions: {len(maths)}, katex errors: {len(bad)}")

for i, c in enumerate(codes):
    text = text.replace(f"CODEPH{i}END", c)

body = markdown.markdown(text, extensions=["tables", "fenced_code", "sane_lists"], output_format="html5")

def put_math(m):
    i = int(m.group(1)); r = rendered[i]
    inner = r["html"] if r["ok"] else f'<span class="mathbad">{html.escape(r["tex"])}</span>'
    return inner
body = re.sub(r"<p>\s*MATHPH(\d+)END\s*</p>", lambda m: '<div class="mathd">' + put_math(m) + "</div>", body)
body = re.sub(r"MATHPH(\d+)END", put_math, body)

# ---- headings: ids + ToC (levels 2,3)
toc, used = [], set()
def head(m):
    lvl, inner = int(m.group(1)), m.group(2)
    plain = re.sub(r"<[^>]+>", "", inner)
    slug = re.sub(r"[^a-z0-9]+", "-", plain.lower()).strip("-")[:60] or "sec"
    base, k = slug, 2
    while slug in used: slug = f"{base}-{k}"; k += 1
    used.add(slug)
    if lvl in (2, 3): toc.append((lvl, slug, plain))
    return f'<h{lvl} id="{slug}">{inner}</h{lvl}>'
body = re.sub(r"<h([1-4])>(.*?)</h\1>", head, body)
toc_html = '<nav class="toc"><h2 class="notoc">Table of Contents</h2>' + "".join(
    f'<div class="t{l}"><a href="#{s}">{t}</a></div>' for l, s, t in toc) + "</nav>"
body = body.replace("<!--TOC-->", toc_html)

# ---- figures and tables
missing = []
nfig = 0
def fig(m):
    global nfig
    alt, src, cap = html.unescape(m.group(1)), m.group(2), m.group(3)
    if not (RES / src).exists(): missing.append(src)
    nfig += 1
    return f'<figure><img src="{src}" alt="{html.escape(alt)}"><figcaption>{cap}</figcaption></figure>'
body = re.sub(r'<p><img alt="([^"]*)" src="([^"]+)"\s*/?></p>\s*<p>(<strong>Figure.*?)</p>', fig, body, flags=re.S)
body = re.sub(r'<p>(<strong>Table [^<]*</strong>.*?)</p>\s*(<table>.*?</table>)', r'<div class="tbl"><p class="tcap">\1</p>\2</div>', body, flags=re.S)
body = body.replace('<table>', '<table>')
ntab = len(re.findall(r'<div class="tbl">', body))
print(f"figures embedded: {nfig}, tables captioned: {ntab}, missing images: {missing}")
stray = re.findall(r'<img [^>]*>', body)
print("img tags total:", len(stray))

css = """
@page { size: A4; margin: 22mm 18mm 22mm 18mm; @bottom-center { content: counter(page); font: 9pt Georgia, serif; color:#444 }
  @top-right { content: string(sect); font: italic 8pt Georgia, serif; color:#777 } }
@page :first { @bottom-center { content: none } @top-right { content: none } }
html { font: 10.5pt/1.42 Georgia, 'Times New Roman', serif; color:#111 }
body { margin:0; hyphens:auto; text-align:justify }
h1 { font: bold 22pt Georgia, serif; text-align:center; margin: 60mm 0 6mm; line-height:1.2 }
h2 { font: bold 15pt Georgia, serif; margin: 7mm 0 3mm; break-after:avoid; border-bottom:1px solid #999; padding-bottom:1mm; string-set: sect content() }
h2.notoc { border:none }
h3 { font: bold 12pt Georgia, serif; margin: 5mm 0 2mm; break-after:avoid }
h4 { font: bold italic 10.5pt Georgia, serif; margin: 4mm 0 1.5mm; break-after:avoid }
h2:not(.notoc) { break-before: page }
p { margin: 0 0 2.4mm; orphans:3; widows:3 }
.mathd { margin: 2mm 0; text-align:center; overflow-x:auto; overflow-y:hidden; font-size:1em; break-inside:avoid }
.katex { font-size:1.05em } .katex-display { margin:0 }
figure { margin: 4mm 0; break-inside:avoid; text-align:center }
figure img { max-width:100%; max-height: 120mm; width:auto; height:auto; border:1px solid #ddd }
figcaption { font-size:8.8pt; text-align:justify; margin-top:1.5mm; line-height:1.32; hyphens:auto }
.tbl { margin: 3mm 0 4mm; break-inside:avoid-page }
.tcap { font-size:9pt; margin:0 0 1mm; text-align:left; break-after:avoid }
table { border-collapse:collapse; width:100%; font-size:8.3pt; line-height:1.25; text-align:left; hyphens:manual }
th, td { border-top:0.5px solid #888; border-bottom:0.5px solid #888; padding:1mm 1.5mm; vertical-align:top; word-break:break-word; overflow-wrap:anywhere }
th { background:#eee; font-weight:bold }
td:not(:first-child), th:not(:first-child) { text-align:right }
code { font: 8.8pt Menlo, monospace; background:#f2f2f2; padding:0 1px; overflow-wrap:anywhere }
pre { font: 8pt/1.3 Menlo, monospace; background:#f6f6f6; padding:2mm; white-space:pre-wrap; break-inside:avoid; overflow-wrap:anywhere }
pre code { background:none }
blockquote { margin: 2mm 8mm; font-size:9.5pt; color:#333; border-left:2px solid #999; padding-left:3mm }
nav.toc { font-size:9.5pt; text-align:left } nav.toc a { color:#111; text-decoration:none }
nav.toc .t2 { font-weight:bold; margin-top:1.6mm } nav.toc .t3 { margin-left:7mm; font-size:9pt }
.meta { text-align:center; font-size:11pt; line-height:1.6 } .pagebreak { break-after:page }
.mathbad { color:#b00; font-family:Menlo,monospace }
ol { padding-left: 7mm } li { margin-bottom: 1.2mm } ol li { text-align:left; font-size:9.3pt }
"""
doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Predicting Cancer Relapse Through Mathematical Modelling of Tumour Dormancy</title>
<link rel="stylesheet" href="file://{ROOT}/tools/katex_css/katex.min.css"><style>{css}</style></head><body>{body}</body></html>"""
html_path.write_text(doc)
cmd = [BRAVE, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={pdf_path}", "--virtual-time-budget=20000",
       "--allow-file-access-from-files", f"file://{html_path}"]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
print("brave rc", r.returncode, r.stderr[-300:])
print("pdf exists:", pdf_path.exists(), pdf_path.stat().st_size if pdf_path.exists() else 0)
