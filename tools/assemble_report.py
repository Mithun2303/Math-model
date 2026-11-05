"""Assemble results/COMPLETE_RESULTS_REPORT.md from part files; resolve [@key] citations against report/references.bib."""
import re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
PARTS = ROOT / "results" / "_report_parts"
ORDER = ["front.md", "s01_06.md", "s07_08.md", "s09.md", "s10_11.md", "s12.md", "s13.md", "s14_18.md", "s19a.md", "s19b.md", "s19c.md", "s20_appendix.md"]

def parse_bib(path):
    txt = path.read_text()
    txt = "\n".join(l for l in txt.splitlines() if not l.lstrip().startswith("%"))
    ents = {}
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", txt, flags=re.S):
        typ, key, body = m.groups()
        f = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}\s*,?\s*\n", body + "\n", flags=re.S):
            f[fm.group(1).lower()] = re.sub(r"\s+", " ", fm.group(2)).strip()
        f["_type"] = typ
        ents[key.strip()] = f
    return ents

ACC = {r"{\'a}": "á", r"{\~n}": "ñ", r"{\'e}": "é", r"{\'o}": "ó", r"{\'i}": "í", r"{\"o}": "ö", r"{\"u}": "ü", r"{\'u}": "ú", r"{\'y}": "ý",
       r"{\v{c}}": "č", r"{\'{a}}": "á", r"{\'{e}}": "é", r"{\'{o}}": "ó", r"{\'{i}}": "í", r"{\c{c}}": "ç", r"{\aa}": "å", r"{\o}": "ø"}
def clean(s):
    for k, v in ACC.items():
        s = s.replace(k, v)
    s = s.replace("--", "–").replace("{", "").replace("}", "").replace("\\&", "&").replace("\\", "")
    return s

def surname_list(a):
    names = [x.strip() for x in a.split(" and ")]
    sn = [clean(n.split(",")[0]) for n in names]
    return sn

def cite_label(e):
    sn = surname_list(e.get("author", e.get("title", "?")))
    y = e.get("year", "n.d.")
    if len(sn) == 1: return f"{sn[0]} {y}"
    if len(sn) == 2: return f"{sn[0]} and {sn[1]} {y}"
    return f"{sn[0]} et al. {y}"

def ref_entry(e):
    au = e.get("author", "")
    names = []
    for n in [x.strip() for x in au.split(" and ")]:
        if "," in n:
            s, g = n.split(",", 1); names.append(f"{clean(s.strip())}, {''.join(w[0] + '.' for w in clean(g).replace('-', ' ').split())}")
        else: names.append(clean(n))
    s = "; ".join(names) if names else ""
    out = f"{s} ({e.get('year','n.d.')}). {clean(e.get('title',''))}."
    if e.get("journal"):
        out += f" *{clean(e['journal'])}*"
        if e.get("volume"): out += f" **{e['volume']}**"
        if e.get("number"): out += f"({e['number']})"
        if e.get("pages"): out += f": {clean(e['pages'])}"
        out += "."
    elif e.get("booktitle"): out += f" In: *{clean(e['booktitle'])}*."
    elif e.get("howpublished") or e.get("publisher"): out += f" {clean(e.get('howpublished', e.get('publisher','')))}."
    if e.get("doi"): out += f" doi:{e['doi']}"
    elif e.get("url"): out += f" {e['url']}"
    return out

def main():
    bib = parse_bib(ROOT / "report" / "references.bib")
    cited = []
    def sub(m):
        keys = [k.strip().lstrip("@") for k in m.group(1).split(";")]
        labs = []
        for k in keys:
            if k not in bib:
                print("UNKNOWN CITATION KEY:", k, file=sys.stderr); labs.append(f"??{k}??"); continue
            if k not in cited: cited.append(k)
            labs.append(cite_label(bib[k]))
        return "(" + "; ".join(labs) + ")"
    chunks = []
    for name in ORDER:
        p = PARTS / name
        if p.exists():
            chunks.append(p.read_text().rstrip() + "\n")
        else:
            print("missing part:", name, file=sys.stderr)
    body = "\n\n".join(chunks)
    body = re.sub(r"\[((?:@[\w:.-]+\s*;?\s*)+)\]", sub, body)
    # reference section placeholder
    refs = sorted(cited, key=lambda k: cite_label(bib[k]).lower())
    reflist = "\n\n".join(f"{i+1}. {ref_entry(bib[k])}  " for i, k in enumerate(refs))
    reflist = "\n".join(f"{i+1}. {ref_entry(bib[k])}" for i, k in enumerate(refs))
    body = body.replace("<!--REFERENCES-->", reflist)
    out = ROOT / "results" / "COMPLETE_RESULTS_REPORT.md"
    out.write_text(body)
    print("wrote", out, len(body.split()), "words;", len(cited), "cited of", len(bib), "bib entries; uncited:", sorted(set(bib) - set(cited)))
main()
