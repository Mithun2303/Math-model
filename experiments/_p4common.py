"""Shared helpers for the Phase-4 experiments: config, provenance stamp, figure footer with a SYNTHETIC / CLINICAL marker."""
import hashlib
import sys
import textwrap
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _style import *  # noqa: F401,F403
import numpy as np
import pandas as pd
import yaml

CFG_PATH = ROOT / "configs" / "phase4.yaml"
CFG = yaml.safe_load(open(CFG_PATH))
CFG_HASH = hashlib.sha256(CFG_PATH.read_bytes()).hexdigest()[:12]
T0 = time.time()
SYN = "SYNTHETIC DATA from a hypothetical model (known ground truth, assumed observation model); not clinical evidence."
CLIN = "CLINICAL (descriptive): cBioPortal METABRIC export, RFS as defined there; no model fitted or overlaid."


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def truth(name):
    return {k: v for k, v in CFG["ground_truths"][name].items() if k != "note"}


def save_table(df, name, index=False):
    path = TAB / f"{name}.csv"
    with open(path, "w") as fh:
        fh.write(f"# phase4 config sha256[:12]={CFG_HASH}; synthetic/hypothetical unless the file name says metabric\n")
        df.to_csv(fh, index=index)
    return path


def p4_footer(fig, extra="", kind="SYN"):
    """Footer that marks the figure SYNTHETIC or CLINICAL and is wrapped to the figure width."""
    base = SYN if kind == "SYN" else CLIN
    txt = (extra + " " if extra else "") + f"{base} [config {CFG_HASH}]"
    lines = textwrap.wrap(txt, max(60, int(fig.get_figwidth() * 19)))
    if len(lines) > 1:
        fig.subplots_adjust(bottom=min(0.45, fig.subplotpars.bottom + 0.035 * (len(lines) - 1)))
    fig.text(0.01, 0.005, "\n".join(lines), fontsize=6.5, color=("#9b1c1c" if kind == "SYN" else "#1c4f9b"), ha="left", va="bottom")
    fig.text(0.995, 0.995, "SYNTHETIC" if kind == "SYN" else "CLINICAL (descriptive)", fontsize=8, fontweight="bold",
             color=("#9b1c1c" if kind == "SYN" else "#1c4f9b"), ha="right", va="top")
