"""Shared helpers for the Phase-3 experiments (figure style, provenance stamping, saving)."""
import hashlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _style import *  # noqa: F401,F403
import numpy as np
import pandas as pd

from dormancy.uncertainty import load_phase3_config

CFG_PATH = ROOT / "configs" / "phase3.yaml"
CFG = load_phase3_config(CFG_PATH)
CFG_HASH = hashlib.sha256(CFG_PATH.read_bytes()).hexdigest()[:12]
HYP = "Hypothetical simulated quantities (illustrative priors/ranges); not patient risk or clinical calibration."
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def save_table(df, name, index=False):
    """Write results/tables/<name>.csv; first line is a '# config=<hash>' provenance comment."""
    path = TAB / f"{name}.csv"
    with open(path, "w") as fh:
        fh.write(f"# phase3 config sha256[:12]={CFG_HASH}; hypothetical simulated quantities\n")
        df.to_csv(fh, index=index)
    return path


def p3_footer(fig, extra=""):
    """Footer with the hypothetical-assumption disclaimer; wrapped to the figure width so it is never cut off."""
    import textwrap
    txt = (extra + " " if extra else "") + f"{HYP} [config {CFG_HASH}]"
    width = max(60, int(fig.get_figwidth() * 19))
    lines = textwrap.wrap(txt, width)
    if len(lines) > 1:                          # make room for the extra line(s) after tight_layout
        fig.subplots_adjust(bottom=min(0.45, fig.subplotpars.bottom + 0.035 * (len(lines) - 1)))
    fig.text(0.01, 0.005, "\n".join(lines), fontsize=6.5, color=INK2, ha="left", va="bottom")
