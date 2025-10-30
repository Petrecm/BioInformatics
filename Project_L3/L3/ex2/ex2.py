#!/usr/bin/env python3
"""
Windowed nucleotide composition & melting temperature (Tm) from FASTA.

What it does
------------
- Reads a (multi-)FASTA.
- Slides a fixed window along each sequence (size/step set in constants).
- For each window:
  * computes %A, %C, %G, %T (non-ACGT characters are ignored),
  * Tm via Wallace rule:        Tm = 2*(A+T) + 4*(G+C)
  * Tm via salt-adjusted rule:  Tm = 81.5 + 16.6*log10([Na+]) + 0.41*%GC - 600/len
- Saves a CSV with all values and a PNG figure (both Tm curves on the same axes).
- One CSV + one PNG per sequence in the FASTA.

Notes
-----
- If a window has no A/C/G/T, values are recorded as NaN.
- The x-axis uses the 1-based center of each window.
"""

from __future__ import annotations

import os
import math
from collections import Counter
from typing import Iterator, Tuple, List

import pandas as pd
import matplotlib.pyplot as plt


# ========= adjustable settings (edit here) =========
FASTA_FILE = "../../dna.fasta"   # input FASTA
WIN_SIZE   = 9                # window length in bases
WIN_STEP   = 1                # step size between windows
OUTPUT_DIR = "results"        # where CSV/PNG go
OUT_PREFIX = None             # None → derive from FASTA header; or set a custom string
NA_CONC_M  = 0.001            # sodium concentration [M] for the complex Tm
# ===================================================


# ---------- FASTA parsing ----------

def fasta_reader(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) for each record in a FASTA file."""
    header: str | None = None
    chunk: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if header is not None:
                    yield header, "".join(chunk).upper()
                header = s[1:].strip()
                chunk = []
            else:
                chunk.append(s)
    if header is not None:
        yield header, "".join(chunk).upper()


# ---------- core calculations ----------

def tm_wallace(seq: str) -> float:
    """Wallace rule on A/C/G/T only; NaN if none."""
    s = "".join(b for b in seq.upper() if b in "ACGT")
    if not s:
        return float("nan")
    a, c, g, t = s.count("A"), s.count("C"), s.count("G"), s.count("T")
    return 2 * (a + t) + 4 * (g + c)


def tm_salt_adjusted(seq: str, na_molar: float) -> float:
    """Salt-adjusted empirical Tm."""
    if na_molar <= 0:
        raise ValueError("[Na+] must be positive (mol/L)")
    s = "".join(b for b in seq.upper() if b in "ACGT")
    if not s:
        return float("nan")
    n = len(s)
    gc_pct = 100.0 * (s.count("G") + s.count("C")) / n
    return 81.5 + 16.6 * math.log10(na_molar) + 0.41 * gc_pct - (600.0 / n)


def slide_compute(
    seq: str,
    win: int,
    step: int,
    na_molar: float
) -> pd.DataFrame:
    """
    Return a DataFrame with columns:
    position (center, 1-based), A, C, G, T, Tm_simple, Tm_complex
    """
    valid = set("ACGT")
    n = len(seq)

    records = {
        "position": [],
        "A": [], "C": [], "G": [], "T": [],
        "Tm_simple": [], "Tm_complex": []
    }

    for start in range(0, n - win + 1, step):
        frag = seq[start:start + win]
        clean = "".join(b for b in frag if b in valid)
        counts = Counter(clean)
        denom = sum(counts.values())

        if denom == 0:
            a = c = g = t = math.nan
            tms = tmc = math.nan
        else:
            a = 100.0 * counts.get("A", 0) / denom
            c = 100.0 * counts.get("C", 0) / denom
            g = 100.0 * counts.get("G", 0) / denom
            t = 100.0 * counts.get("T", 0) / denom
            tms = tm_wallace(clean)
            tmc = tm_salt_adjusted(clean, na_molar)

        center_1based = start + (win // 2) + 1

        records["position"].append(center_1based)
        records["A"].append(a); records["C"].append(c)
        records["G"].append(g); records["T"].append(t)
        records["Tm_simple"].append(tms); records["Tm_complex"].append(tmc)

    return pd.DataFrame(records)


# ---------- plotting ----------

def plot_both_tm(df: pd.DataFrame, title: str, out_file: str) -> None:
    """Draw Tm (Wallace) and Tm (salt-adjusted) on one figure and save PNG."""
    plt.figure(figsize=(10, 5), dpi=120)
    plt.plot(df["position"], df["Tm_simple"], "--", label="Tm (Wallace)")
    plt.plot(df["position"], df["Tm_complex"], label="Tm (salt-adjusted)")
    plt.xlabel("Sequence position (bp)")
    plt.ylabel("Melting temperature (°C)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


# ---------- utils ----------

def safe_name(s: str) -> str:
    """Make a filesystem-friendly stem."""
    s = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in s)
    return s[:80] or "sequence"


# ---------- main workflow ----------

def main() -> None:
    # quick checks on constants
    if WIN_SIZE <= 0:
        raise SystemExit("WIN_SIZE must be > 0")
    if WIN_STEP <= 0:
        raise SystemExit("WIN_STEP must be > 0")
    if NA_CONC_M <= 0:
        raise SystemExit("NA_CONC_M must be > 0")
    if not os.path.exists(FASTA_FILE):
        raise SystemExit(f"FASTA not found: {FASTA_FILE}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for idx, (header, seq) in enumerate(fasta_reader(FASTA_FILE), start=1):
        if len(seq) < WIN_SIZE:
            print(f"[skip] seq {idx} '{header}': length {len(seq)} < window {WIN_SIZE}")
            continue

        df = slide_compute(seq, WIN_SIZE, WIN_STEP, NA_CONC_M)

        stem = OUT_PREFIX or safe_name(header or f"seq{idx}")
        csv_path = os.path.join(OUTPUT_DIR, f"{stem}.csv")
        png_path = os.path.join(OUTPUT_DIR, f"{stem}.png")

        df.to_csv(csv_path, index=False)
        ttl = f"Tm (Wallace & salt-adjusted) — W={WIN_SIZE}, step={WIN_STEP}, [Na+]={NA_CONC_M} M\n{header}"
        plot_both_tm(df, ttl, png_path)

        print(f"[ok] wrote: {csv_path}  |  {png_path}")


if __name__ == "__main__":
    main()
