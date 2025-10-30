#!/usr/bin/env python3
"""
Sliding-window Tm with thresholds and two plots per sequence.

Outputs per sequence:
  1) <prefix>.csv  — window table (A,C,G,T %, Tm_simple, Tm_complex)
  2) <prefix>_signals_with_thresholds.png
  3) <prefix>_regions_above_threshold.png
Also writes: results/threshold_area_summary.csv
"""

from __future__ import annotations
import os
import math
from typing import Iterator, Tuple, List, Dict
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======= SETTINGS (unchanged semantics) =======
FASTA_PATH = "../../dna.fasta"
WINDOW     = 9
STEP       = 1
OUTDIR     = "results"
PREFIX     = None
NA_MOLAR   = 0.001

THRESHOLDS: Dict[str, float] = {
    "Tm_simple":  25.0,
    "Tm_complex": -10.0,
}
# ==============================================


# -------- FASTA I/O --------
def fasta_iter(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) for each record in a FASTA file."""
    head: str | None = None
    buf: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(buf).upper()
                head = line[1:].strip()
                buf = []
            else:
                buf.append(line)
    if head is not None:
        yield head, "".join(buf).upper()


# -------- Tm formulas --------
def tm_wallace(seq: str) -> float:
    """Tm = 2*(A+T) + 4*(G+C); ignore non-ACGT."""
    s = "".join(b for b in seq.upper() if b in "ACGT")
    if not s:
        return float("nan")
    a, c, g, t = s.count("A"), s.count("C"), s.count("G"), s.count("T")
    return 2 * (a + t) + 4 * (g + c)


def tm_salt_adjusted(seq: str, na_molar: float) -> float:
    """Tm = 81.5 + 16.6*log10([Na+]) + 0.41*%GC - 600/len (on ACGT only)."""
    if na_molar <= 0:
        raise ValueError("[Na+] must be positive (mol/L).")
    s = "".join(b for b in seq.upper() if b in "ACGT")
    if not s:
        return float("nan")
    n = len(s)
    gc_pct = 100.0 * (s.count("G") + s.count("C")) / n
    return 81.5 + 16.6 * math.log10(na_molar) + 0.41 * gc_pct - (600.0 / n)


# -------- sliding-window metrics --------
def make_window_table(seq: str, win: int, step: int, na_molar: float) -> pd.DataFrame:
    """
    Return DataFrame with:
      win, position (center, 1-based), A,C,G,T (percent), Tm_simple, Tm_complex
    """
    valid = set("ACGT")
    n = len(seq)

    rows = {
        "win": [], "position": [],
        "A": [], "C": [], "G": [], "T": [],
        "Tm_simple": [], "Tm_complex": []
    }

    wid = 0
    for start in range(0, n - win + 1, step):
        wid += 1
        frag = seq[start:start + win]
        clean = "".join(b for b in frag if b in valid)
        cnt = Counter(clean)
        denom = sum(cnt.values())

        if denom == 0:
            a = c = g = t = math.nan
            tms = tmc = math.nan
        else:
            a = 100.0 * cnt.get("A", 0) / denom
            c = 100.0 * cnt.get("C", 0) / denom
            g = 100.0 * cnt.get("G", 0) / denom
            t = 100.0 * cnt.get("T", 0) / denom
            tms = tm_wallace(clean)
            tmc = tm_salt_adjusted(clean, na_molar)

        center = start + (win // 2) + 1  # 1-based
        rows["win"].append(wid)
        rows["position"].append(center)
        rows["A"].append(a); rows["C"].append(c); rows["G"].append(g); rows["T"].append(t)
        rows["Tm_simple"].append(tms); rows["Tm_complex"].append(tmc)

    return pd.DataFrame(rows)


# -------- plotting --------
def plot_signals_with_thresholds(df: pd.DataFrame, th: Dict[str, float], title: str, out_png: str) -> None:
    x = df["position"].to_numpy(float)
    y1 = df["Tm_simple"].to_numpy(float)
    y2 = df["Tm_complex"].to_numpy(float)

    plt.figure(figsize=(10, 5), dpi=120)
    plt.plot(x, y1, label="Tm_simple",  color="tab:blue",   linewidth=2.0)
    plt.plot(x, y2, label="Tm_complex", color="tab:orange", linewidth=2.0)
    plt.axhline(th["Tm_simple"],  color="tab:blue",   linestyle="--", linewidth=1.5,
                label=f"threshold simple = {th['Tm_simple']:g} °C")
    plt.axhline(th["Tm_complex"], color="tab:orange", linestyle="--", linewidth=1.5,
                label=f"threshold complex = {th['Tm_complex']:g} °C")
    plt.xlabel("Window center position (bp)")
    plt.ylabel("Melting temperature (°C)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def _mask_to_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Convert boolean mask to 1-based (start, end) inclusive runs of True."""
    runs: List[Tuple[int, int]] = []
    start = None
    for i, v in enumerate(mask, start=1):
        if v and start is None:
            start = i
        elif (not v) and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask)))
    return runs


def plot_regions(df: pd.DataFrame, th: Dict[str, float], out_png: str) -> None:
    wins = df["win"].to_numpy()
    y1 = df["Tm_simple"].to_numpy(float)
    y2 = df["Tm_complex"].to_numpy(float)

    m1 = y1 > th["Tm_simple"]
    m2 = y2 > th["Tm_complex"]

    r1 = [(s - 0.5, e - s + 1) for s, e in _mask_to_runs(m1)]
    r2 = [(s - 0.5, e - s + 1) for s, e in _mask_to_runs(m2)]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 4), dpi=120)

    ax1.broken_barh(r1, (0.25, 0.5), facecolors="tab:blue")
    ax1.set_ylim(0, 1); ax1.set_yticks([0.5]); ax1.set_yticklabels(["P1"])
    ax1.set_title(f"P1 - Regions Above Threshold ({th['Tm_simple']:.1f}°C)")

    ax2.broken_barh(r2, (0.25, 0.5), facecolors="tab:orange")
    ax2.set_ylim(0, 1); ax2.set_yticks([0.5]); ax2.set_yticklabels(["P2"])
    ax2.set_title(f"P2 - Regions Above Threshold ({th['Tm_complex']:.1f}°C)")
    ax2.set_xlabel("Window Number")

    nwin = int(wins.max()) if len(wins) else 0
    ax2.set_xlim(0, max(1, nwin) + 1)

    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


# -------- utilities --------
def area_above(x: np.ndarray, y: np.ndarray, thr: float) -> float:
    """Integral of (y - thr) where positive (trapezoid)."""
    return float(np.trapz(np.clip(y - thr, 0, None), x))


def clean_stem(s: str) -> str:
    s = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in s)
    return s[:80] or "sequence"


# -------- main --------
def main() -> None:
    # sanity checks
    if WINDOW <= 0 or STEP <= 0:
        raise SystemExit("WINDOW and STEP must be positive.")
    if NA_MOLAR <= 0:
        raise SystemExit("NA_MOLAR must be positive (mol/L).")
    for key in ("Tm_simple", "Tm_complex"):
        if key not in THRESHOLDS:
            raise SystemExit(f"THRESHOLDS missing '{key}'")
    if not os.path.exists(FASTA_PATH):
        raise SystemExit(f"FASTA not found: {FASTA_PATH}")

    os.makedirs(OUTDIR, exist_ok=True)
    summary: List[Dict[str, float]] = []

    for i, (hdr, seq) in enumerate(fasta_iter(FASTA_PATH), start=1):
        if len(seq) < WINDOW:
            print(f"[skip] seq {i} '{hdr}' too short ({len(seq)} < {WINDOW})")
            continue

        df = make_window_table(seq, WINDOW, STEP, NA_MOLAR)
        stem = PREFIX or clean_stem(hdr or f"seq{i}")

        csv_path = os.path.join(OUTDIR, f"{stem}.csv")
        df.to_csv(csv_path, index=False)

        fig1 = os.path.join(OUTDIR, f"{stem}_signals_with_thresholds.png")
        ttl = f"Tm signals with thresholds (W={WINDOW}, step={STEP}, [Na+]={NA_MOLAR} M)\n{hdr}"
        plot_signals_with_thresholds(df, THRESHOLDS, ttl, fig1)

        fig2 = os.path.join(OUTDIR, f"{stem}_regions_above_threshold.png")
        plot_regions(df, THRESHOLDS, fig2)

        print(f"[ok] {hdr} -> {csv_path}, {fig1}, {fig2}")

        x = df["position"].to_numpy(float)
        a1 = area_above(x, df["Tm_simple"].to_numpy(float),  THRESHOLDS["Tm_simple"])
        a2 = area_above(x, df["Tm_complex"].to_numpy(float), THRESHOLDS["Tm_complex"])
        summary.append({
            "sequence": hdr,
            "window": WINDOW,
            "step": STEP,
            "na_molar": NA_MOLAR,
            "threshold_simple": THRESHOLDS["Tm_simple"],
            "threshold_complex": THRESHOLDS["Tm_complex"],
            "area_above_simple": a1,
            "area_above_complex": a2,
        })

    if summary:
        out_summary = os.path.join(OUTDIR, "threshold_area_summary.csv")
        pd.DataFrame(summary).to_csv(out_summary, index=False)
        print(f"[ok] wrote {out_summary}")


if __name__ == "__main__":
    main()
