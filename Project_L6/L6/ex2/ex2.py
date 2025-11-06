# ex2.py — EcoRI digest comparison for 10 influenza genomes (local-only, influenza*.fasta only)

import math, sys, re
from pathlib import Path
from typing import List
import matplotlib.pyplot as plt

OUT_DIR = Path("gels")
MAX_FILES = 10

# Go 2 folders back: ex2 → L6 → Project_L6
ROOT = Path(__file__).resolve().parents[2]

def numeric_suffix(name: str) -> int:
    m = re.search(r'(\d+)$', name)
    return int(m.group(1)) if m else 0

def find_fastas() -> List[Path]:
    """Find up to 10 files matching influenza*.fa* in Project_L6 (ignore gene.fasta)."""
    files = [p for p in ROOT.glob("influenza*.fa*") if p.is_file()]
    # sort by trailing number so 1..10 is in human order
    files.sort(key=lambda p: numeric_suffix(p.stem))
    return files[:MAX_FILES]

def read_fasta_seq(path: Path) -> str:
    """Read FASTA and return a clean A/C/G/T sequence (concatenated)."""
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append("".join(ch for ch in line.strip().upper() if ch in "ACGT"))
    joined = "".join(seq)
    if not joined:
        raise ValueError(f"{path.name} has no valid DNA bases.")
    return joined

def digest_ecori(seq: str) -> List[int]:
    """EcoRI: GAATTC, cut G^AATTC (linear digest)."""
    motif = "GAATTC"
    cut_offset = 1
    cuts, i = [], 0
    while True:
        j = seq.find(motif, i)
        if j == -1:
            break
        cuts.append(j + cut_offset)
        i = j + 1
    if not cuts:
        return [len(seq)]
    frags, prev = [], 0
    for c in cuts:
        frags.append(c - prev)
        prev = c
    frags.append(len(seq) - prev)
    return frags

def migration_distance(bp, a=100.0, b=20.0):
    return max(5.0, a - b * math.log10(max(bp, 1)))

def plot_gel_single(fragment_lengths: List[int], title: str, out_path: Path):
    # prep bands
    bands = [{"bp": bp, "y": migration_distance(bp)} for bp in fragment_lengths]
    bands.sort(key=lambda x: x["y"])
    # de-overlap right labels
    min_gap = 6.0
    y_labels = []
    for b in bands:
        y = b["y"] if not y_labels else max(b["y"], y_labels[-1] + min_gap)
        y_labels.append(y)
    for b, yl in zip(bands, y_labels):
        b["yl"] = yl

    fig, ax = plt.subplots(figsize=(4.2, 9), dpi=220)
    lane_left, lane_right = 0.46, 0.54
    ax.fill_between([lane_left, lane_right], [0, 0], [110, 110])
    for b in bands:
        ax.hlines(b["y"], lane_left + 0.015, lane_right - 0.015, linewidth=4)
    for b in bands:
        ax.text(lane_right + 0.07, b["yl"], f"{b['bp']} bp", va="center", ha="left", fontsize=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 110)
    ax.set_xlabel("Lane")
    ax.set_ylabel("Migration distance (arbitrary units)")
    ax.set_title(title)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

def plot_comparison(names: List[str], counts: List[int], out_path: Path):
    fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
    ax.bar(range(len(names)), counts)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("EcoRI fragment count")
    ax.set_title("EcoRI Digest — Fragment Count per Influenza Genome")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

def main():
    OUT_DIR.mkdir(exist_ok=True)
    fasta_files = find_fastas()
    if not fasta_files:
        print("[error] No files matching 'influenza*.fasta' found in Project_L6.")
        sys.exit(1)

    print(f"[info] Using {len(fasta_files)} FASTA file(s) in Project_L6:")
    for f in fasta_files:
        print("   -", f.name)

    summary = []
    for fp in fasta_files:
        seq = read_fasta_seq(fp)
        frags = digest_ecori(seq)
        summary.append((fp.stem, frags))
        title = f"{fp.stem} — EcoRI digest (n={len(frags)})"
        out_img = OUT_DIR / f"{fp.stem}_gel.png"
        plot_gel_single(frags, title, out_img)
        print(f"[ok] {fp.stem}: {len(frags)} fragments → {out_img.name}")

    # summary + winner
    summary.sort(key=lambda x: len(x[1]), reverse=True)
    print("\n=== EcoRI Digest Summary (influenza only) ===")
    for name, frags in summary:
        print(f"{name}: {len(frags)} fragments")
    winner = summary[0]
    print(f"\nMost fragments: {winner[0]} ({len(winner[1])})")

    # comparison bar
    names = [n for n, _ in summary]
    counts = [len(f) for _, f in summary]
    plot_comparison(names, counts, OUT_DIR / "comparison_counts.png")
    print(f"[ok] Saved comparison bar chart → {OUT_DIR/'comparison_counts.png'}")

if __name__ == "__main__":
    main()
