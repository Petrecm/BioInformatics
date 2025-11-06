import math, random, sys
from pathlib import Path
import matplotlib.pyplot as plt

def load_fasta(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing '{path.name}' in parent folder.")
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"): 
                continue
            seq.append("".join(ch for ch in line.strip().upper() if ch in "ACGT"))
    s = "".join(seq)
    if not s:
        raise ValueError("FASTA has no A/C/G/T bases.")
    return s

def sample_fragments(seq: str, n=10, min_bp=100, max_bp=3000, seed=42):
    random.seed(seed)
    L = len(seq)
    out = []
    for i in range(n):
        start = random.randint(0, max(0, L - min_bp))
        max_here = min(max_bp, L - start)
        length = random.randint(min_bp, max_here)
        out.append({
            "i": i + 1,
            "start": start,
            "end": start + length,
            "bp": length,
            "seq": seq[start:start+length]
        })
    return out

def dist(bp, a=100.0, b=20.0):
    return max(5.0, a - b * math.log10(bp))

def plot_gel(frags, out="gel.png"):
    bands = [{"i":f["i"], "bp":f["bp"], "y":dist(f["bp"])} for f in frags]
    bands.sort(key=lambda x: x["y"])

    min_gap = 6.0
    label_y = []
    for b in bands:
        y = b["y"] if not label_y else max(b["y"], label_y[-1] + min_gap)
        label_y.append(y)
    for b, ylab in zip(bands, label_y):
        b["y_label"] = ylab

    fig, ax = plt.subplots(figsize=(4.2, 9), dpi=220)
    lane_left, lane_right = 0.46, 0.54
    ax.fill_between([lane_left, lane_right], [0, 0], [110, 110])

    for b in bands:
        ax.hlines(b["y"], lane_left + 0.015, lane_right - 0.015, linewidth=4)
    for b in bands:
        ax.text(lane_right + 0.07, b["y_label"], f"{b['bp']} bp",
                va="center", ha="left", fontsize=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 110)
    ax.set_xlabel("Lane")
    ax.set_ylabel("Migration distance (arbitrary units)")
    ax.set_title("Simulated Gel Electrophoresis – 10 DNA Fragments")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] gel saved: {out}")

def main():
    # go exactly 2 folders back: ex1 → L6 → Project_L6
    fasta_path = Path(__file__).resolve().parents[2] / "gene.fasta"

    try:
        seq = load_fasta(fasta_path)
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

    frags = sample_fragments(seq)
    print(f"[info] Loaded gene.fasta ({len(seq)} bp)")
    print("(index) start-end  length_bp")
    for f in frags:
        print(f"  {f['i']:02d}    {f['start']:4d}-{f['end']:4d}   {f['bp']:4d} bp")

    plot_gel(frags, out="gel.png")

if __name__ == "__main__":
    main()
