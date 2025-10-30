#!/usr/bin/env python3
# Pure-Python codon analysis (no Biopython). Requires: matplotlib

import argparse
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional

# --- matplotlib (graceful if missing) ---
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_OK = True
except Exception:
    MATPLOTLIB_OK = False

# --- RNA codon -> 1-letter AA (Stop="*") ---
RNA2AA = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L","UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*","UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L","CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M","ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K","AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V","GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E","GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}
AA1_TO_AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu","G":"Gly",
    "H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","S":"Ser",
    "T":"Thr","W":"Trp","Y":"Tyr","V":"Val","*":"Stop"
}

@dataclass
class GenomeStats:
    name: str
    length: int
    gc_percent: float
    frame: int
    codons_counted: int
    possible_codons: int
    coverage_percent: float

# --- FASTA ---
def read_fasta_concat(path: str) -> str:
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)

# --- normalize ---
def clean_dna(s: str) -> str:
    s = s.upper()
    return "".join(b for b in s if b in "ACGT")

def dna_to_rna(s: str) -> str:
    return s.replace("T","U")

# --- counts ---
def count_codons_rna_frame(rna: str, frame: int=0) -> Counter:
    c = Counter()
    valid = set("ACGU")
    start = frame % 3
    for i in range(start, len(rna)-2, 3):
        cod = rna[i:i+3]
        if all(b in valid for b in cod):
            c[cod] += 1
    return c

def aa_counts(codon_counts: Counter, include_stop: bool=False) -> Counter:
    out = Counter()
    for cod, n in codon_counts.items():
        aa = RNA2AA.get(cod, "?")
        if aa == "*" and not include_stop:
            continue
        out[aa] += n
    return out

# --- stats ---
def gc_percent(dna: str) -> float:
    if not dna: return 0.0
    gc = sum(1 for b in dna if b in "GC")
    return 100.0 * gc / len(dna)

def build_stats(name: str, dna: str, frame: int, codon_counts: Counter) -> GenomeStats:
    rna = dna_to_rna(dna)
    counted = sum(codon_counts.values())
    possible = max((len(rna) - (frame % 3)) // 3, 0)
    cov = 100.0 * counted / possible if possible else 0.0
    return GenomeStats(name, len(dna), gc_percent(dna), frame, counted, possible, cov)

# --- compare ---
def common_codons_norm_top(c1: Counter, c2: Counter, k: int=10) -> List[Tuple]:
    n1, n2 = sum(c1.values()), sum(c2.values())
    if n1==0 or n2==0: return []
    rows = []
    for cod in (set(c1) & set(c2)):
        f1, f2 = c1[cod]/n1, c2[cod]/n2
        score = (f1*f2)**0.5
        rows.append((cod, c1[cod], c2[cod], f1, f2, score))
    rows.sort(key=lambda x: x[-1], reverse=True)
    return rows[:k]

# --- helpers ---
def fmt_top(counter: Counter, k: int=3) -> str:
    return ", ".join(f"{x}:{n}" for x,n in counter.most_common(k))

def codon_label_with_aa(codon: str) -> str:
    aa1 = RNA2AA.get(codon, "?")
    aa3 = AA1_TO_AA3.get(aa1, "?")
    return f"{codon} ({aa3})"

# --- plotting ---
def annotate_counts(ax, bars):
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width()/2, h, f"{int(h)}", ha="center", va="bottom", fontsize=9)

def plot_top(counter: Counter, title: str, png_path: Optional[Path], top_k: int=10, show: bool=False):
    if not MATPLOTLIB_OK:
        print(f"[plot] matplotlib missing, skip: {title}")
        return
    top = counter.most_common(top_k)
    if not top:
        print(f"[plot] no data for {title}")
        return
    labels_raw, counts = zip(*top)
    labels = [codon_label_with_aa(c) for c in labels_raw]
    plt.figure(figsize=(12,6))
    bars = plt.bar(range(len(labels)), counts)
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.title(title); plt.xlabel("Codon (Amino acid)"); plt.ylabel("Count")
    annotate_counts(plt.gca(), bars); plt.tight_layout()
    if png_path:
        plt.savefig(png_path, dpi=150); print(f"[plot] saved {png_path}")
    if show: plt.show()
    plt.close()

def plot_common_comparison(c1: Counter, c2: Counter, title: str, png_path: Optional[Path], top_m: int=10, show: bool=False):
    if not MATPLOTLIB_OK:
        print(f"[plot] matplotlib missing, skip: {title}")
        return
    rows = common_codons_norm_top(c1, c2, k=top_m)
    if not rows:
        print("[plot] no common codons to compare")
        return
    codons = [r[0] for r in rows]
    f1 = [r[3]*100 for r in rows]  # %
    f2 = [r[4]*100 for r in rows]
    labels = [codon_label_with_aa(c) for c in codons]
    x = list(range(len(labels))); w = 0.4
    plt.figure(figsize=(12,6))
    b1 = plt.bar([xi - w/2 for xi in x], f1, width=w, label="SARS-CoV-2")
    b2 = plt.bar([xi + w/2 for xi in x], f2, width=w, label="Influenza")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.ylabel("Normalized frequency (%)"); plt.title(title); plt.legend()
    annotate_counts(plt.gca(), b1); annotate_counts(plt.gca(), b2)
    plt.tight_layout()
    if png_path:
        plt.savefig(png_path, dpi=150); print(f"[plot] saved {png_path}")
    if show: plt.show()
    plt.close()

# --- CLI ---
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Codon frequency: SARS-CoV-2 vs Influenza (no Biopython)")
    p.add_argument("--covid", default=r"C:\FILS\BIOINF\Project_L4\coronavirus.fasta")
    p.add_argument("--flu",   default=r"C:\FILS\BIOINF\Project_L4\influenza.fasta")
    p.add_argument("--frame", type=int, default=0, choices=[0,1,2])
    p.add_argument("--topk",  type=int, default=10)
    p.add_argument("--outdir", default=".")
    p.add_argument("--show", action="store_true")
    return p.parse_args()

def main() -> int:
    a = parse_args()
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    covid_p = Path(a.covid); flu_p = Path(a.flu)
    if not covid_p.exists(): print(f"[err] missing {covid_p}"); return 2
    if not flu_p.exists():   print(f"[err] missing {flu_p}");   return 2

    covid_dna = clean_dna(read_fasta_concat(str(covid_p)))
    flu_dna   = clean_dna(read_fasta_concat(str(flu_p)))
    covid_rna = dna_to_rna(covid_dna)
    flu_rna   = dna_to_rna(flu_dna)

    covid_cod = count_codons_rna_frame(covid_rna, a.frame)
    flu_cod   = count_codons_rna_frame(flu_rna, a.frame)

    sars_stats = build_stats("SARS-CoV-2", covid_dna, a.frame, covid_cod)
    flu_stats  = build_stats("Influenza",  flu_dna,   a.frame, flu_cod)

    # (a) & (b): top 10 charts
    plot_top(covid_cod, f"Top {a.topk} Codons — SARS-CoV-2", outdir/"covid_top.png", a.topk, a.show)
    plot_top(flu_cod,   f"Top {a.topk} Codons — Influenza",  outdir/"influenza_top.png", a.topk, a.show)

    # (c): comparison chart (normalized % side-by-side)
    print("\n(c) Common codons (geom. mean of normalized freq):")
    rows = common_codons_norm_top(covid_cod, flu_cod, k=10)
    if not rows:
        print("  (none)")
    else:
        for cod, c1, c2, f1, f2, sc in rows:
            print(f"  {cod}: COVID {c1} ({f1:.3%}), Flu {c2} ({f2:.3%}), score={sc:.6f}")
    plot_common_comparison(
        covid_cod, flu_cod,
        "Common top codons (normalized frequency, %)",
        outdir/"comparison_common_top.png",
        top_m=10, show=a.show
    )

    # (d): top 3 amino acids
    covid_aa = aa_counts(covid_cod, include_stop=False)
    flu_aa   = aa_counts(flu_cod,   include_stop=False)
    print("\n(d) Top 3 amino acids (SARS-CoV-2):", [(AA1_TO_AA3[k], n) for k,n in covid_aa.most_common(3)])
    print("(d) Top 3 amino acids (Influenza):",   [(AA1_TO_AA3[k], n) for k,n in flu_aa.most_common(3)])

    # summary
    print("\n--- SUMMARY ---")
    for st in (sars_stats, flu_stats):
        print(f"{st.name}: len={st.length:,} bp | GC%={st.gc_percent:.2f} | frame={st.frame} "
              f"| codons={st.codons_counted:,}/{st.possible_codons:,} ({st.coverage_percent:.2f}% coverage)")
    print(f"Top {a.topk} SARS-CoV-2 codons: {fmt_top(covid_cod, a.topk)}")
    print(f"Top {a.topk} Influenza  codons: {fmt_top(flu_cod,  a.topk)}")
    print("Saved PNGs:", outdir/"covid_top.png", "|", outdir/"influenza_top.png", "|", outdir/"comparison_common_top.png")
    return 0

if __name__ == "__main__":
    sys.exit(main())
