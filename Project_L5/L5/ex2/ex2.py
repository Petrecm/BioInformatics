# ex2.py
# For 10 viral genomes: take random reads, rebuild them,
# measure assembly time and GC%, and plot Assembly time vs GC%.

import os, glob, random, time
import matplotlib.pyplot as plt

random.seed(1)  # keep results repeatable

# read a FASTA file
def load_fasta(path):
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))

# find overlap between two reads (suffix of a and prefix of b)
def overlap(a, b, k=10):
    m = min(len(a), len(b))
    for t in range(m, k - 1, -1):
        if a[-t:] == b[:t]:
            return t
    return 0

# simple greedy assembly
def assemble(reads):
    seq = reads[0]
    for r in reads[1:]:
        ov = overlap(seq, r, 10)
        if ov >= 10:
            seq += r[ov:]
        else:
            seq += r
    return seq

# compute GC percentage
def gc_percent(seq):
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)

# main part
folder = "viruses"
files = sorted(glob.glob(os.path.join(folder, "*.fasta")))

times_ms = []
gcs = []
labels = []

for i, f in enumerate(files, 1):
    seq = load_fasta(f)

    # take 2000 random reads of ~100–150 bases
    reads = [seq[random.randint(0, len(seq) - 120):][:random.randint(100, 150)]
             for _ in range(2000)]

    # measure assembly time in milliseconds
    start = time.perf_counter()
    _ = assemble(reads)
    end = time.perf_counter()
    t_ms = (end - start) * 1000.0

    # compute GC%
    gc = gc_percent(seq)

    times_ms.append(t_ms)
    gcs.append(gc)
    labels.append(f"V{i}")

    print(f"{os.path.basename(f)}: len={len(seq)} bp | GC={gc:.2f}% | time={t_ms:.1f} ms")

# plot GC% (X) vs time (Y)
plt.figure(figsize=(7,5))
plt.scatter(gcs, times_ms)
for x, y, lab in zip(gcs, times_ms, labels):
    plt.text(x + 0.2, y, lab, fontsize=8)

plt.xlabel("GC content (%)")
plt.ylabel("Assembly time (ms)")
plt.title("Assembly time vs GC% (10 viral genomes)")
plt.grid(True)
plt.tight_layout()
plt.savefig("Screenshot1_Ex2.jpg", dpi=150)
# plt.show()

print("\nSaved: Screenshot1_Ex2.jpg")