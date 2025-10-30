# ex1.py
import random

# A) load fasta file
def load_fasta(path):
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))

seq = load_fasta("sequence.fasta")
print("Loaded sequence length:", len(seq))

# B) make 2000 random reads of 100–150 bases
random.seed(1)
reads = []
for _ in range(2000):
    L = random.randint(100, 150)
    i = random.randint(0, len(seq) - L)
    reads.append(seq[i:i+L])

print("Created 2000 random samples (reads). Example:\n", reads[0])

# C) rebuild by overlapping reads with at least 10 bases
def overlap(a, b, k=10):
    m = min(len(a), len(b))
    for t in range(m, k-1, -1):
        if a[-t:] == b[:t]:
            return t
    return 0

assembled = reads[0]
for r in reads[1:]:
    ov = overlap(assembled, r, 10)
    if ov >= 10:
        assembled += r[ov:]
    else:
        assembled += r

print("Length of assembled sequence:", len(assembled))
