from itertools import product
from collections import Counter
import csv, os

seq = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".upper()
letters = ['A','C','G','T']

def make_all(k):
    res = []
    for p in product(letters, repeat=k):
        res.append(''.join(p))
    return res

def count_kmers(s, k):
    c = Counter()
    for i in range(len(s)-k+1):
        kmer = s[i:i+k]
        c[kmer] += 1
    return c

def table(s, k):
    combs = make_all(k)
    counts = count_kmers(s,k)
    total = len(s)-k+1
    data = []
    for x in sorted(combs):
        n = counts.get(x,0)
        p = round(n/total*100,2)
        data.append((x,n,total,p))
    return data


print("sequence length", len(seq))
for k in (2,3):
    t = table(seq,k)
    print("k =",k)
    for r in t:
        print(r)
    print()

if not os.path.exists("out"):
    os.mkdir("out")

for k in (2,3):
    fn = "out/" + ("di" if k==2 else "tri") + "nucleoide.csv"
    with open(fn,"w",newline="") as f:
        w = csv.writer(f)
        w.writerow(["kmer","count","total","percent"])
        for r in table(seq,k):
            w.writerow(r)
    print("saved",fn)
