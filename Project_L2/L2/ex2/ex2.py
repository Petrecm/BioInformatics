seq = "ATTGTCCCAATCTGTTG"

def kmer_frequency(sequence: str, k: int):
    """
    Calculate the relative frequency (percentage) of every k-mer
    found in the given sequence.
    """
    counts = {}
    n = len(sequence)
    total = 0

    # slide window through sequence
    for i in range(n - k + 1):
        frag = sequence[i:i + k]
        counts[frag] = counts.get(frag, 0) + 1
        total += 1

    # convert to percentage strings with 4 decimals
    freq = {
        mer: f"{(cnt / total) * 100:.4}%" for mer, cnt in counts.items()
    }

    # return sorted alphabetically by k-mer
    return dict(sorted(freq.items()))

print(kmer_frequency(seq, 2))
print(kmer_frequency(seq, 3))
