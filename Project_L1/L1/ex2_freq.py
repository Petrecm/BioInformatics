#Compute the relative frequencies of the symbols in sequence S.

from collections import Counter

def relative_frequencies(sequence: str) -> dict:
    """
    Return a dict {symbol: freq} where freq is in [0,1].
    """
    s = sequence.upper()
    n = len(s)
    if n == 0:
        return {}
    counts = Counter(s)
    return {k: counts[k] / n for k in sorted(counts)}

if __name__ == "__main__":
    S = "ATTTCGCCGATA"
    freqs = relative_frequencies(S)

    print("Sequence S:", S)
    print("Relative frequencies:")
    for sym, f in freqs.items():
        print(f"  {sym}: {f:.6f}")
