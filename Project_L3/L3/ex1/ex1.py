import math

def base_freq(sequence: str, base: str) -> float:
    """
    Returns frequency (0–1) of a specific base within a sequence.
    """
    if not sequence:
        return 0.0
    seq = sequence.upper()
    b = base.upper()
    return seq.count(b) / len(seq)


def base_count(sequence: str, base: str) -> int:
    """
    Returns the number of times a base appears in a sequence.
    """
    return sequence.upper().count(base.upper())


def melting_temp_simple(sequence: str) -> float:
    """
    Compute melting temperature using the Wallace rule:
    Tm = 2*(A+T) + 4*(G+C)
    """
    seq = sequence.upper()
    a, t, g, c = seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")
    tm = 2 * (a + t) + 4 * (g + c)
    return tm


def melting_temp_detailed(sequence: str, na_conc: float = 0.001) -> float:
    """
    More complex empirical formula for oligonucleotide melting temperature.
    """
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        raise ValueError("Sequence is empty.")
    if na_conc <= 0:
        raise ValueError("Sodium concentration must be positive (mol/L).")

    gc_percent = 100 * (seq.count("G") + seq.count("C")) / n
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / n)
    return tm


# Example usage
dna_seq = "ATTTCGCCGATA"
print(f"Simple melting temperature: {melting_temp_simple(dna_seq)}°C")
print(f"Complex melting temperature: {melting_temp_detailed(dna_seq):.2f}°C")
