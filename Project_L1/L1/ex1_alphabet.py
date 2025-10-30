#Find the alphabet unique symbols of a sequence S.S

def alphabet(sequence: str) -> set:
    """
    Return the set of unique (uppercase) symbols in the sequence.
    """
    return set(sequence.upper())

if __name__ == "__main__":
    S = "ATTTCGCCGATA"
    alph = alphabet(S)

    print("Sequence S:", S)
    print("Alphabet   :", "{" + ", ".join(sorted(alph)) + "}")
