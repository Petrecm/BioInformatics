CODON = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}
AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu","G":"Gly",
    "H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","S":"Ser",
    "T":"Thr","W":"Trp","Y":"Tyr","V":"Val","*":"Stop"
}
STOPS = {"UAA","UAG","UGA"}

def tidy(seq: str) -> str:
    s = "".join(seq.upper().split())   
    return s.replace("T", "U")          

def translate_first_orf(seq: str):
    rna = tidy(seq)
    start = rna.find("AUG")
    if start == -1:
        return "", "", "No AUG (start) found."
    aa, i = [], start
    while i + 2 < len(rna):
        codon = rna[i:i+3]
        if codon in STOPS:
            break
        aa1 = CODON.get(codon)
        if aa1 is None:
            return "", "", f"Invalid codon {codon} at {i+1}-{i+3}."
        aa.append(aa1)
        i += 3
    one = "".join(aa)
    three = "-".join(AA3[x] for x in one) if one else ""
    note = "OK (AUG→STOP)" if i + 2 < len(rna) and rna[i:i+3] in STOPS else "No in-frame STOP."
    return one, three, note

def show(label: str, seq: str):
    aa1, aa3, info = translate_first_orf(seq)
    print(f"\n== {label} ==")
    print("Sequence (5'→3'):", seq)
    print("Info:", info)
    print("AA (1-letter):", aa1 if aa1 else "(empty)")
    print("AA (3-letter):", aa3 if aa3 else "(empty)")

if __name__ == "__main__":
   
    demo_seq = "GGCUUCCGCGGCUUUGGCC AUG GCU GAA CCU UUC CCC GGC ACU UAA CCGGGCUUCGGCCG"
    show("Hardcoded Sequence", demo_seq)

    
    try:
        user = input("\nPaste your DNA/RNA (or press Enter to skip): ").strip()
        if user:
            show("Your input", user)
    except KeyboardInterrupt:
        pass




