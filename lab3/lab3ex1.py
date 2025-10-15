import math

def compute_tm1(dna):
    dna = dna.upper()
    a = dna.count('A')
    t = dna.count('T')
    g = dna.count('G')
    c = dna.count('C')
    tm1 = 4*(g+c)+2*(a+t)
    return tm1

def compute_tm2(dna):
    dna = dna.upper()
    a = dna.count('A')
    t = dna.count('T')
    g = dna.count('G')
    c = dna.count('C')
    gc_percent = ((g + c) / len(dna)) * 100
    Nat = 0.01
    tm2=81.5+16.6*(math.log10(Nat))+.41*gc_percent-600/len(dna)
    return tm2

if __name__ == "__main__":
    seq = input("Enter DNA sequence: ").strip()
    print("tm1=", compute_tm1(seq), " C")
    print("tm2=", compute_tm2(seq), " C")
