import sys
import math
import matplotlib.pyplot as plt

def read_fasta(sequence):
    seq = []
    with open(sequence) as f:
        for line in f:
            if not line.startswith('>'):
                seq.append(line.strip())
    return ''.join(seq)

def tm(sequence):
    a = sequence.count('A')
    t = sequence.count('T')
    g = sequence.count('G')
    c = sequence.count('C')
    return 4*(g+c)+2*(a+t)

def tm2(sequence):
    sequence = sequence.upper()
    a = sequence.count('A')
    t = sequence.count('T')
    g = sequence.count('G')
    c = sequence.count('C')
    gc_percent = ((g+c)/len(sequence))*100
    Nat = 0.01
    return 81.5+16.6*math.log10(Nat)+0.41*gc_percent-600/len(sequence)

def sliding_window_tm(sequence, window_size=8):
    results = []
    for i in range(len(sequence)-window_size+1):
        window_seq = sequence[i:i+window_size]
        temp1 = tm(window_seq)
        temp2 = tm2(window_seq)
        results.append((i+1, window_seq, temp1, temp2))
    return results

def main():
    if len(sys.argv)!=2:
        print("Usage: python lab3ex2.py <input.fasta>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)
    window_size = 8
    results = sliding_window_tm(sequence, window_size)
    
    print("Pos\tSequence\tTm1(C)\tTm2(C)")
    for pos, seq, temp1, temp2 in results:
        print(f"{pos}\t{seq}\t{temp1}\t{temp2:.2f}")
    
    # prepare data for plotting
    positions = [pos for pos, _, _, _ in results]
    tm1_values = [temp1 for _, _, temp1, _ in results]
    tm2_values = [temp2 for _, _, _, temp2 in results]
    
    # plot
    plt.figure(figsize=(10, 5))
    plt.plot(positions, tm1_values, marker='o', color='royalblue', label='Tm1 (simple formula)')
    plt.plot(positions, tm2_values, marker='s', color='darkorange', label='Tm2 (Wallace rule)')
    plt.title("Melting Temperature (Tm) Across DNA Sequence")
    plt.xlabel("Position (start of 8-mer window)")
    plt.ylabel("Tm (C)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    # save and show
    plt.savefig("tm_plot.png", dpi=300)
    print("\n✅ Plot saved as 'tm_plot.png'")
    plt.show()

if __name__ == "__main__":
    main()
