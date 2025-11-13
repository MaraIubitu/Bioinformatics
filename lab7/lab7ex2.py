from collections import Counter
import matplotlib.pyplot as plt
import random
import json
from urllib import request, parse
from contextlib import closing

def fetch_ncbi_sequence(min_len=1000, max_len=3000, retmax=100, api_key=None, seed=None):
    if seed is not None:
        random.seed(seed)
    try:
        term = f'{min_len}:{max_len}[SLEN] AND biomol_genomic[PROP] NOT mitochondrial[Title] NOT chloroplast[Title]'
        params = {
            'db': 'nuccore',
            'term': term,
            'retmode': 'json',
            'retmax': str(retmax),
        }
        if api_key:
            params['api_key'] = api_key
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + parse.urlencode(params)
        with closing(request.urlopen(url, timeout=10)) as resp:
            data = json.load(resp)
        ids = data.get('esearchresult', {}).get('idlist', [])
        if not ids:
            raise RuntimeError('No IDs returned from NCBI search.')
        picked = random.choice(ids)
        params = {
            'db': 'nuccore',
            'id': picked,
            'rettype': 'fasta',
            'retmode': 'text',
        }
        if api_key:
            params['api_key'] = api_key
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + parse.urlencode(params)
        with closing(request.urlopen(url, timeout=10)) as resp:
            fasta = resp.read().decode('utf-8', errors='ignore')
        header, seq = parse_fasta(fasta)
        if not seq:
            raise RuntimeError('Fetched sequence is empty.')
        return {'accession': header, 'sequence': seq}
    except Exception as e:
        print(f"Error fetching from NCBI: {e}")
        print("Falling back to synthetic random sequence...")
        L = random.randint(min_len, max_len)
        seq = random_dna(L)
        return {'accession': 'SYNTHETIC_RANDOM_SEQ', 'sequence': seq}

def parse_fasta(fasta_text):
    lines = [l.strip() for l in fasta_text.splitlines() if l.strip()]
    if not lines or not lines[0].startswith('>'):
        return ('UNKNOWN', '')
    header = lines[0][1:]
    seq = ''.join(lines[1:]).upper()
    seq = ''.join(c for c in seq if c in 'ACGTN')
    return (header, seq)

def random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def download_influenza_genomes(count=10):
    """Download genome sequences from NCBI"""
    sequences = []
    print(f"Fetching {count} genome sequences...")
    
    for i in range(count):
        result = fetch_ncbi_sequence(min_len=1000, max_len=15000, retmax=100, seed=None)
        sequences.append(result['sequence'])
        print(f"  Downloaded genome {i+1}/{count}: {result['accession'][:50]}...")
    
    return sequences

def find_most_frequent_repeat(sequence, min_length=3, max_length=10):
    """Find the most frequent repeated substring in a sequence"""
    repeat_counts = Counter()
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            substring = sequence[i:i+length]
            if substring.count('N') == 0:  # Skip substrings with ambiguous bases
                repeat_counts[substring] += 1
    
    # Get the most common repeat that appears more than once
    most_common = [(seq, count) for seq, count in repeat_counts.most_common(20) if count > 1]
    return most_common[0] if most_common else (None, 0)

# Download genomes
print("Downloading influenza genomes...")
genomes = download_influenza_genomes(10)

# Find most frequent repeats for each genome
results = []
for i, genome in enumerate(genomes):
    repeat, count = find_most_frequent_repeat(genome)
    results.append((i+1, repeat, count))
    print(f"Genome {i+1}: Most frequent repeat '{repeat}' appears {count} times")

# Plot the results
genome_numbers = [r[0] for r in results]
frequencies = [r[2] for r in results]

plt.figure(figsize=(10, 6))
plt.bar(genome_numbers, frequencies, color='steelblue', edgecolor='black')
plt.xlabel('Genome Number')
plt.ylabel('Frequency of Most Common Repeat')
plt.title('Most Frequent Repetitions in 10 Influenza Genomes')
plt.xticks(genome_numbers)
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.show()