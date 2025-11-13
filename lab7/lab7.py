from collections import defaultdict
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
        if not seq or not (min_len <= len(seq) <= max_len):
            raise RuntimeError('Fetched sequence outside desired length or empty.')
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

def find_repetitive_sequences(dna_sequence, min_length=3, max_length=6, min_repetitions=2):
    repetitive_sequences = defaultdict(list)
    
    for pattern_length in range(min_length, max_length + 1):
        for i in range(len(dna_sequence) - pattern_length + 1):
            pattern = dna_sequence[i:i + pattern_length]
            
            positions = []
            for j in range(len(dna_sequence) - pattern_length + 1):
                if dna_sequence[j:j + pattern_length] == pattern:
                    positions.append(j)
            
            if len(positions) >= min_repetitions:
                if pattern not in repetitive_sequences:
                    repetitive_sequences[pattern] = positions
    
    return repetitive_sequences

def display_results(dna_sequence, repetitive_sequences):
    print(f"\nDNA Sequence Length: {len(dna_sequence)} nucleotides\n")
    print(f"Found {len(repetitive_sequences)} repetitive sequences:\n")
    
    sorted_patterns = sorted(
        repetitive_sequences.items(),
        key=lambda x: (len(x[1]), len(x[0])),
        reverse=True
    )
    
    for pattern, positions in sorted_patterns[:20]:
        print(f"Pattern: {pattern}")
        print(f"Length: {len(pattern)} bp")
        print(f"Repetitions: {len(positions)} times")
        print(f"Positions: {positions[:10]}{'...' if len(positions) > 10 else ''}")
        print("-" * 50)

if __name__ == "__main__":
    print("DNA Repetitive Sequence Detector")
    print("=" * 50)
    
    result = fetch_ncbi_sequence(min_len=1000, max_len=3000)
    print(f"Using sequence: {result['accession']}")
    print(f"Sequence length: {len(result['sequence'])} nucleotides\n")
    
    DNA_SEQUENCE = result['sequence']
    
    if DNA_SEQUENCE is None:
        print("Failed to fetch DNA sequence. Exiting.")
        exit(1)
    
    results = find_repetitive_sequences(DNA_SEQUENCE, min_length=3, max_length=6, min_repetitions=2)
    
    display_results(DNA_SEQUENCE, results)
    
    print(f"\nTotal unique repetitive patterns found: {len(results)}")
