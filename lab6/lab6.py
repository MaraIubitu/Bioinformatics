import random
import math
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
    except Exception:
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

def sample_fragments(seq, n=10, min_len=100, max_len=300, seed=None):
    if seed is not None:
        random.seed(seed)
    L = len(seq)
    frags = []
    for i in range(n):
        frag_len = random.randint(min_len, max_len)
        if frag_len >= L:
            frag_len = L - 1
        start = random.randint(0, L - frag_len)
        frag_seq = seq[start:start + frag_len]
        frags.append({
            'id': f'F{i+1}',
            'start': start,
            'length': frag_len,
            'sequence': frag_seq
        })
    return frags

def migrate_position(bp, min_bp, max_bp, height):
    bp = max(1, bp)
    log_min = math.log10(min_bp)
    log_max = math.log10(max_bp)
    val = (log_max - math.log10(bp)) / (log_max - log_min)
    y = int(round(val * (height - 1)))
    return max(0, min(height - 1, y))

def ascii_gel(fragment_lengths, ladder_bp=None, gel_height=40, lane_width=7, smear=False, seed=None):
    if seed is not None:
        random.seed(seed)
    if ladder_bp is None:
        ladder_bp = [100, 200, 300, 400, 500, 700, 1000, 1500, 2000, 2500, 3000]
    all_bp = ladder_bp + fragment_lengths
    min_bp = max(50, min(all_bp))
    max_bp = max(all_bp)
    lanes = 1 + len(fragment_lengths)
    pad = 2
    width = pad * 2 + lanes * lane_width
    height = gel_height
    grid = [[' ' for _ in range(width)] for _ in range(height)]
    for lane in range(lanes):
        cx = pad + lane * lane_width + lane_width // 2
        for dx in range(-lane_width // 2 + 1, lane_width // 2):
            x = cx + dx
            if 0 <= x < width:
                grid[0][x] = '-'
    def draw_band(bp, lane_idx, char='='):
        cx = pad + lane_idx * lane_width + lane_width // 2
        y = migrate_position(bp, min_bp, max_bp, height)
        smear_span = 0
        if smear:
            smear_span = random.randint(0, 1)
        for dy in range(-smear_span, smear_span + 1):
            yy = y + dy
            if 0 <= yy < height:
                for dx in range(-2, 3):
                    x = cx + dx
                    if 0 <= x < width:
                        grid[yy][x] = char
    for bp in ladder_bp:
        if min_bp <= bp <= max_bp:
            draw_band(bp, 0, char='#')
    for i, bp in enumerate(fragment_lengths, start=1):
        draw_band(bp, i, char='=')
    lines = [''.join(row) for row in grid]
    return '\n'.join(lines)

def main():
    seq_info = fetch_ncbi_sequence(min_len=1000, max_len=3000, seed=42)
    seq = seq_info['sequence']
    fragments = sample_fragments(seq, n=10, min_len=100, max_len=300, seed=42)
    samples = fragments
    frag_lengths = [f['length'] for f in samples]
    gel = ascii_gel(frag_lengths, gel_height=40, lane_width=7, smear=False, seed=42)
    print(f"Sequence: {seq_info['accession']}")
    print(f"Sequence length: {len(seq)} bp")
    print("Sample fragments (lane -> length bp):")
    for i, f in enumerate(samples, start=1):
        print(f"  Lane {i}: {f['length']} bp (start {f['start']})")
    print("\nASCII gel (Lane 0 = Ladder, Lanes 1-10 = Samples):\n")
    print(gel)
    print("\nLegend:")
    print("  = : sample fragment bands")
    print("  # : ladder bands")
    print("  - : wells at top")

if __name__ == '__main__':
    main()