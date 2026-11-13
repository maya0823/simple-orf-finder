#!/usr/bin/env python3
"""
Simple ORF finder and DNA utilities.
Pure Python, standard library only.
"""

import sys


CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def clean_sequence(seq: str) -> str:
    """Keep only A, C, G, T and uppercase."""
    return ''.join(c for c in seq.upper() if c in 'ACGT')


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def gc_content(seq: str) -> float:
    """Return GC content as a percentage."""
    if not seq:
        return 0.0
    gc = sum(1 for base in seq if base in 'GC')
    return 100.0 * gc / len(seq)


def translate(seq: str) -> str:
    """Translate a DNA sequence into amino acids. Stops at first stop codon."""
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def find_orfs(seq: str, min_length: int = 30) -> list:
    """
    Find open reading frames starting with ATG and ending with a stop codon.
    Returns list of (start, end, protein, frame) tuples.
    frame is 0, 1, or 2 for the forward strand.
    """
    orfs = []
    n = len(seq)
    for frame in range(3):
        i = frame
        while i + 2 < n:
            codon = seq[i:i+3]
            if codon == 'ATG':
                # Look for the next stop
                j = i + 3
                while j + 2 < n:
                    stop = seq[j:j+3]
                    if stop in STOP_CODONS:
                        orf_seq = seq[i:j+3]
                        if len(orf_seq) >= min_length:
                            protein = translate(orf_seq)
                            orfs.append((i, j+3, protein, frame))
                        break
                    j += 3
            i += 3
    return orfs


def main():
    if len(sys.argv) > 1:
        raw = ' '.join(sys.argv[1:])
    else:
        print("Enter a DNA sequence (or press Enter for a demo sequence):")
        raw = input().strip()
        if not raw:
            raw = "ATGAAACCCGGGTTTAAATAGATGCCCGGGTAACCCGGGTTTAAATAA"

    seq = clean_sequence(raw)
    if not seq:
        print("No valid DNA bases found.")
        return

    print(f"Sequence length: {len(seq)}")
    print(f"GC content: {gc_content(seq):.1f}%")
    print(f"Reverse complement: {reverse_complement(seq)}")
    print()

    orfs = find_orfs(seq)
    if not orfs:
        print("No ORFs found (try a longer sequence or lower min length).")
        return

    print(f"Found {len(orfs)} ORF(s):")
    for start, end, protein, frame in orfs:
        print(f"  Frame {frame}, positions {start}-{end}, length {end-start} bp")
        print(f"  Protein: {protein}")
        print()


if __name__ == "__main__":
    main()
