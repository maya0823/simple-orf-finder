# simple-orf-finder

A small pure-Python tool that takes a DNA sequence, computes basic stats, the reverse complement, and looks for open reading frames (ORFs) that start with ATG and end with a stop codon.

It is meant as a quick way to explore string processing and a tiny bit of molecular biology without any external libraries.

## What you can learn

- How to walk a sequence in different reading frames
- Building a simple codon table and translating DNA to protein
- Cleaning and validating biological sequence data

## Learning from mistakes

- Forgetting that ORFs can start in any of the three frames (or on the reverse strand) is the most common miss.
- Stop codons only count when they are in-frame; an out-of-frame TAA does not end the ORF.

## How to run

Requires Python 3. No dependencies.

```bash
python main.py ATGAAACCCGGGTTTAAATAG
```

Or run without arguments and paste a sequence (or press Enter for a short demo).

```bash
python main.py
```

The script prints length, GC content, reverse complement, and any ORFs longer than 30 bp.

## Notes

This only searches the forward strand and uses a simple start-to-first-stop logic. Real annotation is much more involved.

MIT License
