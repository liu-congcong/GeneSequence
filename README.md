# GeneSequence

Input FASTA and GFF3 files, and output a FASTA file containing all transcripts, cds sequences and translated proteins.

## Installation

Download all files and compile as follows:

```shell
gcc *.c -o gene_sequence
```

## Usage

```shell
gene_sequence -fasta FASTA -gff GFF > OUTPUT.FASTA
```
