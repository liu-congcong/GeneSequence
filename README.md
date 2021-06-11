# GeneSequence

Input FASTA and GFF3 files, and output a FASTA file containing all transcripts, cds sequences and translated proteins.

## Installation

Clone and compile as follows:

```shell
git clone https://github.com/liu-congcong/GeneSequence.git
cd GeneSequence
gcc *.c -o gene_sequence -lm
```

## Usage

```shell
gene_sequence -fasta FASTA -gff GFF > OUTPUT.FASTA
```
