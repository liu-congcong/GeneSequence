# GeneSequence

Generate processed RNA sequences, CDS sequences and translated protein sequences of all transcripts.

* >X_**transcript**: processed RNA sequence of transcript X

* >X_**cds**: CDS sequence of transcript X

* >X_**protein**: protein sequence of transcript X

## Installation

```shell
git clone https://github.com/liu-congcong/GeneSequence.git
cd GeneSequence
gcc *.c -o gene_sequence -lm
```

## Usage

```shell
gene_sequence -fasta FASTA -gff GFF > OUTPUT.FASTA
```
