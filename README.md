# GeneSequence

Generate processed RNA sequences, CDS sequences and translated protein sequences of all transcripts.

## Installation

```shell
git clone https://github.com/liu-congcong/GeneSequence.git
cd GeneSequence
gcc *.c -o gene_sequence -lm
```

## Usage

* Transcript sequences

```shell
gene_sequence -fasta FASTA -gff GFF -type transcript > TRANSCRIPT.FASTA
```

* CDS sequences

```shell
gene_sequence -fasta FASTA -gff GFF -type cds > CDS.FASTA
```

* Protein sequences

```shell
gene_sequence -fasta FASTA -gff GFF -type protein > PROTEIN.FASTA
```
