# GeneSequence

Generate processed RNA sequences, CDS sequences and translated protein sequences of all transcripts.

## Installation

```shell
git clone https://github.com/liu-congcong/GeneSequence.git
cd GeneSequence
make && make clean
```

## Usage

* Transcript sequences

```shell
GeneSequence -fasta FASTA -gff GFF -type transcript > TRANSCRIPT.FASTA
```

* CDS sequences

```shell
GeneSequence -fasta FASTA -gff GFF -type cds > CDS.FASTA
```

* Protein sequences

```shell
GeneSequence -fasta FASTA -gff GFF -type protein > PROTEIN.FASTA
```
