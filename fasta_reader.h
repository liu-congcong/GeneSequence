#ifndef __FASTA_READER_H__

#define __FASTA_READER_H__

#define LINE 10240

#define FILE_NAME 10240

typedef struct Sequence{
    char *id;
    char *sequence;
    unsigned long length;
    struct Sequence *next;
} Sequence;

unsigned long read_fasta_file(char *, Sequence ***, unsigned long);

int free_fasta_hash(Sequence **, unsigned long);

char *find_sequence(Sequence **, unsigned long, char *);

#endif