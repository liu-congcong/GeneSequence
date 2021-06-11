#ifndef __FASTA_READER_H__

#define __FASTA_READER_H__

#define LINE 10240

#define FILE_NAME 10240

typedef struct Sequence{
    char *id;
    char *sequence;
    size_t length;
    struct Sequence *next;
} Sequence;

size_t read_fasta_file(char *, Sequence ***, size_t);

int free_fasta_hash(Sequence **, size_t);

char *find_sequence(Sequence **, size_t, char *);

#endif