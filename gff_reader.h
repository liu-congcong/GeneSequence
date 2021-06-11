#ifndef __GFF_READER_H__

#define __GFF_READER_H__

#define LINE 10240

#define FILE_NAME 10240

typedef struct Element{ // exon or cds
    short type;
    short phase;
    size_t start;
    size_t end;
    struct Element *next;
} Element;

typedef struct Transcript{
    char *transcript;
    char *ref_name;
    char strand;
    size_t exon_number;
    size_t cds_number;
    struct Transcript *next;
    Element *element;
} Transcript;

int read_gff_file(char *, Transcript ***, size_t);

int free_gff_hash(Transcript **, size_t);

#endif