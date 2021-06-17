#ifndef __GFF_READER_H__

#define __GFF_READER_H__

#define LINE 10240

#define FILE_NAME 10240

typedef struct Element{ // exon or cds
    char type;
    short phase;
    unsigned long start;
    unsigned long end;
    struct Element *next;
} Element;

typedef struct Transcript{
    char *transcript;
    char *ref_name;
    char strand;
    unsigned long exon_number;
    unsigned long cds_number;
    struct Transcript *next;
    Element *element;
} Transcript;

int read_gff_file(char *, Transcript ***, unsigned long);

int free_gff_hash(Transcript **, unsigned long);

#endif