#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "gff_reader.h"
#include "hash.h"

static int create_hash(Transcript ***hash, unsigned long hash_size)
{
    *hash = malloc(sizeof(Transcript *) * hash_size);
    memset(*hash, 0, sizeof(Transcript *) * hash_size); // Set all pointers to NULL.
    return 0;
}

static int add2hash(Transcript **hash, unsigned long hash_size, char *transcript, char *type, char strand, char *ref_name, unsigned long start, unsigned long end, char phase)
{
    unsigned long hash_value = ElfHash(transcript) % hash_size;
    // printf("transcript id: %s, hash value: %lu, type: %s\n", transcript, hash_value, type);
    Transcript *transcript_node = NULL;
    Transcript *transcript_node_ = hash[hash_value];

    while (transcript_node_)
    {
        if (!strcmp(transcript_node_->transcript, transcript))
        {
            transcript_node = transcript_node_; // If transcript_node has been defined, transcript_node = !NULL.
            break;
        }
        else if (transcript_node_->next)
        {
            transcript_node_ = transcript_node_->next;
        }
        else
        {
            break;
        };
    };

    if (!transcript_node) // Add Transcript node.
    {
        transcript_node = malloc(sizeof(Transcript));
        transcript_node->transcript = malloc(strlen(transcript) + 1);
        strcpy(transcript_node->transcript, transcript);
        transcript_node->ref_name = malloc(strlen(ref_name) + 1);
        strcpy(transcript_node->ref_name, ref_name);
        transcript_node->strand = strand;
        transcript_node->exon_number = 0;
        transcript_node->cds_number = 0;
        transcript_node->next = NULL;
        transcript_node->element = NULL;

        if (transcript_node_)
        {
            transcript_node_->next = transcript_node;
        }
        else
        {
            hash[hash_value] = transcript_node;
        };
    };

    // Add exon or cds.
    Element *element_node = malloc(sizeof(Element)); // New element node.
    if (strcmp(type, "CDS"))
    {
        element_node->type = 0; // exon
        transcript_node->exon_number += 1;
    }
    else
    {
        element_node->type = 1; // cds
        transcript_node->cds_number += 1;
    };
    if (phase == '0' || phase == '.')
    {
        element_node->phase = 0;
    }
    else if (phase == '1')
    {
        element_node->phase = 1;
    }
    else
    {
        element_node->phase = 2;
    };
    element_node->start = start;
    element_node->end = end;
    element_node->next = NULL;

    Element *element_node_ = transcript_node->element;
    if (element_node_)
    {
        while (element_node_->next)
        {
            element_node_ = element_node_->next;
        };
        element_node_->next = element_node;
    }
    else
    {
        transcript_node->element = element_node;
    };
    return 0;
}

int read_gff_file(char *file, Transcript ***hash, size_t hash_size)
{
    char buffer[LINE];
    char ref_name[LINE]; // 0
    char type[LINE]; // 2
    unsigned long start; // 3
    unsigned long end; // 4
    char strand; // 6
    char phase; // 7
    char attributes[LINE]; // 8
    char transcript[LINE]; // 8 Parent=
    char *sep = NULL;
    char *attributes_pointer = NULL;

    size_t buffer_size = 0;
    bool new_line = 1;
    const char *format = "%s\t%*[^\t]\t%s\t%lu\t%lu\t%*s\t%c\t%c\t%[^\n]";
    
    create_hash(hash, hash_size);

    FILE *open_file = fopen(file, "r");
    while (fgets(buffer, LINE, open_file))
    {
        if ((buffer[0] != '#') && new_line)
        {
            sscanf(buffer, format, ref_name, type, &start, &end, &strand, &phase, attributes);
            if (!strcmp(type, "CDS") || !strcmp(type, "exon"))
            {
                attributes_pointer = attributes;
                while ((sep = strsep(&attributes_pointer, ";")))
                {
                    if(!strncmp(sep, "Parent=", 7))
                    {
                        strcpy(transcript, sep + 7);
                        break;
                    };
                };
                // printf("type: %s, transcript: %s, strand: %c, ref_name: %s, start: %lu, end: %lu, attributes: %s\n", type, transcript, strand, ref_name, start, end, attributes);
                add2hash(*hash, hash_size, transcript, type, strand, ref_name, start, end, phase);
            };
        };

        buffer_size = strlen(buffer);
        new_line = buffer[buffer_size - 1] == '\n' ? 1 : 0;
    };
    fclose(open_file);
    return 0;
}

int free_gff_hash(Transcript **hash, size_t hash_size)
{
    Transcript *transcript_node = NULL;
    Transcript *transcript_node_ = NULL;
    Element *element_node = NULL;
    Element *element_node_ = NULL;
    for (size_t hash_index = 0; hash_index < hash_size; hash_index++)
    {
        transcript_node = hash[hash_index];
        while (transcript_node)
        {
            transcript_node_ = transcript_node->next;
            free(transcript_node->transcript);
            free(transcript_node->ref_name);
            element_node = transcript_node->element;
            while (element_node)
            {
                element_node_ = element_node->next;
                free(element_node);
                element_node = element_node_;
            };
            free(transcript_node);
            transcript_node = transcript_node_;
        };
    };
    free(hash);
    return 0;
}