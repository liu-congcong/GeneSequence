#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "fasta_reader.h"
#include "gff_reader.h"
#include "hash.h"

int compare (const void *x1, const void *x2)
{
    return (*(unsigned long *)x1 < *(unsigned long *)x2) ? -1 : 1;
}

int revcom(char *string, unsigned long string_length)
{
    char temp;
    for (unsigned long index = 0; index < ceil(string_length / 2.0); index++)
    {
        temp = string[index];
        switch (string[string_length - 1 - index])
        {
            case 'A':
                string[index] = 'T';
                break;
            case 'C':
                string[index] = 'G';
                break;
            case 'G':
                string[index] = 'C';
                break;
            case 'T':
                string[index] = 'A';
                break;
            case 'a':
                string[index] = 't';
                break;
            case 'c':
                string[index] = 'g';
                break;
            case 'g':
                string[index] = 'c';
                break;
            case 't':
                string[index] = 'a';
                break;
            default:
                string[index] = string[string_length - 1 - index];
                break;
        }
        switch (temp)
        {
            case 'A':
                string[string_length - 1 - index] = 'T';
                break;
            case 'C':
                string[string_length - 1 - index] = 'G';
                break;
            case 'G':
                string[string_length - 1 - index] = 'C';
                break;
            case 'T':
                string[string_length - 1 - index] = 'A';
                break;
            case 'a':
                string[string_length - 1 - index] = 't';
                break;
            case 'c':
                string[string_length - 1 - index] = 'g';
                break;
            case 'g':
                string[string_length - 1 - index] = 'c';
                break;
            case 't':
                string[string_length - 1 - index] = 'a';
                break;
            default:
                string[string_length - 1 - index] = temp;
                break;
        }
    }
    return 0;
}

int output(Transcript **transcript_hash, Sequence **sequence_hash, unsigned long hash_size, unsigned long max_sequence_length, int type)
{
    char *buffer = malloc(sizeof(char) * (max_sequence_length + 1));
    char *protein = malloc(sizeof(char) * (max_sequence_length / 3 + 1));
    char cds_prefix[3] = "NN\0";
    char protein_prefix[2] = "x\0";
    char codon_hash2amino_acid[65] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLFx";

    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (Transcript *transcript_node = transcript_hash[hash_index]; transcript_node; transcript_node = transcript_node->next)
        {
            char *sequence = find_sequence(sequence_hash, hash_size, transcript_node->ref_name);
            if (sequence)
            {
                unsigned long *exon_positions = malloc(sizeof(unsigned long) * 2 * transcript_node->exon_number);
                unsigned long exon_index = 0;
                unsigned long *cds_positions = malloc(sizeof(unsigned long) * 2 * transcript_node->cds_number);
                unsigned long cds_index = 0;
                unsigned long min_cds_start = max_sequence_length;
                unsigned long max_cds_start = 0;
                short phase = 0;
                for (Element *element_node = transcript_node->element; element_node; element_node = element_node->next)
                {
                    if ((element_node->type == 'c') && type) // cds
                    {
                        cds_positions[cds_index] = element_node->start;
                        cds_positions[cds_index + 1] = element_node->end;
                        if ((transcript_node->strand == '+') && (min_cds_start > cds_positions[cds_index]))
                        {
                            phase = element_node->phase;
                            min_cds_start = cds_positions[cds_index];
                        }
                        else if ((transcript_node->strand == '-') && (max_cds_start < cds_positions[cds_index]))
                        {
                            phase = element_node->phase;
                            max_cds_start = cds_positions[cds_index];
                        }
                        cds_index += 2;
                    }
                    else if ((element_node->type == 'e') && (!type))
                    {
                        exon_positions[exon_index] = element_node->start;
                        exon_positions[exon_index + 1] = element_node->end;
                        exon_index += 2;
                    }
                }
                if (!type)
                {
                    qsort(exon_positions, 2 * transcript_node->exon_number, sizeof(unsigned long), compare);
                }
                else
                {
                    qsort(cds_positions, 2 * transcript_node->cds_number, sizeof(unsigned long), compare);
                }
                // printf("transcript: %s, #exon: %lu, #cds: %lu\n", transcript_node->transcript, transcript_node->exon_number, transcript_node->cds_number);

                unsigned long element_length;
                unsigned long buffer_offset;

                if (!type) // transcript
                {
                    /* Exon */
                    buffer_offset = 0;
                    for (exon_index = 0; exon_index < 2 * transcript_node->exon_number; exon_index += 2)
                    {
                        // printf("exon: %lu %lu\n", exon_positions[exon_index], exon_positions[exon_index + 1]);
                        element_length = exon_positions[exon_index + 1] - exon_positions[exon_index] + 1;
                        memcpy(buffer + buffer_offset, sequence + exon_positions[exon_index] - 1, element_length);
                        buffer_offset += element_length;
                    }
                    buffer[buffer_offset] = 0;
                    if (transcript_node->strand == '-')
                    {
                        revcom(buffer, buffer_offset);
                    }
                    printf(">%s\n", transcript_node->transcript);
                    printf("%s\n", buffer);
                }
                else if (type && transcript_node->cds_number)
                {
                    /* CDS */
                    buffer_offset = 0;
                    for (cds_index = 0; cds_index < 2 * transcript_node->cds_number; cds_index += 2)
                    {
                        // printf("cds: %lu %lu\n", cds_positions[cds_index], cds_positions[cds_index + 1]);
                        element_length = cds_positions[cds_index + 1] - cds_positions[cds_index] + 1;
                        memcpy(buffer + buffer_offset, sequence + cds_positions[cds_index] - 1, element_length);
                        buffer_offset += element_length;
                    }
                    buffer[buffer_offset] = 0;
                    if (transcript_node->strand == '-')
                    {
                        revcom(buffer, buffer_offset);
                    }
                    if (type == 1)
                    {
                        cds_prefix[(3 - phase) % 3] = 0;
                        printf(">%s\n", transcript_node->transcript);
                        printf("%s%s\n", cds_prefix, buffer);
                    }
                    else
                    {
                        protein_prefix[phase ? 1 : 0] = 0;
                        unsigned long protein_offset = 0;
                        while (protein_offset < buffer_offset / 3)
                        {
                            unsigned long codon_hash_value = CodonHash(buffer + protein_offset * 3);
                            protein[protein_offset] = codon_hash2amino_acid[codon_hash_value < 64 ? codon_hash_value : 64];
                            protein_offset++;
                        }
                        protein[protein_offset] = 0;
                        printf(">%s\n", transcript_node->transcript);
                        printf("%s%s\n", protein_prefix, protein);
                    }
                }
                free(exon_positions);
                free(cds_positions);
            }
        }
    }
    free(buffer);
    free(protein);
    return 0;
}

int print_help()
{
    printf("Usage:\nGeneSequence -fasta FASTA -gff GFF -type {transcript | cds | protein}.\n");
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    char fasta[FILE_NAME];
    char gff[FILE_NAME];
    int type; // 0: transcript, 1: cds, 2: protein
    int nesscessary_parameters = 0;
    for (int arg_index = 1; arg_index < argc - 1; arg_index += 2)
    {
        if (!strcmp(argv[arg_index], "-fasta"))
        {
            strncpy(fasta, argv[arg_index + 1], FILE_NAME - 1);
            fasta[FILE_NAME - 1] = 0;
            nesscessary_parameters++;
        }
        else if (!strcmp(argv[arg_index], "-gff"))
        {
            strncpy(gff, argv[arg_index + 1], FILE_NAME - 1);
            gff[FILE_NAME - 1] = 0;
            nesscessary_parameters++;
        }
        else if (!strcmp(argv[arg_index], "-type"))
        {
            if (!strcmp(argv[arg_index + 1], "transcript"))
            {
                type = 0;
                nesscessary_parameters++;
            }
            else if (!strcmp(argv[arg_index + 1], "cds"))
            {
                type = 1;
                nesscessary_parameters++;
            }
            else if (!strcmp(argv[arg_index + 1], "protein"))
            {
                type = 2;
                nesscessary_parameters++;
            }
        }
    }
    if (nesscessary_parameters != 3)
    {
        print_help();
    }
    //char *fasta = "fasta";
    //char *gff = "gff";
    Transcript **transcript_hash = NULL;
    Sequence **sequence_hash = NULL;
    /*
    53ul, 97ul, 193ul, 389ul, 769ul, 1543ul, 3079ul, 6151ul, 12289ul
    24593ul, 49157ul, 98317ul, 196613ul, 393241ul, 786433ul, 1572869ul,
    3145739ul, 6291469ul, 12582917ul, 25165843ul, 50331653ul, 100663319ul,
    201326611ul, 402653189ul, 805306457ul, 1610612741ul, 3221225473ul, 4294967291ul
    */
    unsigned long hash_size = 12582917ul;
    unsigned long max_sequence_length = read_fasta_file(fasta, &sequence_hash, hash_size);
    read_gff_file(gff, &transcript_hash, hash_size);
    output(transcript_hash, sequence_hash, hash_size, max_sequence_length, type);
    free_fasta_hash(sequence_hash, hash_size);
    free_gff_hash(transcript_hash, hash_size);
    return 0;
}