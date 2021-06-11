#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "fasta_reader.h"
#include "gff_reader.h"
#include "hash.h"

int compare (const void *number1, const void *number2)
{
    return (*(size_t *)number1 < *(size_t *)number2) ? -1 : 1;
}

int revcom(char *string, size_t string_length)
{
    char temp;
    for (size_t index = 0; index < ceil(string_length / 2.0); index++)
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
        };
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
        };
    };
    return 0;
}

int output(Transcript **transcript_hash, Sequence **sequence_hash, unsigned long hash_size, size_t max_sequence_length)
{
    char *buffer = malloc(sizeof(char) * (max_sequence_length + 1));
    char *protein = malloc(sizeof(char) * (max_sequence_length / 3 + 1));
    char cds_prefix[3] = "NN\0";
    char protein_prefix[2] = "x\0";
    char codon_hash2amino_acid[65] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLFx";

    for (size_t hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (Transcript *transcript_node = transcript_hash[hash_index]; transcript_node; transcript_node = transcript_node->next)
        {
            char *sequence = find_sequence(sequence_hash, hash_size, transcript_node->ref_name);
            size_t *exon_positions = malloc(sizeof(size_t) * 2 * transcript_node->exon_number);
            size_t exon_index = 0;
            size_t *cds_positions = malloc(sizeof(size_t) * 2 * transcript_node->cds_number);
            size_t cds_index = 0;
            size_t min_cds_start = max_sequence_length;
            size_t max_cds_start = 0;
            short phase = 0;
            for (Element *element_node = transcript_node->element; element_node; element_node = element_node->next)
            {
                if (element_node->type) // cds
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
                    };
                    cds_index += 2;
                }
                else
                {
                    exon_positions[exon_index] = element_node->start;
                    exon_positions[exon_index + 1] = element_node->end;
                    exon_index += 2;
                };
            };
            qsort(exon_positions, 2 * transcript_node->exon_number, sizeof(size_t), compare);
            qsort(cds_positions, 2 * transcript_node->cds_number, sizeof(size_t), compare);
            // printf("transcript: %s, #exon: %lu, #cds: %lu\n", transcript_node->transcript, transcript_node->exon_number, transcript_node->cds_number);
            
            size_t buffer_offset = 0;
            size_t element_length;

            /* Exon */
            for (exon_index = 0; exon_index < 2 * transcript_node->exon_number; exon_index += 2)
            {
                // printf("exon: %lu %lu\n", exon_positions[exon_index], exon_positions[exon_index + 1]);
                element_length = exon_positions[exon_index + 1] - exon_positions[exon_index] + 1;
                memcpy(buffer + buffer_offset, sequence + exon_positions[exon_index] - 1, element_length);
                buffer_offset += element_length;
            };
            buffer[buffer_offset] = 0;
            if (transcript_node->strand == '-')
            {
                revcom(buffer, buffer_offset);
            };
            printf(">%s_transcript\n", transcript_node->transcript);
            printf("%s\n", buffer);
            free(exon_positions);

            if (transcript_node->cds_number)
            {
                /* CDS */
                buffer_offset = 0;
                for (cds_index = 0; cds_index < 2 * transcript_node->cds_number; cds_index += 2)
                {
                    // printf("cds: %lu %lu\n", cds_positions[cds_index], cds_positions[cds_index + 1]);
                    element_length = cds_positions[cds_index + 1] - cds_positions[cds_index] + 1;
                    memcpy(buffer + buffer_offset, sequence + cds_positions[cds_index] - 1, element_length);
                    buffer_offset += element_length;
                };
                buffer[buffer_offset] = 0;
                if (transcript_node->strand == '-')
                {
                    revcom(buffer, buffer_offset);
                };
                cds_prefix[(3 - phase) % 3] = 0;
                printf(">%s_cds\n", transcript_node->transcript);
                printf("%s%s\n", cds_prefix, buffer);

                protein_prefix[!phase ? 0 : 1] = 0;
                size_t protein_offset = 0;
                while (protein_offset < buffer_offset / 3)
                {
                    unsigned long codon_hash_value = CodonHash(buffer + protein_offset * 3);
                    protein[protein_offset] = codon_hash2amino_acid[codon_hash_value < 64 ? codon_hash_value : 64];
                    protein_offset++;
                };
                protein[protein_offset] = 0;
                printf(">%s_protein\n", transcript_node->transcript);
                printf("%s%s\n", protein_prefix, protein);
                free(cds_positions);
            };
        };
    };
    free(buffer);
    free(protein);
    return 0;
}

int print_help()
{
    printf("Usage:\ngene_sequence -fasta FASTA -gff GFF\n");
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    char fasta[FILE_NAME];
    char gff[FILE_NAME];
    if (argc != 5)
    {
        print_help();
    }
    int ARGS = 0;
    for (int arg_index = 0; arg_index < argc - 1; arg_index++)
    {
        if (!strcmp(argv[arg_index], "-fasta"))
        {
            strncpy(fasta, argv[arg_index + 1], FILE_NAME - 1);
            fasta[FILE_NAME - 1] = 0;
            ARGS++;
        }
        else if (!strcmp(argv[arg_index], "-gff"))
        {
            strncpy(gff, argv[arg_index + 1], FILE_NAME - 1);
            gff[FILE_NAME - 1] = 0;
            ARGS++;
        };
    };
    if (ARGS != 2)
    {
        print_help();
    };
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
    size_t hash_size = 12582917ul;
    size_t max_sequence_length = read_fasta_file(fasta, &sequence_hash, hash_size);
    read_gff_file(gff, &transcript_hash, hash_size);
    output(transcript_hash, sequence_hash, hash_size, max_sequence_length);
    free_fasta_hash(sequence_hash, hash_size);
    free_gff_hash(transcript_hash, hash_size);
    return 0;
}