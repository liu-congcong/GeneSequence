#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "fasta_reader.h"
#include "hash.h"

static int create_hash(Sequence ***hash, unsigned long hash_size)
{
    *hash = malloc(sizeof(Sequence *) * hash_size);
    memset(*hash, 0, sizeof(Sequence *) * hash_size); // Set all pointers to NULL.
    return 0;
}

static int add2hash(Sequence **hash, unsigned long hash_size, char *sequence_id, unsigned long sequence_id_length, char *sequence, unsigned long sequence_length)
{
    unsigned long hash_value = ElfHash(sequence_id) % hash_size;
    Sequence *new_node = malloc(sizeof(Sequence));
    new_node->id = malloc(sequence_id_length + 1);
    strcpy(new_node->id, sequence_id);
    new_node->sequence = malloc(sequence_length + 1);
    strcpy(new_node->sequence, sequence);
    new_node->length = sequence_length;
    new_node->next = hash[hash_value];
    hash[hash_value] = new_node;
    return 0;
}

unsigned long read_fasta_file(char *file, Sequence ***hash, unsigned long hash_size)
{
    char buffer[LINE];
    unsigned long buffer_size = 0;
    unsigned long sequence_length = 0;
    char sequence_id[1024];
    unsigned long sequence_id_length = 0;
    char *sequence_id_blank = NULL;
    char *sequence = NULL;
    unsigned long max_sequence_length = 0;

    FILE *open_file = fopen(file, "r");
    fpos_t file_position;
    fgetpos(open_file, &file_position);
    while (fgets(buffer, LINE, open_file))
    {
        buffer_size = strlen(buffer);
        if (buffer[buffer_size - 1] == '\n')
        {
            buffer[buffer_size - 1] = 0;
            buffer_size--;
        };
        if (buffer[0] == '>')
        {
            max_sequence_length = max_sequence_length > sequence_length ? max_sequence_length : sequence_length;
            sequence_length = 0;
        }
        else{
            sequence_length += buffer_size;
        };
    };
    max_sequence_length = max_sequence_length > sequence_length ? max_sequence_length : sequence_length;

    sequence = malloc(sizeof(char) * (max_sequence_length + 1));
    *sequence = 0;
    create_hash(hash, hash_size);

    fsetpos(open_file, &file_position);
    while (fgets(buffer, LINE, open_file))
    {
        buffer_size = strlen(buffer);
        if (buffer[buffer_size - 1] == '\n')
        {
            buffer[buffer_size - 1] = 0;
            buffer_size--;
        };
        if (buffer[0] == '>')
        {
            if (*sequence != 0)
            {
                add2hash(*hash, hash_size, sequence_id, sequence_id_length, sequence, sequence_length);
                *sequence = 0;
            };

            sequence_id_blank = strchr(buffer, ' ');
            if (sequence_id_blank)
            {
                *sequence_id_blank = 0;
            };
            sequence_id_length = buffer_size - 1;
            strcpy(sequence_id, buffer + 1); // >sequence_id\0
            sequence_id[1023] = 0;
            sequence_length = 0;
        }
        else
        {
            strcat(sequence + sequence_length, buffer);
            sequence_length += buffer_size;
        };
    };
    if (*sequence != 0)
    {
        add2hash(*hash, hash_size, sequence_id, sequence_id_length, sequence, sequence_length);
        *sequence = 0;
    };
    free(sequence);
    fclose(open_file);
    return max_sequence_length;
}

int free_fasta_hash(Sequence **hash, unsigned long hash_size)
{
    Sequence *temp = NULL;
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        while (hash[hash_index])
        {
            temp = hash[hash_index]->next;
            free(hash[hash_index]->id);
            free(hash[hash_index]->sequence);
            free(hash[hash_index]);
            hash[hash_index] = temp;
        };
    };
    free(hash);
    return 0;
}

char *find_sequence(Sequence **hash, unsigned long hash_size, char *sequence_id)
{
    unsigned long hash_value = ElfHash(sequence_id) % hash_size;
    Sequence *sequence_node = hash[hash_value];
    while (sequence_node)
    {
        if (!strcmp(sequence_node->id, sequence_id))
        {
            break;
        };
        sequence_node = sequence_node->next;
    };
    return sequence_node ? sequence_node->sequence : NULL;
}