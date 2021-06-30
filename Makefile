CC = cc -std=gnu99

GeneSequence: GeneSequence.c fasta_reader.o fasta_reader.h gff_reader.o gff_reader.h hash.o hash.h
	${CC} -Wall -g -o GeneSequence GeneSequence.c fasta_reader.o gff_reader.o hash.o -lm

fasta_reader.o: fasta_reader.c fasta_reader.h hash.o hash.h
	${CC} -Wall -c -o fasta_reader.o fasta_reader.c

gff_reader.o: gff_reader.c gff_reader.h hash.o hash.h
	${CC} -Wall -c -o gff_reader.o gff_reader.c

hash.o: hash.c hash.h
	${CC} -Wall -c -o hash.o hash.c

clean:
	rm *.o