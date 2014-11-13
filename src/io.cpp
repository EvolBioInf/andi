/**
 * @file
 * @brief This file contains the definitions for various io methods.
 */
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "io.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <kseq.h>

// was: KSEQ_INIT(FILE*, read)

KSEQ_INIT(int, read)

#ifdef __cplusplus
}
#endif

/**
 * @brief Joins all sequences from a file into a single long sequence.
 *
 * Apart from reading all sequences from a file, this function also
 * merges them into one long sequence.
 *
 * "I didn't learn joined up handwriting for nothing, you know."
 * ~ Gilderoy Lockhart
 *
 * @param in - The file pointer to read from.
 * @param dsa - An array that holds found sequences.
 * @param name - The name of the file to be used for the name of the sequence.
 */
void joinedRead( FILE *in, dsa_t *dsa, const char *name){
	if( !in || !dsa || !name) return;

	dsa_t *single = dsa_new();
	readFile( in, single);
	
	if( dsa_size( single) == 0 ){
		return;
	}
	
	seq_t joined = dsa_join( single);

	/* In join mode we try to be clever about the sequence name. Given the file
	 * path we extract just the file name. ie. path/file.ext -> file
	 * This obviously fails on Windows.
	 */

	const char *left = strrchr( name, '/'); // find the last path separator
	left = (left == NULL) ? name : left + 1; // left is the position one of to the right of the path separator
	
	const char *dot = strchrnul( left, '.'); // find the extension
	
	joined.name = strndup( left, dot-left ); // copy only the file name, not its path or extension
	dsa_push( dsa, joined);
	dsa_free( single);
}


/**
 * This function reads sequences from a file.
 * @param in - The file pointer to read from.
 * @param dsa - An array that holds found sequences.
 */
void readFile( FILE *in, dsa_t *dsa){
	if( !in || !dsa) return;
	int l;
	int check;
	seq_t top = { NULL, NULL, 0, 0, NULL, 0.0};
	
	kseq_t *seq = kseq_init(fileno(in));
	
	while( ( l = kseq_read(seq)) >= 0){
		check = seq_init( &top, seq->seq.s, seq->name.s);

		// skip broken sequences
		if( check != 0) continue;
		
		dsa_push( dsa, top);
	}
	
	kseq_destroy(seq);
}


/**
 * @brief Prints the distance matrix.
 *
 * This function pretty prints the distance matrix. For small distances
 * scientific notation is used.
 * @param D - The distance matrix
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 */
void printDistMatrix( double *D, seq_t *sequences, size_t n){

	int use_scientific = 0;
	int failed = 0;
	size_t i,j;

	for( i=0; i<n && (!use_scientific || !failed); i++){
		for( j=0; j<n; j++){
			if( D(i,j) > 0 && D(i,j) < 0.001 ){
				use_scientific = 1;
			}
			if( isnan(D(i,j))){
				failed = 1;
			}
		}
	}

	if( failed){
		warnx("Some distance computations failed and are reported as nan. "
			"Please refer to the documentation for further details.");
	}

	printf("%lu\n", n);
	for( i=0;i<n;i++){
		// Print exactly nine characters of the name. Pad with spaces if necessary.
		printf("%-9.9s", sequences[i].name);
		
		for( j=0;j<n;j++){
			double avg = (D(i,j) + D(j,i))/2;
			if( use_scientific){
				printf(" %1.4e", avg);
			} else {
				printf(" %1.4f", avg);
			}
		}
		printf("\n");
	}
}

void printCovMatrix( data_t *D, size_t n){
	size_t i,j;
	printf("\nCoverage:\n");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			printf("%1.4e ", D(i,j).coverage);
		}
		printf("\n");
	}
}

