#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include "global.h"
#include "process.h"

#include "sequence.h"


int FLAGS = 0; /* global */
int CORES = 1; /* global */
int STRATEGY = S_SIMPLE; /* global */


#ifdef __cplusplus
extern "C" {
#endif
#include <kseq.h>

// was: KSEQ_INIT(FILE*, read)

KSEQ_INIT(int, read)

#ifdef __cplusplus
}
#endif


#define MAX_SEQUENCES 1000

void usage(void);
int readFile( FILE *in, seq_t *nextSequence);
void version(void);

int main( int argc, char *argv[]){
	int c, i;
	int version_flag = 0;
	
	static struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{0,0,0,0}
	};
	
	// parse arguments
	while( 1 ){
	
		int option_index = 0;
		
		c = getopt_long( argc, argv, "s:vDhrc:", long_options, &option_index);
		
		if( c == -1){
			break;
		} 
	
		switch (c){
			case 0:
				break;
			case 'h':
				usage();
				break;
			case 'r':
				FLAGS |= F_RAW;
				break;
			case 'D':
				FLAGS |= F_SINGLE;
				break;
			case 'v':
				FLAGS |= FLAGS & F_VERBOSE ? F_EXTRA_VERBOSE : F_VERBOSE;
				break;
			case 's':
				if ( strcmp( optarg, "simple") == 0 ){
					STRATEGY = S_SIMPLE;
				} else if ( strcmp( optarg, "inc") == 0 ){
					STRATEGY = S_INC;	
				} else if ( strcmp( optarg, "window") == 0 ){
					STRATEGY = S_WINDOW;
				} else if ( strcmp( optarg, "check") == 0 ){
					STRATEGY = S_CHECK;
				} else if ( strcmp( optarg, "sophisticated") == 0 ){
					STRATEGY = S_SOPHISTICATED;
				} else{
					fprintf(stderr, "error, unknown strategy: %s\n", optarg);
					usage();
				}
				break;
#ifdef _OPENMP
			case 'c':
				CORES = atoi( optarg);
				break;
#endif
			case '?': /* intentional fallthrough */
			default:
				usage();
				break;
		}
	}
	
	if( version_flag ){
		version();
	}
	
	seq_t *sequences = (seq_t*) malloc( MAX_SEQUENCES * sizeof(seq_t));
	assert( sequences );
	
	FILE *in = NULL;
	int n = 0;
	
	// if no files are supplied, read from stdin
	if( optind == argc){
		in = stdin;
		n = readFile( in, sequences );
	}
	
	// parse all files
	for( i=optind; i< argc; i++){
		in = fopen( argv[i], "r");
		if( !in) continue;
		
		n += readFile( in, sequences + n);
		
		assert( n < MAX_SEQUENCES);
		fclose(in);
	}
	
	fprintf( stderr, "Comparing %d sequences\n", n);
	fflush( stderr);
	
	// compute distance matrix
	if( n >= 2){
		printDistMatrix(sequences, n);
	}

	for( i=0; i<n; i++){
		freeSeq( &sequences[i]);
	}
	free( sequences);
	return 0;
}


/**
 * This function reads sequences from a file.
 * @param in The file pointer to read from.
 * @param nextSequence A pointer to the next free beginning of a sequence.
 * @return The number of found sequences.
 */
int readFile( FILE *in, seq_t *nextSequence){
	int n = 0, l, check;
	
	kseq_t *seq = kseq_init(fileno(in));
	
	while( ( l = kseq_read(seq)) >= 0){

		check = asprintf( &nextSequence->S, "%s", seq->seq.s);
		if( check == -1 ){
			continue;
		}
		
		check = asprintf( &nextSequence->name, "%s", seq->name.s);
		if( check == -1 ){
			free( nextSequence->S);
			continue;
		}

		nextSequence->RS = NULL;
		nextSequence++;
		n++;
	}
	
	kseq_destroy(seq);
	return n;
}

/**
 * Prints the usage to stdout. Does not return.
 */
void usage(void){
	const char str[]= {
		"Usage: np [OPTIONS] FILES...\n"
		"\tFILES... can be any sequence of fasta files. If no files are supplied, stdin is used instead.\n"
		"Options:\n"
		"\t[-r  raw distances; default: corrected]\n"
		"\t[-D  single stranded comparison]\n"
		"\t[-v  verbose]\n"
		"\t[-s <simple|inc> strategy]\n"
#ifdef _OPENMP
		"\t[-c <INT>  The number of cores to be used]\n"
#endif
		"\t[-h  display this help]\n"
	};

	printf("%s", str);
	exit(1);
}

void version(void){
	const char str[]= {
		"np v" VERSION  "\n"
		"Copyright (C) 2014 Fabian Kloetzl\n"
		"GPLv3+ license\n"
	};
	printf("%s", str);
	exit(1);
}


