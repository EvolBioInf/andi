/**
 * @file
 *
 * This is the main file. It contains functions to parse the commandline arguments,
 * read files etc.
 * 
 * @brief The main file
 * @author Fabian Klötzl
 
 * @section License
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include "global.h"
#include "process.h"

#include "sequence.h"

/* Global variables */
int FLAGS = 0;
int THREADS = 1;
double RANDOM_ANCHOR_PROP = 0.95;


#ifdef __cplusplus
extern "C" {
#endif
#include <kseq.h>

// was: KSEQ_INIT(FILE*, read)

KSEQ_INIT(int, read)

#ifdef __cplusplus
}
#endif

void usage(void);
void readFile( FILE *in, dsa_t *dsa);
void joinedRead( FILE *in, dsa_t *dsa, char *name);
void version(void);

/**
 * @brief The main function.
 *
 * The main function reads and parses the commandline arguments. Depending on
 * the set flags it reads the input files and forwards the contained sequences
 * to processing.
 */
int main( int argc, char *argv[]){
	int c, i;
	int version_flag = 0;
	
	static struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{"raw", no_argument, NULL, 'r'},
		{"verbose", no_argument, NULL, 'v'},
		{"join", no_argument, NULL, 'j'},
		{0,0,0,0}
	};
	
	// parse arguments
	while( 1 ){
	
		int option_index = 0;
		
		c = getopt_long( argc, argv, "jvhrt:p:s", long_options, &option_index);
		
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
			case 'v':
				FLAGS |= FLAGS & F_VERBOSE ? F_EXTRA_VERBOSE : F_VERBOSE;
				break;
			case 'p':
				{
					double prop = atof( optarg);
					if( prop < 0.0 || prop > 1.0 ){
						const char str[] = {
							"Warning: A propability should be a value between 0 and 1; "
							"Ignoring -p argument.\n"
						};
						fprintf( stderr, "%s", str);
						break;
					}
					RANDOM_ANCHOR_PROP = prop;
				}
				break;
			case 'j':
				FLAGS |= F_JOIN;
				break;
#ifdef _OPENMP
			case 't':
				{
					int threads = atoi( optarg);
					if( threads < 1 ){
						fprintf( stderr, 
							"Warning: The number of threads should be positive; Ignoring -t argument.\n"
						);
						break;
					}
					THREADS = threads;
				}
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
	
	dsa_t *dsa = init_dsa();
	FILE *in = NULL;
	
	if( FLAGS & F_JOIN ){
		// atleast one file name must be given
		if( optind == argc ){
			fprintf( stderr, "Error: in join mode at least one filename needs to be supplied\n");
			exit(EXIT_FAILURE);
		}
		
		// only one file
		if( optind + 1 == argc ){
			in = stdin;
			joinedRead( in, dsa, strdup("stdin"));
		}
		
		// parse all files
		for( i=optind; i< argc; i++){
			in = fopen( argv[i], "r");
			if( !in) continue;
			
			/* In join mode we try to be clever about the sequence name. Given the file
			 * path we extract just the file name. ie. path/file.ext -> file
			 */
			char *filename = argv[i];
			char *left = strrchr( filename, '/') + 1;
			if( left == NULL ){
				left = filename;
			}
			char *dot = strchrnul( left, '.');
			
			filename = strndup( left, dot-left );

			joinedRead( in, dsa, filename);

			fclose(in);
		}
		
	} else {
		// if no files are supplied, read from stdin
		if( optind == argc){
			in = stdin;
			readFile( in, dsa );
		}
	
		// parse all files
		for( i=optind; i< argc; i++){
			in = fopen( argv[i], "r");
			if( !in) continue;

			readFile( in, dsa);

			fclose(in);
		}
	}
	
	size_t n = size_dsa( dsa);
	
	if( FLAGS & F_VERBOSE){
		fprintf( stderr, "Comparing %lu sequences\n", n);
		fflush( stderr);
	}
	
	seq_t *sequences = data_dsa( dsa);
	// compute distance matrix
	if( n >= 2){
		calcDistMatrix(sequences, n);
	} else {
		fprintf( stderr, "I am truly sorry, but with less than two sequences there is nothing to compare.\n");
	}

	free_dsa( dsa);
	return 0;
}

/**
 * @brief Joins all sequnces from a file into a single long sequence.
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
void joinedRead( FILE *in, dsa_t *dsa, char *name){
	if( !in || !dsa) return;

	dsa_t *single = init_dsa();
	readFile( in, single);
	
	if( size_dsa( single) == 0 ){
		return;
	}
	
	seq_t joined = join_dsa( single);
	joined.name = name;
	push_dsa( dsa, joined);
	free_dsa( single);
}


/**
 * This function reads sequences from a file.
 * @param in - The file pointer to read from.
 * @param dsa - An array that holds found sequences.
 */
void readFile( FILE *in, dsa_t *dsa){
	if( !in || !dsa) return;
	int l;
	seq_t top = { NULL, NULL, 0, 0, NULL, 0.0};
	
	kseq_t *seq = kseq_init(fileno(in));
	
	while( ( l = kseq_read(seq)) >= 0){
		top.S = strdup( seq->seq.s);
		top.name = strdup( seq->name.s);
		top.len = strlen( top.S);
		
		if( !top.S || !top.name ){
			free_seq( &top);
			continue;
		}
		
		push_dsa( dsa, top);
	}
	
	kseq_destroy(seq);
}

/**
 * Prints the usage to stdout and then exits successfully.
 */
void usage(void){
	const char str[]= {
		"Usage: andi [-jrv] [-p FLOAT] FILES...\n"
		"\tFILES... can be any sequence of fasta files. If no files are supplied, stdin is used instead.\n"
		"Options:\n"
		"  -j, --join      Treat all sequences from one file as a single genome\n"
		"  -p <FLOAT>      Certainty that a pair of anchors was not found by chance; default: 0.95\n"
		"  -r, --raw       Calculates raw distances; default: Jukes-Cantor corrected\n"
		"  -v, --verbose   Prints additional information\n"
#ifdef _OPENMP
		"  -t <INT>        The number of threads to be used; default: 1\n"
#endif
		"  -h, --help      display this help and exit\n"
		"      --version   output version information and acknowledgements\n"
	};

	printf("%s", str);
	exit(EXIT_SUCCESS);
}

/**
 * This function just prints the version string and then aborts
 * the program. It conforms to the [GNU Coding Standard](http://www.gnu.org/prep/standards/html_node/_002d_002dversion.html#g_t_002d_002dversion).
 */
void version(void){
	const char str[]= {
		"andi " VERSION  "\n"
		"Copyright (C) 2014 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"
		"Acknowledgements:\n"
		"1) Andi: Haubold, B. Klötzl, F. and Pfaffelhuber, P. (2014)."
		" Fast and accurate estimation of evolutionary distances between genomes. In preparation.\n"
		"2) Algorithm: Ohlebusch, E. (2013). Bioinformatics Algoritms. Sequence Analysis, "
		"Genome Rearrangements, and Phylogenetic Reconstruction. pp 118f.\n"
		"3) RMQ: Fischer, J. and Heun, V. (2007). A new succinct representation of "
		"RMQ-information and improvements in the enhanced suffix array. "
		"Chen, B. Paterson, M., and Zhang, G. (Eds): ESCAPE 2007, LCNS 4614, pp. 459-470.\n"
		"4) SA construction: Mori, Y. (2005). Short description of improved two-stage suffix "
		"sorting alorithm. http://homepage3.nifty.com/wpage/software/itssort.txt\n"
	};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}


