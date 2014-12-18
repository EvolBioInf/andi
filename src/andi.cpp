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
#include <getopt.h>
#include <errno.h>
#include "global.h"
#include "process.h"
#include "io.h"
#include "sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Global variables */
int FLAGS = 0;
int THREADS = 1;
double RANDOM_ANCHOR_PROP = 0.05;

void usage(void);
void version(void);

/**
 * @brief The main function.
 *
 * The main function reads and parses the commandline arguments. Depending on
 * the set flags it reads the input files and forwards the contained sequences
 * to processing.
 */
int main( int argc, char *argv[]){
	int c;
	int version_flag = 0;
	
	static struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{"raw", no_argument, NULL, 'r'},
		{"verbose", no_argument, NULL, 'v'},
		{"join", no_argument, NULL, 'j'},
		{"low-memory", no_argument, NULL, 'm'},
		{0,0,0,0}
	};
	
	// parse arguments
	while( 1 ){
	
		int option_index = 0;
		
		c = getopt_long( argc, argv, "jvhrt:p:m", long_options, &option_index);
		
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
					errno = 0;
					char *end;
					double prop = strtod( optarg, &end);

					if( errno || end == optarg || *end != '\0'){
						warnx(
							"Expected a floating point number for -p argument, but '%s' was given. "
							"Skipping argument.", optarg
						);
						break;
					}

					if( prop < 0.0 || prop > 1.0 ){
						warnx(
							"A probability should be a value between 0 and 1; "
							"Ignoring -p %f argument.", prop
						);
						break;
					}

					RANDOM_ANCHOR_PROP = prop;
				}
				break;
			case 'm':
				FLAGS |= F_LOW_MEMORY;
				break;
			case 'j':
				FLAGS |= F_JOIN;
				break;
#ifdef _OPENMP
			case 't':
				{
					errno = 0;
					char *end;
					long unsigned int threads = strtoul( optarg, &end, 10);

					if( errno || end == optarg || *end != '\0'){
						warnx(
							"Expected a number for -t argument, but '%s' was given. "
							"Ignoring -t argument.", optarg
						);
						break;
					}

					if( threads > (long unsigned int) omp_get_num_procs() ){
						warnx(
							"The number of threads to be used, is greater then the number of available processors; "
							"Ignoring -t %lu argument.", threads
						);
						break;
					}

					THREADS = threads;
				}
				break;
#endif
			case '?': /* intentional fall-through */
			default:
				usage();
				break;
		}
	}
	
	if( version_flag ){
		version();
	}

	argc -= optind;
	argv += optind;

	// at least one file name must be given
	if( FLAGS & F_JOIN && argc == 0 ){
		errx(1, "In join mode at least one filename needs to be supplied.");
	}
	
	dsa_t *dsa = dsa_new();
	FILE *in = NULL;
	const char *name;

	// parse all files
	int minfiles = FLAGS & F_JOIN ? 2 : 1;
	for( ; ; minfiles-- ){
		if( !*argv){
			if( minfiles <= 0) break;

			// if no files are supplied, read from stdin
			in = stdin;
			name = "stdin";
		} else {
			name = *argv++;
			in = fopen( name, "r");
			if( !in) {
				warn("%s", name);
				continue;
			}
		}

		if( FLAGS & F_JOIN){
			joinedRead( in, dsa, name);
		} else {
			readFile( in, dsa);
		}

		fclose(in);
	}


	size_t n = dsa_size( dsa);
	
	if( FLAGS & F_VERBOSE){
		fprintf( stderr, "Comparing %lu sequences\n", n);
		fflush( stderr);
	}
	
	seq_t *sequences = dsa_data( dsa);
	// compute distance matrix
	if( n >= 2){
		calcDistMatrix(sequences, n);
	} else {
		warnx("I am truly sorry, but with less than two sequences (%lu given) there is nothing to compare.", n);
	}

	dsa_free( dsa);
	return 0;
}

/**
 * Prints the usage to stdout and then exits successfully.
 */
void usage(void){
	const char str[]= {
		"Usage: andi [-jrv] [-p FLOAT] FILES...\n"
		"\tFILES... can be any sequence of FASTA files. If no files are supplied, stdin is used instead.\n"
		"Options:\n"
		"  -j, --join        Treat all sequences from one file as a single genome\n"
		"  -m, --low-memory  Use less memory at the cost of speed\n"
		"  -p <FLOAT>        Significance of an anchor pair; default: 0.05\n"
		"  -r, --raw         Calculates raw distances; default: Jukes-Cantor corrected\n"
		"  -v, --verbose     Prints additional information\n"
#ifdef _OPENMP
		"  -t <INT>          The number of threads to be used; default: 1\n"
#endif
		"  -h, --help        Display this help and exit\n"
		"      --version     Output version information and acknowledgments\n"
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
		"Acknowledgments:\n"
		"1) Andi: Haubold, B. Klötzl, F. and Pfaffelhuber, P. (2014). "
		"Fast and accurate estimation of evolutionary distances between closely related genomes\n"
		"2) Algorithm: Ohlebusch, E. (2013). Bioinformatics Algorithms. Sequence Analysis, "
		"Genome Rearrangements, and Phylogenetic Reconstruction. pp 118f.\n"
		"3) RMQ: Fischer, J. and Heun, V. (2007). A new succinct representation of "
		"RMQ-information and improvements in the enhanced suffix array. "
		"Chen, B. Paterson, M., and Zhang, G. (Eds): ESCAPE 2007, LCNS 4614, pp. 459-470.\n"
		"4) SA construction: Mori, Y. (2005). Short description of improved two-stage suffix "
		"sorting algorithm. http://homepage3.nifty.com/wpage/software/itssort.txt\n"
	};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}


