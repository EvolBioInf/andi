/**
 * @file
 *
 * This is the main file. It contains functions to parse the commandline
 * arguments, read files etc.
 *
 * @brief The main file
 * @author Fabian Klötzl
 *
 * @section License
 *
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

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
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
long unsigned int BOOTSTRAP = 0;
double RANDOM_ANCHOR_PROP = 0.05;
gsl_rng *RNG = NULL;
int MODEL = M_JC;

void usage(void);
void version(void);

/**
 * @brief The main function.
 *
 * The main function reads and parses the commandline arguments. Depending on
 * the set flags it reads the input files and forwards the contained sequences
 * to processing.
 */
int main(int argc, char *argv[]) {
	int c;
	int version_flag = 0;

	struct option long_options[] = {{"version", no_argument, &version_flag, 1},
									{"help", no_argument, NULL, 'h'},
									{"verbose", no_argument, NULL, 'v'},
									{"join", no_argument, NULL, 'j'},
									{"low-memory", no_argument, NULL, 'l'},
									{"threads", required_argument, NULL, 't'},
									{"bootstrap", required_argument, NULL, 'b'},
									{"model", required_argument, NULL, 'm'},
									{0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	// parse arguments
	while (1) {

		int option_index = 0;

		c = getopt_long(argc, argv, "jvht:p:m:b:l", long_options,
						&option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: break;
			case 'h': usage(); break;
			case 'v':
				FLAGS |= FLAGS & F_VERBOSE ? F_EXTRA_VERBOSE : F_VERBOSE;
				break;
			case 'p': {
				errno = 0;
				char *end;
				double prop = strtod(optarg, &end);

				if (errno || end == optarg || *end != '\0') {
					warnx(
						"Expected a floating point number for -p argument, but "
						"'%s' was given. Skipping argument.",
						optarg);
					break;
				}

				if (prop < 0.0 || prop > 1.0) {
					warnx("A probability should be a value between 0 and 1; "
						  "Ignoring -p %f argument.",
						  prop);
					break;
				}

				RANDOM_ANCHOR_PROP = prop;
				break;
			}
			case 'l': FLAGS |= F_LOW_MEMORY; break;
			case 'j': FLAGS |= F_JOIN; break;
			case 't': {
#ifdef _OPENMP
				errno = 0;
				char *end;
				long unsigned int threads = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0') {
					warnx("Expected a number for -t argument, but '%s' was "
						  "given. Ignoring -t argument.",
						  optarg);
					break;
				}

				if (threads > (long unsigned int)omp_get_num_procs()) {
					warnx(
						"The number of threads to be used, is greater than the "
						"number of available processors; Ignoring -t %lu "
						"argument.",
						threads);
					break;
				}

				THREADS = threads;
#else
				warnx(
					"This version of andi was built without OpenMP and thus "
					"does not support multi threading. Ignoring -t argument.");
#endif
				break;
			}
			case 'b': {
				errno = 0;
				char *end;
				long unsigned int bootstrap = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0' || bootstrap == 0) {
					warnx(
						"Expected a positive number for -b argument, but '%s' "
						"was given. Ignoring -b argument.",
						optarg);
					break;
				}

				BOOTSTRAP = bootstrap - 1;
				break;
			}
			case 'm': {
				// valid options are 'RAW' and 'JC'
				if (strcasecmp(optarg, "RAW") == 0) {
					MODEL = M_RAW;
				} else if (strcasecmp(optarg, "JC") == 0) {
					MODEL = M_JC;
				} else if (strcasecmp(optarg, "KIMURA") == 0) {
					MODEL = M_KIMURA;
				} else {
					warnx("Ignoring argument for --model. Expected Raw, JC or "
						  "Kimura");
				}
				break;
			}
			case '?': /* intentional fall-through */
			default: usage(); break;
		}
	}

	if (version_flag) {
		version();
	}

	argc -= optind;
	argv += optind;

	// at least one file name must be given
	if (FLAGS & F_JOIN && argc == 0) {
		errx(1, "In join mode at least one filename needs to be supplied.");
	}

	dsa_t dsa;
	if (dsa_init(&dsa)) {
		errx(errno, "Out of memory.");
	}

	const char *file_name;

	// parse all files
	int minfiles = FLAGS & F_JOIN ? 2 : 1;
	for (;; minfiles--) {
		if (!*argv) {
			if (minfiles <= 0) break;

			// if no files are supplied, read from stdin
			file_name = "-";
		} else {
			file_name = *argv++;
		}

		if (FLAGS & F_JOIN) {
			read_fasta_join(file_name, &dsa);
		} else {
			read_fasta(file_name, &dsa);
		}
	}

	size_t n = dsa_size(&dsa);

	if (n < 2) {
		errx(1,
			 "I am truly sorry, but with less than two sequences (%zu given) "
			 "there is nothing to compare.",
			 n);
	}

	if (FLAGS & F_VERBOSE) {
		fprintf(stderr, "Comparing %zu sequences\n", n);
		fflush(stderr);
	}

	RNG = gsl_rng_alloc(gsl_rng_default);
	if (!RNG) {
		err(1, "RNG allocation failed.");
	}

	// seed the random number generator with the current time
	gsl_rng_set(RNG, time(NULL));


	// The sequence length is validated in seq_init.
	if (FLAGS & F_SHORT) {
		warnx("One of the given input sequences is shorter than a thousand "
			  "nucleotides. This may result in inaccurate distances. Try an "
			  "alignment instead.");
	}

	// Warn about non ACGT residues.
	if (FLAGS & F_NON_ACGT) {
		warnx("The input sequences contained characters other than acgtACGT. "
			  "These were automatically stripped to ensure correct results.");
	}

	// compute distance matrix
	calculate_distances(dsa_data(&dsa), n);

	dsa_free(&dsa);
	gsl_rng_free(RNG);
	return 0;
}

/**
 * Prints the usage to stdout and then exits successfully.
 */
void usage(void) {
	const char str[] = {
		"Usage: andi [-jlv] [-b INT] [-p FLOAT] [-m MODEL] [-t INT] FILES...\n"
		"\tFILES... can be any sequence of FASTA files. If no files are "
		"supplied, stdin is used instead.\n"
		"Options:\n"
		"  -b, --bootstrap <INT> \n"
		"                    Print additional bootstrap matrices\n"
		"  -j, --join        Treat all sequences from one file as a single "
		"genome\n"
		"  -l, --low-memory  Use less memory at the cost of speed\n"
		"  -m, --model <Raw|JC|Kimura>\n"
		"                    Pick an evolutionary model; default: JC\n"
		"  -p <FLOAT>        Significance of an anchor pair; default: 0.05\n"
		"  -v, --verbose     Prints additional information\n"
#ifdef _OPENMP
		"  -t, --threads <INT> \n"
		"                    The number of threads to be used; by default, all "
		"available processors are used\n"
#endif
		"  -h, --help        Display this help and exit\n"
		"      --version     Output version information and acknowledgments\n"};

	printf("%s", str);
	exit(EXIT_SUCCESS);
}

/**
 * This function just prints the version string and then aborts
 * the program. It conforms to the [GNU Coding
 * Standard](http://www.gnu.org/prep/standards/html_node/_002d_002dversion.html#g_t_002d_002dversion).
 */
void version(void) {
	const char str[] = {
		"andi " VERSION "\n"
		"Copyright (C) 2014 - 2016 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"
		"Acknowledgments:\n"
		"1) Andi: Haubold, B. Klötzl, F. and Pfaffelhuber, P. (2015). "
		"Fast and accurate estimation of evolutionary distances between "
		"closely related genomes\n"
		"2) Algorithms: Ohlebusch, E. (2013). Bioinformatics Algorithms. "
		"Sequence Analysis, Genome Rearrangements, and Phylogenetic "
		"Reconstruction. pp 118f.\n"
		"3) SA construction: Mori, Y. (2005). Short description of improved "
		"two-stage suffix sorting algorithm. "
		"http://homepage3.nifty.com/wpage/software/itssort.txt\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
