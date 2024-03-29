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

#include "global.h"
#include "io.h"
#include "process.h"
#include "sequence.h"
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Global variables */
int FLAGS = 0;
int THREADS = 1;
long unsigned int BOOTSTRAP = 0;
double ANCHOR_P_VALUE = 0.025;
gsl_rng *RNG = NULL;
int MODEL = M_JC;

void usage(int);
void version(void);

/**
 * @brief The main function.
 *
 * The main function reads and parses the commandline arguments. Depending on
 * the set flags it reads the input files and forwards the contained sequences
 * to processing. Also it verifies the input for correctness and issues warnings
 * and errors.
 */
int main(int argc, char *argv[]) {
	struct option long_options[] = {
		{"version", no_argument, NULL, 0},
		{"truncate-names", no_argument, NULL, 0},
		{"file-of-filenames", required_argument, NULL, 0},
		{"progress", optional_argument, NULL, 0},
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

	enum { P_AUTO, P_NEVER, P_ALWAYS } progress = P_AUTO;

	struct string_vector file_names;
	string_vector_init(&file_names);

	// parse arguments
	while (1) {

		int option_index = 0;
		int c = getopt_long(argc, argv, "jvht:p:m:b:l", long_options,
							&option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: {
				const char *option_str = long_options[option_index].name;
				if (strcasecmp(option_str, "version") == 0) {
					version();
				}
				if (strcasecmp(option_str, "truncate-names") == 0) {
					FLAGS |= F_TRUNCATE_NAMES;
				}
				if (strcasecmp(option_str, "file-of-filenames") == 0) {
					read_into_string_vector(optarg, &file_names);
				}
				if (strcasecmp(option_str, "progress") == 0) {
					if (!optarg || strcasecmp(optarg, "always") == 0) {
						progress = P_ALWAYS;
					} else if (strcasecmp(optarg, "auto") == 0) {
						progress = P_AUTO;
					} else if (strcasecmp(optarg, "never") == 0) {
						progress = P_NEVER;
					} else {
						warnx("invalid argument to --progress '%s'. Expected "
							  "one of 'auto', 'always', or 'never'.",
							  optarg);
					}
				}
				break;
			}
			case 'h': usage(EXIT_SUCCESS); break;
			case 'v':
				FLAGS |= FLAGS & F_VERBOSE ? F_EXTRA_VERBOSE : F_VERBOSE;
				break;
			case 'p': {
				errno = 0;
				char *end;
				double prop = strtod(optarg, &end);

				if (errno || end == optarg || *end != '\0') {
					soft_errx(
						"Expected a floating point number for -p argument, but "
						"'%s' was given. Skipping argument.",
						optarg);
					break;
				}

				if (prop <= 0.0 || prop >= 1.0) {
					soft_errx("A probability should be a value between 0 and "
							  "1, exclusive; Ignoring -p %f argument.",
							  prop);
					break;
				}

				ANCHOR_P_VALUE = prop;
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
					soft_errx(
						"Expected a positive number for -b argument, but '%s' "
						"was given. Ignoring -b argument.",
						optarg);
					break;
				}

				BOOTSTRAP = bootstrap - 1;
				break;
			}
			case 'm': {
				if (strcasecmp(optarg, "RAW") == 0) {
					MODEL = M_RAW;
				} else if (strcasecmp(optarg, "JC") == 0) {
					MODEL = M_JC;
				} else if (strcasecmp(optarg, "KIMURA") == 0) {
					MODEL = M_KIMURA;
				} else if (strcasecmp(optarg, "LOGDET") == 0) {
					MODEL = M_LOGDET;
				} else {
					soft_errx("Ignoring argument for --model. Expected Raw, "
							  "JC, Kimura or LogDet");
				}
				break;
			}
			case '?': /* intentional fall-through */
			default: usage(EXIT_FAILURE); break;
		}
	}

	argc -= optind;
	argv += optind;

	// copy command line arguments into vector
	// std::copy, anyone?
	for (size_t i = 0; i < (unsigned int)argc; i++) {
		string_vector_push_back(&file_names, argv[i]);
	}

	// at least one file name must be given
	if (FLAGS & F_JOIN && string_vector_size(&file_names) == 0) {
		errx(1, "In join mode at least one filename needs to be supplied.");
	}

	size_t minfiles = FLAGS & F_JOIN ? 2 : 1;
	if (string_vector_size(&file_names) < minfiles) {
		// not enough files passed via arguments
		if (!isatty(STDIN_FILENO)) {
			// read from stdin in pipe
			string_vector_push_back(&file_names, "-");
		} else {
			// print a helpful message on './andi' without args
			usage(EXIT_FAILURE);
		}
	}

	// parse fasta files
	dsa_t dsa;
	dsa_init(&dsa);
	for (size_t i = 0; i < string_vector_size(&file_names); i++) {
		char *file_name = string_vector_at(&file_names, i);
		if (FLAGS & F_JOIN) {
			read_fasta_join(file_name, &dsa);
		} else {
			read_fasta(file_name, &dsa);
		}
	}

	string_vector_free(&file_names);

	size_t n = dsa_size(&dsa);

	if (n < 2) {
		errx(1,
			 "I am truly sorry, but with less than two sequences (%zu given) "
			 "there is nothing to compare.",
			 n);
	}

	RNG = gsl_rng_alloc(gsl_rng_default);
	if (!RNG) {
		err(1, "RNG allocation failed.");
	}

	// seed the random number generator with the current time
	// TODO: enable seeding for reproducibility
	gsl_rng_set(RNG, time(NULL));

	// Warn about non ACGT residues.
	if (FLAGS & F_NON_ACGT) {
		warnx("The input sequences contained characters other than acgtACGT. "
			  "These were automatically stripped to ensure correct results.");
	}

	// validate sequence correctness
	const seq_t *seq = dsa_data(&dsa);
	for (size_t i = 0; i < n; ++i, seq++) {
		if ((FLAGS & F_TRUNCATE_NAMES) && strlen(seq->name) > 10) {
			warnx("The sequence name '%s' is longer than ten characters. It "
				  "will be truncated in the output to '%.10s'.",
				  seq->name, seq->name);
		}

		const size_t LENGTH_LIMIT = (INT_MAX - 1) / 2;
		if (seq->len > LENGTH_LIMIT) {
			errx(1, "The sequence %s is too long. The technical limit is %zu.",
				 seq->name, LENGTH_LIMIT);
		}

		if (seq->len == 0) {
			errx(1, "The sequence %s is empty.", seq->name);
		}

		if (seq->len < 1000) {
			FLAGS |= F_SHORT;
		}
	}

	if (FLAGS & F_SHORT) {
		soft_errx(
			"One of the given input sequences is shorter than a thousand "
			"nucleotides. This may result in inaccurate distances. Try an "
			"alignment instead.");
	}

	// determine whether to print a progress bar
	if (progress == P_AUTO) {
		progress = isatty(STDERR_FILENO) ? P_ALWAYS : P_NEVER;
	}
	if (progress == P_ALWAYS) {
		FLAGS |= F_PRINT_PROGRESS;
	}

	// compute distance matrix
	calculate_distances(dsa_data(&dsa), n);

	dsa_free(&dsa);
	gsl_rng_free(RNG);

	return FLAGS & F_SOFT_ERROR ? EXIT_FAILURE : EXIT_SUCCESS;
}

/**
 * @brief Prints the usage and then exits.
 */
void usage(int status) {
	const char str[] = {
		"Usage: andi [OPTIONS...] FILES...\n"
		"\tFILES... can be any sequence of FASTA files.\n"
		"\tUse '-' as file name to read from stdin.\n"
		"Options:\n"
		"  -b, --bootstrap=INT  Print additional bootstrap matrices\n"
		"      --file-of-filenames=FILE  Read additional filenames from FILE; "
		"one per line\n"
		"  -j, --join           Treat all sequences from one file as a single "
		"genome\n"
		"  -l, --low-memory     Use less memory at the cost of speed\n"
		"  -m, --model=MODEL    Pick an evolutionary model of 'Raw', 'JC', "
		"'Kimura', 'LogDet'; default: JC\n"
		"  -p FLOAT             Significance of an anchor; default: 0.025\n"
		"      --progress=WHEN  Print a progress bar 'always', 'never', or "
		"'auto'; default: auto\n"
#ifdef _OPENMP
		"  -t, --threads=INT    Set the number of threads; by default, all "
		"processors are used\n"
#endif
		"      --truncate-names Truncate names to ten characters\n"
		"  -v, --verbose        Prints additional information\n"
		"  -h, --help           Display this help and exit\n"
		"      --version        Output version information and "
		"acknowledgments\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, "%s", str);
	exit(status);
}

/**
 * @brief This function just prints the version string and then aborts
 * the program. It conforms to the [GNU Coding
 * Standard](http://www.gnu.org/prep/standards/html_node/_002d_002dversion.html#g_t_002d_002dversion).
 */
void version(void) {
	const char str[] = {
		"andi " VERSION "\n"
		"Copyright (C) 2014 - 2020 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"
		"Acknowledgments:\n"
		"1) Andi: Haubold, B. Klötzl, F. and Pfaffelhuber, P. (2015). "
		"Fast and accurate estimation of evolutionary distances between "
		"closely related genomes, Bioinformatics.\n"
		"2) Algorithms: Ohlebusch, E. (2013). Bioinformatics Algorithms. "
		"Sequence Analysis, Genome Rearrangements, and Phylogenetic "
		"Reconstruction. pp 118f.\n"
		"3) SA construction: Mori, Y. (2005). libdivsufsort, unpublished.\n"
		"4) Bootstrapping: Klötzl, F. and Haubold, B. (2016). Support Values "
		"for Genome Phylogenies, Life 6.1.\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
