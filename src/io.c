/**
 * @file
 * @brief This file contains the definitions for various io methods.
 */
#define _GNU_SOURCE
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <pfasta.h>
#include <compat-string.h>

#include "global.h"
#include "io.h"

/**
 * @brief Joins all sequences from a file into a single long sequence.
 *
 * Apart from reading all sequences from a file, this function also
 * merges them into one long sequence.
 *
 * "I didn't learn joined up handwriting for nothing, you know."
 * ~ Gilderoy Lockhart
 *
 * @param file_name - The name of the file to be used for reading. The name is
 *  also used to infer the sequence name.
 * @param dsa - (output parameter) An array that holds found sequences.
 */
void read_fasta_join(const char *file_name, dsa_t *dsa) {
	if (!file_name || !dsa) return;

	dsa_t single;
	dsa_init(&single);
	read_fasta(file_name, &single);

	if (dsa_size(&single) == 0) {
		return;
	}

	seq_t joined = dsa_join(&single);

	/* In join mode we try to be clever about the sequence name. Given the file
	 * path we extract just the file name. ie. path/file.ext -> file
	 * This obviously fails on Windows.
	 */

	const char *left = strrchr(file_name, '/'); // find the last path separator
	left = (left == NULL) ? file_name : left + 1;
	// left is the position one of to the right of the path separator

	const char *dot = strchrnul(left, '.'); // find the extension

	// copy only the file name, not its path or extension
	joined.name = strndup(left, dot - left);
	CHECK_MALLOC(joined.name);

	dsa_push(dsa, joined);
	dsa_free(&single);
}

/**
 * @brief This function reads sequences from a file.
 * @param file_name - The file to read.
 * @param dsa - (output parameter) An array that holds found sequences.
 */
void read_fasta(const char *file_name, dsa_t *dsa) {
	if (!file_name || !dsa) return;

	int file_descriptor =
		strcmp(file_name, "-") ? open(file_name, O_RDONLY) : STDIN_FILENO;

	if (file_descriptor < 0) {
		warn("%s", file_name);
		return;
	}

	int l;
	int check;

	seq_t top = {};
	pfasta_file pf;

	if ((l = pfasta_parse(&pf, file_descriptor)) != 0) {
		warnx("%s: %s", file_name, pfasta_strerror(&pf));
		goto fail;
	}

	pfasta_seq ps;
	while ((l = pfasta_read(&pf, &ps)) == 0) {
		check = seq_init(&top, ps.seq, ps.name);

		// skip broken sequences
		if (check != 0) continue;

		dsa_push(dsa, top);
		pfasta_seq_free(&ps);
	}

	if (l < 0) {
		warnx("%s: %s", file_name, pfasta_strerror(&pf));
		pfasta_seq_free(&ps);
	}

fail:
	pfasta_free(&pf);
	close(file_descriptor);
}

/**
 * @brief Prints the distance matrix.
 *
 * This function pretty prints the distance matrix. For small distances
 * scientific notation is used.
 *
 * @param D - The distance matrix
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 * @param warnings - Print warnings? Set to 0 for bootstrapped matrices.
 */
void print_distances(const struct model *D, const seq_t *sequences, size_t n,
					 int warnings) {
	size_t i, j;
	int use_scientific = 0;

	double *DD = malloc(n * n * sizeof(*DD));
	CHECK_MALLOC(DD);

#define DD(X, Y) (DD[(X)*n + (Y)])

	typedef double(estimate_fn)(const model *);
	estimate_fn *estimate;

	switch (MODEL) {
		case M_RAW: estimate = &estimate_RAW; break;
		default:
		/* intentional fall-through. This is just here to silence any
		 * compiler warnings. The real default is set in andi.c.*/
		case M_JC: estimate = &estimate_JC; break;
		case M_KIMURA: estimate = &estimate_KIMURA; break;
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			model datum = D(i, j);

			if (!(FLAGS & F_EXTRA_VERBOSE)) {
				datum = model_average(&D(i, j), &D(j, i));
			}

			double dist = DD(i, j) = i == j ? 0.0 : estimate(&datum);
			double coverage = model_coverage(&datum);

			if (isnan(dist) && warnings) {
				const char str[] = {
					"For the two sequences '%s' and '%s' the distance "
					"computation failed and is reported as nan. "
					"Please refer to the documentation for further details."};
				warnx(str, sequences[i].name, sequences[j].name);
			} else if (dist > 0 && dist < 0.001) {
				use_scientific = 1;
			} else if (i < j && coverage < 0.05 && warnings) {
				const char str[] = {
					"For the two sequences '%s' and '%s' less than 5%% "
					"homology were found (%f and %f, respectively)."};
				warnx(str, sequences[i].name, sequences[j].name,
					  model_coverage(&D(i, j)), model_coverage(&D(j, i)));
			}
		}
	}

	printf("%zu\n", n);
	for (i = 0; i < n; i++) {
		// Print exactly nine characters of the name. Pad with spaces if
		// necessary.
		printf("%-9.9s", sequences[i].name);

		for (j = 0; j < n; j++) {
			// use scientific notation for small numbers
			printf(use_scientific ? " %1.4e" : " %1.4f", DD(i, j));
		}
		printf("\n");
	}

	free(DD);
}

/**
 * @brief Prints the coverage matrix.
 * @param D - The distance matrix
 * @param n - The number of sequences.
 */
void print_coverages(const struct model *D, size_t n) {
	size_t i, j;
	printf("\nCoverage:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%1.4e ", model_coverage(&D(i, j)));
		}
		printf("\n");
	}
}
