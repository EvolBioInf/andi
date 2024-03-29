/**
 * @file
 * @brief This file contains the definitions for various io methods.
 */
#define _GNU_SOURCE
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <compat-stdlib.h>
#include <compat-string.h>
#include <pfasta.h>

#include "global.h"
#include "io.h"

/**
 * @brief Access an element.
 * @param sv - The base vector.
 * @param index - The element index to access.
 * @returns the string at position `index`.
 */
char *string_vector_at(struct string_vector *sv, size_t index) {
	return sv->data[index];
}

/**
 * @brief Access the underlying buffer.
 * @param sv - The base vector.
 * @returns the underlying buffer.
 */
char **string_vector_data(struct string_vector *sv) {
	return sv->data;
}

/**
 * @brief Free all data.
 * @param sv - The base vector.
 */
void string_vector_free(struct string_vector *sv) {
	size_t i = 0;
	for (; i < sv->size; i++) {
		free(sv->data[i]);
	}
	free(sv->data);
}

/**
 * @brief Initialise the vector.
 * @param sv - The base vector.
 */
void string_vector_init(struct string_vector *sv) {
	sv->data = malloc(sizeof(*sv->data) * 4);
	CHECK_MALLOC(sv->data);

	sv->capacity = 4;
	sv->size = 0;
}

/**
 * @brief Adds a copy to the end of the vector.
 * @param sv - The base vector.
 * @param str - The new string to add.
 */
void string_vector_push_back(struct string_vector *sv, const char *str) {
	string_vector_emplace_back(sv, strdup(str));
}

/**
 * @brief Add a file name to the end of the vector, directly.
 * @param sv - The base vector.
 * @param str - The string to emplace.
 */
void string_vector_emplace_back(struct string_vector *sv, char *str) {
	if (sv->size < sv->capacity) {
		sv->data[sv->size++] = str;
	} else {
		char **ptr = reallocarray(sv->data, sv->capacity / 2, 3 * sizeof(*ptr));
		CHECK_MALLOC(ptr);
		sv->data = ptr;
		sv->capacity = (sv->capacity / 2) * 3;
		sv->data[sv->size++] = str;
	}
}

/**
 * @brief Return the number of elements.
 * @param sv - The base vector.
 * @returns the size of the vector.
 */
size_t string_vector_size(const struct string_vector *sv) {
	return sv->size;
}

/**
 * @brief Read a *fof* and add its contents into a vector.
 * @param file_name - The file of file names aka. fof.
 * @param sv - The vector to add file names to.
 */
void read_into_string_vector(const char *file_name, struct string_vector *sv) {
	FILE *file = strcmp(file_name, "-") ? fopen(file_name, "r") : stdin;
	if (!file) {
		soft_err("%s", file_name);
		return;
	}

	while (1) {
		char *str = NULL;
		size_t buffer_size = 0;
		ssize_t check = getline(&str, &buffer_size, file);

		// EOF is set only *after* getline tried to read past it.
		if (check == -1 && feof(file) != 0) {
			free(str);
			break; // EOF
		}

		if (check == -1) {
			soft_err("%s", file_name);
			break;
		}

		char *nl = strchr(str, '\n');
		if (nl) {
			*nl = '\0'; // remove newline character
		}

		// ignore empty lines
		if (strlen(str) == 0) {
			free(str);
			continue;
		}

		string_vector_emplace_back(sv, str);
	}

	int check = fclose(file);
	if (check != 0) {
		soft_err("%s", file_name);
	}
}

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
		soft_err("%s", file_name);
		return;
	}

	struct pfasta_parser pp = pfasta_init(file_descriptor);
	if (pp.errstr) {
		soft_errx("%s: %s", file_name, pp.errstr);
		goto fail;
	}

	seq_t top = {};
	while (!pp.done) {
		struct pfasta_record pr = pfasta_read(&pp);
		if (pp.errstr) {
			soft_errx("%s: %s", file_name, pp.errstr);
			goto fail;
		}

		int check = seq_init(&top, pr.sequence, pr.name);

		// skip broken sequences
		if (check != 0) continue;

		dsa_push(dsa, top);
		pfasta_record_free(&pr);
	}

fail:
	pfasta_free(&pp);
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
		case M_LOGDET: estimate = &estimate_LOGDET; break;
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			model datum = D(i, j);

			if (!(FLAGS & F_EXTRA_VERBOSE)) {
				datum = model_average(&D(i, j), &D(j, i));
			}

			double dist = DD(i, j) = i == j ? 0.0 : estimate(&datum);

			if (dist > 0 && dist < 0.001) {
				use_scientific = 1;
			}

			if (isnan(dist) && warnings) {
				const char str[] = {
					"For the two sequences '%s' and '%s' the distance "
					"computation failed and is reported as nan. "
					"Please refer to the documentation for further details."};
				soft_errx(str, sequences[i].name, sequences[j].name);
			}

			if (!isnan(dist) && i < j && warnings) {
				double coverage1 = model_coverage(&D(i, j));
				double coverage2 = model_coverage(&D(j, i));

				if (coverage1 < 0.2 || coverage2 < 0.2) {
					const char str[] = {
						"For the two sequences '%s' and '%s' very little "
						"homology was found (%f and %f, respectively)."};
					soft_errx(str, sequences[i].name, sequences[j].name,
							  coverage1, coverage2);
				}
			}
		}
	}

	printf("%zu\n", n);
	for (i = 0; i < n; i++) {
		// Print ten characters of the name. Pad with spaces, if
		// necessary. Truncate to exactly ten characters if requested by user.
		printf(FLAGS & F_TRUNCATE_NAMES ? "%-10.10s" : "%-10s",
			   sequences[i].name);

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
