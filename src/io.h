/**
 * @file
 * @brief This header contains function declarations for io procedures.
 */
#ifndef _IO_H_
#define _IO_H_

#include "model.h"
#include "sequence.h"
#include <err.h>
#include <errno.h>
#include <stdio.h>

/**
 * This is a neat hack for dealing with matrices.
 */
#define D(X, Y) (D[(X)*n + (Y)])
#define M(X, Y) (M[(X)*n + (Y)])

void read_fasta(const char *, dsa_t *dsa);
void read_fasta_join(const char *, dsa_t *dsa);

void print_distances(const struct model *, const seq_t *, size_t, int);
void print_coverages(const struct model *, size_t);

/**
 * @brief A dynamically growing structure for file_names.
 */
struct string_vector {
	char **data;
	size_t capacity, size;
};

char *string_vector_at(struct string_vector *, size_t);
char **string_vector_data(struct string_vector *);
void string_vector_free(struct string_vector *);
void string_vector_init(struct string_vector *);
void string_vector_push_back(struct string_vector *, const char *);
void string_vector_emplace_back(struct string_vector *, char *);
size_t string_vector_size(const struct string_vector *);

void read_into_string_vector(const char *, struct string_vector *);

#endif // _IO_H_
