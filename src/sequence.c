/**
 * @file
 * @brief Sequence utilities
 *
 * This file contains utility functions for working with DNA sequences.
 */
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "sequence.h"
#include <compat-stdlib.h>

void normalize(seq_t *S);
double shustring_cum_prob(size_t x, double g, size_t l);
size_t min_anchor_length(double p, double g, size_t l);

/** Create a new dynamic array for sequences. */
int dsa_init(dsa_t *A) {
	// allocate at least 4 slots so the growth by 1.5 below doesn't get stuck
	// at 3 slots.
	A->data = malloc(sizeof(*A->data) * 4);
	CHECK_MALLOC(A->data);

	A->capacity = 4;
	A->size = 0;
	return 0;
}

/** Add a sequence to an array. */
void dsa_push(dsa_t *A, seq_t S) {
	if (A->size < A->capacity) {
		A->data[A->size++] = S;
	} else {
		// use the near-optimal growth factor of 1.5
		seq_t *ptr = reallocarray(A->data, A->capacity / 2, sizeof(seq_t) * 3);
		CHECK_MALLOC(ptr);

		A->capacity = (A->capacity / 2) * 3;
		A->data = ptr;
		A->data[A->size++] = S;
	}
}

/** Frees the array and all sequences stored within. */
void dsa_free(dsa_t *A) {
	size_t i;
	for (i = 0; i < A->size; i++) {
		seq_free(&A->data[i]);
	}

	free(A->data);
	*A = (dsa_t){};
}

/** Returns the number of sequences stored within an array. */
size_t dsa_size(const dsa_t *A) {
	return A->size;
}

/** Get the raw C array. */
seq_t *dsa_data(dsa_t *A) {
	return A->data;
}

/**
 * @brief Convert an array of multiple sequences into a single sequence.
 *
 * This function joins all sequences contained in an array into one
 * long sequence. The sequences are separated by a `!` character. The
 * caller has to free the initial array.
 *
 * @returns A new sequence representation the union of the array.
 */
seq_t dsa_join(dsa_t *A) {
	seq_t joined = {};

	if (A->size == 0) {
		return joined;
	}

	if (A->size == 1) {
		/* If we are to join just one sequence, _move_ its contents. */
		joined = A->data[0];
		A->data[0] = (seq_t){};
		return joined;
	}

	seq_t *data = A->data;
	seq_t *it = data;

	// Compute the total length
	size_t total = 0, i;
	for (i = 0; i < A->size; i++, it++) {
		total += it->len + 1;
	}

	// A single malloc for the whole new sequence
	char *ptr = malloc(total);
	CHECK_MALLOC(ptr);
	char *next = ptr;

	// Copy all old sequences and add a `!` in between

	it = data;
	memcpy(next, it->S, it->len);
	next += it->len;

	for (i = 1, it++; i < A->size; i++, it++) {
		*next++ = '!';
		memcpy(next, it->S, it->len);
		next += it->len;
	}

	// Don't forget the null byte.
	*next = '\0';

	joined.S = ptr;
	joined.len = total - 1; // subtract the null byte

	return joined;
}

/**
 * @brief Frees the memory of a given sequence.
 * @param S - The sequence to free.
 */
void seq_free(seq_t *S) {
	free(S->S);
	free(S->name);
	*S = (seq_t){};
}

/**
 * @brief Compute the reverse complement.
 * @param str The master string.
 * @param len The length of the master string
 * @return The reverse complement. The caller has to free it!
 */
char *revcomp(const char *str, size_t len) {
	if (!str) return NULL;
	char *rev = malloc(len + 1);
	CHECK_MALLOC(rev);

	char *r = rev;
	const char *s = &str[len - 1];
	rev[len] = '\0';

	do {
		char d;

		switch (*s--) {
			case 'A': d = 'T'; break;
			case 'T': d = 'A'; break;
			case 'G': d = 'C'; break;
			case 'C': d = 'G'; break;
			case '!': d = ';'; break; // rosebud
			default: continue;
		}

		*r++ = d;
	} while (--len);

	return rev;
}

/**
 * @brief This function concatenates the reverse complement to a given master
 * string. A `#` sign is used as a separator.
 * @param s The master string.
 * @param len Its length.
 * @return The newly concatenated string.
 */
char *catcomp(char *s, size_t len) {
	if (!s) return NULL;

	char *rev = revcomp(s, len);

	char *temp = realloc(rev, 2 * len + 2);
	CHECK_MALLOC(temp);

	rev = temp;
	rev[len] = '#';

	memcpy(rev + len + 1, s, len + 1);

	return rev;
}

/**
 * @brief Calculates the GC content of a sequence.
 *
 * This function computes the relative amount of G and C in the total sequence.
 */
double calc_gc(const seq_t *S) {
	size_t GC = 0;
	const char *p = S->S;

	for (; *p; p++) {
		if (*p == 'G' || *p == 'C') {
			GC++;
		}
	}

	return (double)GC / S->len;
}

/** @brief Prepares a sequences to be used as the subject in a comparison. */
int seq_subject_init(seq_subject *S, const seq_t *base) {
	S->gc = calc_gc(base);
	S->RS = catcomp(base->S, base->len);
	if (!S->RS) return 1;
	S->RSlen = 2 * base->len + 1;

	S->threshold = min_anchor_length(ANCHOR_P_VALUE, S->gc, S->RSlen);

	return 0;
}

/** @brief Frees some memory unused for when a sequence is only used as query.
 */
void seq_subject_free(seq_subject *S) {
	free(S->RS);
	S->RS = NULL;
	S->RSlen = 0;
	S->gc = 0.0;
}

/** @brief Initializes a sequences
 *
 * @returns 0 iff successful.
 */
int seq_init(seq_t *S, const char *seq, const char *name) {
	if (!S || !seq || !name) {
		return 1;
	}

	*S = (seq_t){.S = strdup(seq), .name = strdup(name)};

	CHECK_MALLOC(S->S);
	CHECK_MALLOC(S->name);

	normalize(S);

	// recalculate the length because `normalize` might have stripped some
	// characters.
	S->len = strlen(S->S);

	return 0;
}

/**
 * @brief Restricts a sequence characters set to ACGT.
 *
 * This function strips a sequence of non ACGT characters and converts acgt to
 * the upper case equivalent. A flag is set if a non-canonical character was
 * encountered.
 */
void normalize(seq_t *S) {
	char *p, *q;
	char local_non_acgt = 0;
	for (p = q = S->S; *p; p++) {
		switch (*p) {
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case '!': *q++ = *p; break;
			case 'a':
			case 'c':
			case 'g':
			case 't': *q++ = toupper((unsigned char)*p); break;
			default: local_non_acgt = 1; break;
		}
	}
	*q = '\0';
	if (local_non_acgt) {
#pragma omp atomic
		FLAGS |= F_NON_ACGT;
	}
}

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The probability with which an anchor will be created under a
 * random model.
 * @param g - The the relative amount of GC in the subject.
 * @param l - The length of the subject (includes revcomp).
 * @returns The minimum length of an anchor.
 */
size_t min_anchor_length(double p, double g, size_t l) {
	size_t x = 1;

	while (shustring_cum_prob(x, g / 2, l) < 1 - p) {
		x++;
	}

	return x;
}

/**
 * @brief Calculates the binomial coefficient of n and k.
 *
 * We could (and probably should) use gsl_sf_lnchoose(xx,kk) for this.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analogue.
 * @returns (n choose k)
 */
size_t binomial_coefficient(size_t n, size_t k) {
	if (n <= 0 || k > n) {
		return 0;
	}

	if (k == 0 || k == n) {
		return 1;
	}

	if (k > n - k) {
		k = n - k;
	}

	size_t res = 1;

	for (size_t i = 1; i <= k; i++) {
		res *= n - k + i;
		res /= i;
	}

	return res;
}

/**
 * @brief Given `x` this function calculates the probability of a shustring
 * with a length less or equal to `x` under a random model. This means, it is
 * the cumulative probability.
 *
 * Let X be the longest shortest unique substring (shustring) at any position.
 * Then this function computes P{X <= x} with respect to the given parameter
 * set. See Haubold et al. (2009). Note that `x` includes the final mismatch.
 * Thus, `x` is `match length + 1`.
 *
 * @param x - The maximum length of a shustring.
 * @param p - The half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The probability of a certain shustring length.
 */
double shustring_cum_prob(size_t x, double p, size_t l) {
	double xx = (double)x;
	double ll = (double)l;
	size_t k;

	double s = 0.0;

	for (k = 0; k <= x; k++) {
		double kk = (double)k;
		double t = pow(p, kk) * pow(0.5 - p, xx - kk);

		s += pow(2, xx) * (t * pow(1 - t, ll)) *
			 (double)binomial_coefficient(x, k);
		if (s >= 1.0) {
			s = 1.0;
			break;
		}
	}

	return s;
}
