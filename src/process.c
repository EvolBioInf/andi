

/**
 * @file
 * @brief This file contains various distance methods.
 */
#include "process.h"
#include "esa.h"
#include "global.h"
#include "io.h"
#include "model.h"
#include "sequence.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int calculate_bootstrap(const struct model *M, const seq_t *sequences,
						size_t n);

typedef _Bool bool;
#define false 0
#define true !false

/**
 * @brief This structure captures properties of an anchor.
 */
struct anchor {
	/** The position on the subject. */
	size_t pos_S;
	/** The position on the query. */
	size_t pos_Q;
	/** The length of the exact match. */
	size_t length;
};

/**
 * @brief This is a structure of assorted variables needed for anchor finding.
 */
struct context {
	const esa_s *C;
	const char *query;
	size_t query_length;
	size_t threshold;
};

/**
 * @brief Compute the length of the longest common prefix of two strings.
 *
 * @param S - One string.
 * @param Q - Another string.
 * @param remaining - The length of one of the strings.
 * @returns the length of the lcp.
 */
static inline size_t lcp(const char *S, const char *Q, size_t remaining) {
	size_t length = 0;
	while (length < remaining && S[length] == Q[length]) {
		length++;
	}
	return length;
}

/**
 * @brief Check whether the last anchor can be extended by a lucky anchor.
 *
 * Anchors are defined to be unique and of a minimum length. The uniqueness
 * requires us to search throw the suffix array for a second appearance of the
 * anchor. However, if a left anchor is already unique, we could be sloppy and
 * drop the uniqueness criterion for the second anchor. This way we can skip the
 * lookup and just compare characters directly. However, for a lucky anchor the
 * match still has to be longer than the threshold.
 *
 * @param ctx - Matching context of various variables.
 * @param last_match - The last anchor.
 * @param this_match - Input/Output variable for the current match.
 * @returns true iff the current match is a lucky anchor.
 */
static inline bool lucky_anchor(const struct context *ctx,
								const struct anchor *last_match,
								struct anchor *this_match) {

	size_t advance = this_match->pos_Q - last_match->pos_Q;
	size_t gap = this_match->pos_Q - last_match->pos_Q - last_match->length;

	size_t try_pos_S = last_match->pos_S + advance;
	if (try_pos_S >= (size_t)ctx->C->len || gap > ctx->threshold) {
		return false;
	}

	this_match->pos_S = try_pos_S;
	this_match->length =
		lcp(ctx->query + this_match->pos_Q, ctx->C->S + try_pos_S,
			ctx->query_length - this_match->pos_Q);

	return this_match->length >= ctx->threshold;
}

/**
 * @brief Check for a new anchor.
 *
 * Given the current context and starting position check if the new match is an
 * anchor. The latter requires uniqueness and a certain minimum length.
 *
 * @param ctx - Matching context of various variables.
 * @param last_match - (unused)
 * @param this_match - Input/Output variable for the current match.
 * @returns true iff an anchor was found.
 */
static inline bool anchor(const struct context *ctx,
						  const struct anchor *last_match,
						  struct anchor *this_match) {

	lcp_inter_t inter = get_match_cached(ctx->C, ctx->query + this_match->pos_Q,
										 ctx->query_length - this_match->pos_Q);

	this_match->pos_S = ctx->C->SA[inter.i];
	this_match->length = inter.l <= 0 ? 0 : inter.l;
	return inter.i == inter.j && this_match->length >= ctx->threshold;
}

/**
 * @brief Divergence estimation using the anchor technique.
 *
 * The dist_anchor() function estimates the divergence between two
 * DNA sequences. The subject is given as an ESA, whereas the query
 * is a simple string. This function then looks for *anchors* -- long
 * substrings that exist in both sequences. Then it manually checks for
 * mutations between those anchors.
 *
 * @param C - The enhanced suffix array of the subject.
 * @param query - The actual query string.
 * @param query_length - The length of the query string. Needed for speed
 * reasons.
 * @param gc - The relative GC content of the subject.
 * @returns A matrix with estimates of base substitutions.
 */
model dist_anchor(const esa_s *C, const char *query, size_t query_length,
				  size_t threshold) {
	struct model ret = {.seq_len = query_length, .counts = {0}};

	struct anchor this_match = {0};
	struct anchor last_match = {0};
	bool last_was_right_anchor = false;

	struct context ctx = {C, query, query_length, threshold};

	// Iterate over the complete query.
	while (this_match.pos_Q < query_length) {

		// Check for lucky anchors and fall back to normal strategy.
		if (lucky_anchor(&ctx, &last_match, &this_match) ||
			anchor(&ctx, &last_match, &this_match)) {
			// We have reached a new anchor.

			size_t end_S = last_match.pos_S + last_match.length;
			size_t end_Q = last_match.pos_Q + last_match.length;
			// Check if this can be a right anchor to the last one.
			if (this_match.pos_S > end_S &&
				this_match.pos_Q - end_Q == this_match.pos_S - end_S) {

				// classify nucleotides in the left qanchor
				model_count_equal(&ret, query + last_match.pos_Q,
								  last_match.length);

				// Count the SNPs in between.
				model_count(&ret, C->S + end_S, query + end_Q,
							this_match.pos_Q - end_Q);
				last_was_right_anchor = true;
			} else {
				if (last_was_right_anchor) {
					// If the last was a right anchor, but with the current one,
					// we cannot extend, then add its length.
					model_count_equal(&ret, query + last_match.pos_Q,
									  last_match.length);
				} else if (last_match.length >= threshold * 2) {
					// The last anchor wasn't neither a left or right anchor.
					// But, it was as long as an anchor pair. So still count it.
					model_count_equal(&ret, query + last_match.pos_Q,
									  last_match.length);
				}

				last_was_right_anchor = false;
			}

			// Cache values for later
			last_match = this_match;
		}

		// Advance
		this_match.pos_Q += this_match.length + 1;
	}

	// Very special case: The sequences are identical
	if (last_match.length >= query_length) {
		model_count_equal(&ret, query, query_length);
		return ret;
	}

	// We might miss a few nucleotides if the last anchor was also a right
	// anchor. The logic is the same as a few lines above.
	if (last_was_right_anchor) {
		model_count_equal(&ret, query + last_match.pos_Q, last_match.length);
	} else if (last_match.length >= threshold * 2) {
		model_count_equal(&ret, query + last_match.pos_Q, last_match.length);
	}

	return ret;
}

/*
 * Include distMatrix and distMatrixLM.
 */
#define FAST
#include "dist_hack.h"

#undef FAST
#include "dist_hack.h"

/**
 * @brief Calculates and prints the distance matrix
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 */
void calculate_distances(seq_t *sequences, size_t n) {
	struct model *M = NULL;

	// The maximum number of sequences is near 457'845'052.
	size_t intermediate = SIZE_MAX / sizeof(*M) / n;
	if (intermediate < n) {
		size_t root = (size_t)sqrt(SIZE_MAX / sizeof(*M));
		err(1, "Comparison is limited to %zu sequences (%zu given).", root, n);
	}

	M = malloc(n * n * sizeof(*M));
	if (!M) {
		err(errno, "Could not allocate enough memory for the comparison "
				   "matrix. Try using --join or --low-memory.");
	}

	// compute the distances
	if (FLAGS & F_LOW_MEMORY) {
		distMatrixLM(M, sequences, n);
	} else {
		distMatrix(M, sequences, n);
	}

	// print the results
	print_distances(M, sequences, n, 1);

	// print additional information.
	if (FLAGS & F_VERBOSE) {
		print_coverages(M, n);
	}

	// create new bootstrapped distance matrices
	if (BOOTSTRAP) {
		int res = calculate_bootstrap(M, sequences, n);
		if (res) {
			soft_errx("Bootstrapping failed.");
		}
	}

	free(M);
}

/** Yet another hack. */
#define B(X, Y) (B[(X)*n + (Y)])

/** @brief Computes a bootstrap from _pairwise_ alignments.
 *
 * Doing bootstrapping for alignments with only two sequences is easy. It boils
 * down to a simple multi-nomial process over the substitution matrix.
 *
 * @param M - the initial distance matrix
 * @param sequences - a list of the sequences, containing their lengths
 * @param n - the number of sequences
 *
 * The number of bootstrapped distance matrices to print is implicitly
 * passed via the global `BOOTSTRAP` variable.
 *
 * @returns 0 iff successful.
 */
int calculate_bootstrap(const struct model *M, const seq_t *sequences,
						size_t n) {
	if (!M || !sequences || !n) {
		return 1;
	}

	// B is the new bootstrap matrix
	struct model *B = malloc(n * n * sizeof(*B));
	CHECK_MALLOC(B);

	// Compute a number of new distance matrices
	while (BOOTSTRAP--) {
		for (size_t i = 0; i < n; i++) {
			for (size_t j = i; j < n; j++) {
				if (i == j) {
					B(i, j) = (struct model){.seq_len = 1.0, .counts = {1.0}};
					continue;
				}

				// Bootstrapping should only be used with averaged distances.
				model datum = model_average(&M(i, j), &M(j, i));
				datum = model_bootstrap(datum);

				B(j, i) = B(i, j) = datum;
			}
		}

		print_distances(B, sequences, n, 0);
	}

	free(B);
	return 0;
}
