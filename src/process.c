

/**
 * @file
 * @brief This file contains various distance methods.
 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "esa.h"
#include "global.h"
#include "io.h"
#include "model.h"
#include "process.h"
#include "sequence.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double shuprop(size_t x, double g, size_t l);
int calculate_bootstrap(const struct model *M, const seq_t *sequences,
						size_t n);

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The probability with which an anchor is allowed to be random.
 * @param g - The the relative amount of GC in the subject.
 * @param l - The length of the subject.
 * @returns The minimum length of an anchor.
 */
size_t minAnchorLength(double p, double g, size_t l) {
	size_t x = 1;

	double prop = 0.0;
	while (prop < 1 - p) {
		prop = shuprop(x, g / 2, l);
		x++;
	}

	return x;
}

/**
 * @brief Calculates the binomial coefficient of n and k.
 *
 * We used to use gsl_sf_lnchoose(xx,kk) for this functionality.
 * After all, why implement something that has already been done?
 * Well, the reason is simplicity: GSL is used for only this one
 * function and the input (n<=20) is not even considered big.
 * Hence its much easier to have our own implementation and ditch
 * the GSL dependency even if that means our code is a tiny bit
 * less optimized and slower.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analog.
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
 * with a length less than `x`.
 *
 * Let X be the longest shortest unique substring (shustring) at any position.
 * Then this function computes P{X <= x} with respect to the given parameter
 * set. See Haubold et al. (2009).
 *
 * @param x - The maximum length of a shustring.
 * @param g - The the half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The probability of a certain shustring length.
 */
double shuprop(size_t x, double p, size_t l) {
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

/**
 * @brief Divergence estimation using the anchor technique.
 *
 * The dist_anchor() function estimates the divergence between two
 * DNA sequences. The subject is given as an ESA, whereas the query
 * is a simple string. This function then looks for *anchors* -- long
 * substrings that exist in both sequences. Then it manually checks for
 * mutations between those anchors.
 *
 * @return An estimate for the number of mutations within homologous regions.
 * @param C - The enhanced suffix array of the subject.
 * @param query - The actual query string.
 * @param query_length - The length of the query string. Needed for speed
 * reasons.
 */
model dist_anchor(const esa_s *C, const char *query, size_t query_length,
				  double gc) {
	struct model ret = {.seq_len = query_length, .counts = {0}};

	lcp_inter_t inter;

	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	// This variable indicates that the last anchor was the right anchor of a
	// pair.
	size_t last_was_right_anchor = 0;

	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;

	size_t num_right_anchors = 0;

#ifdef DEBUG
	size_t num_matches = 0;
	size_t num_anchors = 0;
	size_t num_anchors_in_rc = 0;
	size_t num_right_anchors_in_rc = 0;
	size_t length_anchors = 0;
	double off_num = 0.0;
	double off_dem = 0.0;
#endif

	size_t threshold =
		minAnchorLength(1 - sqrt(1 - RANDOM_ANCHOR_PROP), gc, C->len);

	// Iterate over the complete query.
	while (this_pos_Q < query_length) {
		inter =
			get_match_cached(C, query + this_pos_Q, query_length - this_pos_Q);

#ifdef DEBUG
		num_matches++;
#endif

		this_length = inter.l <= 0 ? 0 : inter.l;

		if (inter.i == inter.j && this_length >= threshold) {
			// We have reached a new anchor.
			this_pos_S = C->SA[inter.i];

#ifdef DEBUG
			num_anchors++;
			length_anchors += this_length;
			if (this_pos_S < (size_t)(C->len / 2)) {
				num_anchors_in_rc++;
			}
#endif

			// Check if this can be a right anchor to the last one.
			if (this_pos_S > last_pos_S &&
				this_pos_Q - last_pos_Q == this_pos_S - last_pos_S) {
				num_right_anchors++;
#ifdef DEBUG
				if (this_pos_S < (size_t)(C->len / 2)) {
					num_right_anchors_in_rc++;
				}
#endif
				// classify nucleotides in the qanchor
				model_count_equal(&ret, query + last_pos_Q, last_length);

				// Count the SNPs in between.
				model_count(&ret, C->S + last_pos_S + last_length,
							query + last_pos_Q + last_length,
							this_pos_Q - last_pos_Q - last_length);
				last_was_right_anchor = 1;
			} else {
#ifdef DEBUG
				double off = fabs((double)(this_pos_Q - last_pos_Q) -
								  (double)(this_pos_S - last_pos_S));
				if (off < 100) {
					off_num += off;
					off_dem++;
				}
#endif
				if (last_was_right_anchor) {
					// If the last was a right anchor, but with the current one,
					// we cannot extend, then add its length.
					model_count_equal(&ret, C->S + last_pos_S, last_length);
				} else if ((last_length / 2) >= threshold) {
					// The last anchor wasn't neither a left or right anchor.
					// But, it was as long as an anchor pair. So still count it.
					model_count_equal(&ret, C->S + last_pos_S, last_length);
				}

				last_was_right_anchor = 0;
			}

			// Cache values for later
			last_pos_Q = this_pos_Q;
			last_pos_S = this_pos_S;
			last_length = this_length;
		}

		// Advance
		this_pos_Q += this_length + 1;
	}

#ifdef DEBUG
	if (FLAGS & F_EXTRA_VERBOSE) {
		const char str[] = {"- threshold: %ld\n"
							"- matches: %lu\n"
							"- anchors: %lu (rc: %lu)\n"
							"- right anchors: %lu (rc: %lu)\n"
							"- avg length: %lf\n"
							"- off: %f (skipped: %.0f)\n"
							"\n"};

#pragma omp critical
		{
			fprintf(stderr, str, threshold, num_matches, num_anchors,
					num_anchors_in_rc, num_right_anchors,
					num_right_anchors_in_rc,
					(double)length_anchors / num_anchors, off_num / off_dem,
					off_dem);
		}
	}
#endif

	// Very special case: The sequences are identical
	if (last_length >= query_length) {
		model_count(&ret, C->S + last_pos_S, query, query_length);
		return ret;
	}

	// We might miss a few nucleotides if the last anchor was also a right
	// anchor.
	if (last_was_right_anchor) {
		model_count(&ret, C->S + last_pos_S, query + last_pos_Q, last_length);
	}

	return ret;
}

/**
 * @brief Computes the distance matrix.
 *
 * The distMatrix() populates the D matrix with computed distances.
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
#define FAST
#include "dist_hack.h"

/**
 * @brief Computes the distance matrix.
 *
 * The distMatrixLM() populates the D matrix with computed distances.
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
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
			warnx("Bootstrapping failed.");
		}
	}

	free(M);
}

/** Yet another hack. */
#define B(X, Y) (B[(X)*n + (Y)])

/** @brief Computes a bootstrap from _pairwise_ aligments.
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
