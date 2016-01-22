/** @file
 * @brief This file is a preprocessor hack for the two functions `distMatrix`
 * and `distMatrixLM`.
 */
#ifdef FAST
#define NAME distMatrix
#define P_OUTER _Pragma("omp parallel for num_threads( THREADS)")
#define P_INNER
#else
#undef NAME
#undef P_OUTER
#undef P_INNER
#define NAME distMatrixLM
#define P_OUTER
#define P_INNER _Pragma("omp parallel for num_threads( THREADS)")
#endif

/** @brief This function calls dist_andi for pairs of subjects and queries, and
 *thereby fills the distance matrix.
 *
 * This function is actually two functions. It is one template that gets
 *compiled into two functions via
 * preprocessor hacks. The reason is DRY (Do not Repeat Yourselves).
 *   The two functions only differ by their name and pragmas; i.e. They run in
 *different parallel modes.
 * `distMatrix` is faster than `distMatrixLM` but needs more memory.
 *
 * @param sequences - The sequences to compare
 * @param n - The number of sequences
 * @param M - A matrix for additional output data
 */
void NAME(struct model *M, seq_t *sequences, size_t n) {
	size_t i;

	//#pragma
	P_OUTER
	for (i = 0; i < n; i++) {
		seq_t *subject = &sequences[i];
		esa_s E;

		if (seq_subject_init(subject) || esa_init(&E, subject)) {
			errx(1, "Failed to create index for %s.", subject->name);
		}

		// now compare every other sequence to i
		size_t j;

		P_INNER
		for (j = 0; j < n; j++) {
			if (j == i) {
				M(i, j) = (struct model){.seq_len = 9, .counts = {9}};
				continue;
			}

			// TODO: Provide a nicer progress indicator.
			if (FLAGS & F_EXTRA_VERBOSE) {
#pragma omp critical
				{ fprintf(stderr, "comparing %zu and %zu\n", i, j); }
			}

			size_t ql = sequences[j].len;

			M(i, j) = dist_anchor(&E, sequences[j].S, ql, subject->gc);
		}

		esa_free(&E);
		seq_subject_free(subject);
	}
}
