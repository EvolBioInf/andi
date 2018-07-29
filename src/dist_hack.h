/** @file
 * @brief This file is a preprocessor hack for the two functions `distMatrix`
 * and `distMatrixLM`.
 */
// clang-format off
#ifdef FAST
#define NAME distMatrix
#define P_OUTER _Pragma("omp parallel for num_threads( THREADS) default(none) shared(progress_counter) firstprivate( ANCHOR_P_VALUE, stderr, M, sequences, n, print_progress)")
#define P_INNER
#else
#undef NAME
#undef P_OUTER
#undef P_INNER
#define NAME distMatrixLM
#define P_OUTER
#define P_INNER _Pragma("omp parallel for num_threads( THREADS) default(none) shared(progress_counter) firstprivate( threshold, stderr, M, sequences, n, print_progress, i, E, subject)")
#endif
// clang-format on

/** @brief This function calls dist_andi for pairs of subjects and queries, and
 * thereby fills the distance matrix.
 *
 * This function is actually two functions. It is one template that gets
 * compiled into two functions via preprocessor hacks. The reason is DRY (Do not
 * Repeat Yourselves).
 * The two functions only differ by their name and pragmas; i.e. They run in
 * different parallel modes.
 * `distMatrix` is faster than `distMatrixLM` but needs more memory.
 *
 * @param sequences - The sequences to compare
 * @param n - The number of sequences
 * @param M - A matrix for additional output data
 */
void NAME(struct model *M, const seq_t *sequences, size_t n) {
	size_t i;

	size_t progress_counter = 0;
	int print_progress = FLAGS & F_PRINT_PROGRESS;

	if (print_progress) {
		fprintf(stderr, "Comparing %zu sequences: %5.1f%% (%zu/%zu)", n, 0.0,
				(size_t)0, n * n - n);
	}

	//#pragma
	P_OUTER
	for (i = 0; i < n; i++) {
		seq_subject subject;
		esa_s E;

		if (seq_subject_init(&subject, &sequences[i]) ||
			esa_init(&E, &subject)) {
			errx(1, "Failed to create index for %s.", sequences[i].name);
		}

		size_t threshold = min_anchor_length(ANCHOR_P_VALUE, subject.gc, subject.RSlen);

		// now compare every other sequence to i
		size_t j;

		P_INNER
		for (j = 0; j < n; j++) {
			if (j == i) {
				M(i, j) = (struct model){.seq_len = 9, .counts = {9}};
				continue;
			}

			size_t ql = sequences[j].len;

			M(i, j) = dist_anchor(&E, sequences[j].S, ql, threshold);

#pragma omp atomic update
			progress_counter++;
		}

		if (print_progress) {
			size_t local_progress_counter;
			size_t num_comparisons = n * n - n;

#pragma omp atomic read
			local_progress_counter = progress_counter;

			double progress =
				100.0 * (double)local_progress_counter / num_comparisons;

#pragma omp critical
			fprintf(stderr, "\rComparing %zu sequences: %5.1f%% (%zu/%zu)", n,
					progress, local_progress_counter, num_comparisons);
		}

		esa_free(&E);
		seq_subject_free(&subject);
	}

	if (print_progress) {
		fprintf(stderr, ", done.\n");
	}
}
