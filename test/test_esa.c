#include "esa.h"
#include "global.h"
#include <glib.h>
#include <stdio.h>
#include <string.h>

int FLAGS = F_NONE;
int THREADS = 1;
double ANCHOR_P_VALUE = 0.025;

extern const int CACHE_LENGTH;

char code3char(ssize_t code) {
	switch (code & 0x7) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		case 4: return '!';
		case 5: return ';';
		case 6: return '#';
	}
	return '\0';
}

typedef struct {
	esa_s *C;
	seq_t *S;
	seq_subject subject;
} esa_fixture;

void assert_equal_lcp(const lcp_inter_t *a, const lcp_inter_t *b) {
	g_assert_cmpint(a->i, ==, b->i);
	g_assert_cmpint(a->j, ==, b->j);
	g_assert_cmpint(a->l, ==, b->l);
}

void assert_equal_cache_nocache(const esa_s *C, const char *str, size_t qlen) {
	lcp_inter_t a = get_match_cached(C, str, qlen);
	lcp_inter_t b = get_match(C, str, qlen);
	assert_equal_lcp(&a, &b);
	g_assert(strncmp(str, C->S + C->SA[a.i], a.l) == 0);
	g_assert(str[a.l] != C->S[a.l + C->SA[a.i]] || str[a.l] == '\0');
}

void setup(esa_fixture *ef, gconstpointer test_data) {
	ef->C = malloc(sizeof(esa_s));
	ef->S = malloc(sizeof(seq_t));

	g_assert(ef->C != NULL);
	g_assert(ef->S != NULL);

	const char *seq = {
		"TACGAGCACTGGTGGAATTGATGTC"
		"CAGTCTTATATGGCGCACCAGGCTG"
		"ATAGTAGTAGCAGTTTGCTTATCTC"
		"ATCGCGTGTTTCCGGATGACAGAGA"
		"TACGTGCACTGGTGGGATTGATGTC"
		"TAGTATTATATGGCGCACCAGGATG"
		"ATAGTAGTAGCAGTTTGCTTATCCC"
		"ATCGCGTGTTTGCGGATGACCGAGA" // format hack
	};

	g_assert(seq_init(ef->S, seq, "S0") == 0);
	seq_subject_init(&ef->subject, ef->S);
	g_assert(ef->subject.RS != NULL);
	int check = esa_init(ef->C, &ef->subject);
	g_assert(check == 0);
}

void setup2(esa_fixture *ef, gconstpointer test_data) {
	ef->C = malloc(sizeof(esa_s));
	ef->S = malloc(sizeof(seq_t));

	g_assert(ef->C != NULL);
	g_assert(ef->S != NULL);

	const char *seq = {
		"TACGAGCACTGGTGGAATTGATGTC"
		"CAGTCTTATATGGCGCACCAGGCTG"
		"ATAGTAGTAGCAGTTTGCTTATCTC"
		"ATCGCGTGTTTCCGGATGACAGAGA"
		"!"
		"TACGTGCACTGGTGGGATTGATGTC"
		"TAGTATTATATGGCGCACCAGGATG"
		"ATAGTAGTAGCAGTTTGCTTATCCC"
		"ATCGCGTGTTTGCGGATGACCGAGA" // format hack
	};

	g_assert(seq_init(ef->S, seq, "S0") == 0);
	seq_subject_init(&ef->subject, ef->S);
	g_assert(ef->subject.RS != NULL);
	int check = esa_init(ef->C, &ef->subject);
	g_assert(check == 0);
}

void teardown(esa_fixture *ef, gconstpointer test_data) {
	esa_free(ef->C);
	free(ef->C);
	seq_free(ef->S);
	free(ef->S);
	seq_subject_free(&ef->subject);
}

extern int count;

void basic(esa_fixture *ef, gconstpointer test_data) {
	esa_s *C = ef->C;
	g_assert(C->SA);

	lcp_inter_t a = get_match_cached(C, "AAGACTGG", 8);
	lcp_inter_t b = get_match(C, "AAGACTGG", 8);
	assert_equal_lcp(&a, &b);
	g_assert(strncmp("AAGACTGG", C->S + C->SA[a.i], 8) == 0);

	a = get_match_cached(C, "AATTAAAA", 8);
	b = get_match(C, "AATTAAAA", 8);
	assert_equal_lcp(&a, &b);
	g_assert(strncmp("AATTAAAA", C->S + C->SA[a.i], a.l) == 0);

	a = get_match_cached(C, "ACCGAGAA", 8);
	b = get_match(C, "ACCGAGAA", 8);
	assert_equal_lcp(&a, &b);
	g_assert(strncmp("ACCGAGAA", C->S + C->SA[a.i], a.l) == 0);

	a = get_match_cached(C, "AAAAAAAAAAAA", 12);
	b = get_match(C, "AAAAAAAAAAAA", 12);
	assert_equal_lcp(&a, &b);
	g_assert(strncmp("AAAAAAAAAAAA", C->S + C->SA[a.i], a.l) == 0);

	// g_assert_cmpint(count, >=, 1 << (2*8));
}

void normq_cached(esa_fixture *ef, gconstpointer test_data) {
	esa_s *C = ef->C;
	g_assert(C->SA);
	lcp_inter_t a, b;

	a = get_match_cached(C, "A", 1);
	b = get_match(C, "A", 1);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "C", 1);
	b = get_match(C, "C", 1);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "CT", 2);
	b = get_match(C, "CT", 2);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "AAGACTGG", 8);
	b = get_match(C, "AAGACTGG", 8);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "AATTAAAA", 8);
	b = get_match(C, "AATTAAAA", 8);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "ACCGAGAA", 8);
	b = get_match(C, "ACCGAGAA", 8);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "AAAAAAAAAAAA", 12);
	b = get_match(C, "AAAAAAAAAAAA", 12);
	assert_equal_lcp(&a, &b);

	a = get_match_cached(C, "!AAAAAAAAAAAA", 12);
	b = get_match(C, "!AAAAAAAAAAAA", 12);
	assert_equal_lcp(&a, &b);
}

size_t MAX_DEPTH = 11;

void prefix_dfs(esa_s *C, char *str, size_t depth);

void prefix(esa_fixture *ef, gconstpointer test_data) {
	esa_s *C = ef->C;
	char str[MAX_DEPTH + 1];
	str[MAX_DEPTH] = '\0';
	prefix_dfs(C, str, 0);
}

void prefix_dfs(esa_s *C, char *str, size_t depth) {
	if (depth < MAX_DEPTH) {
		for (int code = 0; code < 4; ++code) {
			str[depth] = code2char(code);
			prefix_dfs(C, str, depth + 1);
		}
	} else {
		assert_equal_cache_nocache(C, str, depth);
	}
}

int main(int argc, char *argv[]) {
	g_test_init(&argc, &argv, NULL);
	g_test_add("/esa/basic", esa_fixture, NULL, setup, basic, teardown);
	g_test_add("/esa/sample cache", esa_fixture, NULL, setup, normq_cached,
			   teardown);
	g_test_add("/esa/sample cache 2", esa_fixture, NULL, setup2, normq_cached,
			   teardown);
	g_test_add("/esa/full cache", esa_fixture, NULL, setup, prefix, teardown);
	g_test_add("/esa/full cache 2", esa_fixture, NULL, setup2, prefix,
			   teardown);

	return g_test_run();
}
