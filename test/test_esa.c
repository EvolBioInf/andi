#include <glib.h>
#include "global.h"
#include <stdio.h>
#include <string.h>
#include "esa.h"
#include <stdlib.h>

int FLAGS = F_NONE;
int THREADS = 1;

extern const int CACHE_LENGTH;

char code3char( ssize_t code){
	switch( code & 0x7){
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
	esa_t *C;
	seq_t *S;
} esa_fixture;

void assert_equal_lcp( const lcp_inter_t *a, const lcp_inter_t *b){
	g_assert_cmpint( a->i, ==, b->i);
	g_assert_cmpint( a->j, ==, b->j);
	g_assert_cmpint( a->l, ==, b->l);
}

void assert_equal_cache_nocache( const esa_t *C, const char *str, size_t qlen){
	lcp_inter_t a = get_match_cached(C, str, qlen);
	lcp_inter_t b = get_match(C, str, qlen);
	assert_equal_lcp( &a, &b);
}

void test_esa_setup( esa_fixture *ef, gconstpointer test_data){
	ef->C = (esa_t *) malloc( sizeof(esa_t));
	ef->S = (seq_t *) malloc( sizeof(seq_t));

	g_assert( ef->C != NULL);
	g_assert( ef->S != NULL);

	const char *seq = {
		"TACGAGCACTGGTGGAATTGATGTC"
		"CAGTCTTATATGGCGCACCAGGCTG"
		"ATAGTAGTAGCAGTTTGCTTATCTC"
		"ATCGCGTGTTTCCGGATGACAGAGA"
		"TACGTGCACTGGTGGGATTGATGTC"
		"TAGTATTATATGGCGCACCAGGATG"
		"ATAGTAGTAGCAGTTTGCTTATCCC"
		"ATCGCGTGTTTGCGGATGACCGAGA"
	};

	g_assert( seq_init( ef->S, seq, "S0" ) == 0);
	seq_subject_init( ef->S);
	g_assert( ef->S->RS != NULL);
	int check = esa_init( ef->C, ef->S);
	g_assert( check == 0);
}

void test_esa_setup2( esa_fixture *ef, gconstpointer test_data){
	ef->C = (esa_t *) malloc( sizeof(esa_t));
	ef->S = (seq_t *) malloc( sizeof(seq_t));

	g_assert( ef->C != NULL);
	g_assert( ef->S != NULL);

	const char *seq = {
		"TACGAGCACTGGTGGAATTGATGTC"
		"CAGTCTTATATGGCGCACCAGGCTG"
		"ATAGTAGTAGCAGTTTGCTTATCTC"
		"ATCGCGTGTTTCCGGATGACAGAGA"
		"!"
		"TACGTGCACTGGTGGGATTGATGTC"
		"TAGTATTATATGGCGCACCAGGATG"
		"ATAGTAGTAGCAGTTTGCTTATCCC"
		"ATCGCGTGTTTGCGGATGACCGAGA"
	};

	g_assert( seq_init( ef->S, seq, "S0" ) == 0);
	seq_subject_init( ef->S);
	g_assert( ef->S->RS != NULL);
	int check = esa_init( ef->C, ef->S);
	g_assert( check == 0);
}

void test_esa_teardown( esa_fixture *ef, gconstpointer test_data){
	esa_free(ef->C);
	free(ef->C);
	seq_free(ef->S);
	free(ef->S);
}

extern int count;

void test_esa_basic( esa_fixture *ef, gconstpointer test_data){
	esa_t *C = ef->C;
	g_assert( C->SA);

	lcp_inter_t a = get_match_cached(C, "AAGACTGG", 8);
	lcp_inter_t b = get_match(C, "AAGACTGG", 8);
	assert_equal_lcp( &a, &b);

	a = get_match_cached(C, "AATTAAAA", 8);
	b = get_match(C, "AATTAAAA", 8);
	assert_equal_lcp( &a, &b);

	a = get_match_cached(C, "ACCGAGAA", 8);
	b = get_match(C, "ACCGAGAA", 8);
	assert_equal_lcp( &a, &b);

	a = get_match_cached(C, "AAAAAAAAAAAA", 12);
	b = get_match(C, "AAAAAAAAAAAA", 12);
	assert_equal_lcp( &a, &b);

	//g_assert_cmpint(count, >=, 1 << (2*8));
}


size_t MAX_DEPTH = 9;

void test_esa_prefix_dfs( esa_t *C, char *str, size_t depth);

void test_esa_prefix( esa_fixture *ef, gconstpointer test_data){
	esa_t *C = ef->C;
	char str[MAX_DEPTH+1];
	str[MAX_DEPTH] = '\0';
	test_esa_prefix_dfs( C, str, 0);
}

void test_esa_prefix_dfs( esa_t *C, char *str, size_t depth){
	if( depth < MAX_DEPTH){
		for( int code = 0; code < 4; ++code){
			str[depth] = code2char(code);
			test_esa_prefix_dfs( C, str, depth + 1);
		}
	} else {
		assert_equal_cache_nocache(C, str, depth);
	}
}

int main(int argc, char *argv[])
{
	g_test_init( &argc, &argv, NULL);
	g_test_add("/esa/basic", esa_fixture, NULL, test_esa_setup, test_esa_basic, test_esa_teardown);
	g_test_add("/esa/cache", esa_fixture, NULL, test_esa_setup, test_esa_prefix, test_esa_teardown);

	return g_test_run();
}

