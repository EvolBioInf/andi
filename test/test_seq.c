#include <glib.h>
#include "global.h"
#include <stdio.h>
#include <string.h>
#include "sequence.h"

int FLAGS = F_NONE;

void test_seq_basic(){

	seq_t S;

	seq_init( &S, "ACGT", "name");

	g_assert_cmpstr(S.S, ==, "ACGT");
	g_assert_cmpstr(S.name, ==, "name");
	g_assert_cmpuint(S.len, ==, 4);

	seq_free( &S);
}

void test_seq_full(){

	seq_t S;
	seq_subject subject;

	seq_init( &S, "ACGTTGCA", "name");
	int check = seq_subject_init( &subject, &S);

	g_assert_cmpint(check, ==, 0);

	g_assert_cmpstr(subject.RS, ==, "TGCAACGT#ACGTTGCA");
	g_assert_cmpuint(subject.RSlen, ==, 8*2+1);
	g_assert( subject.gc == 0.5);

	seq_subject_free( &subject);
	seq_free( &S);
}

void test_seq_nonacgt(){
	seq_t S;
	seq_subject subject;

	seq_init( &S, "11ACGTNN7682394689NNTGCA11", "name");
	seq_subject_init( &subject, &S);

	g_assert_cmpstr(S.S, ==, "ACGTTGCA");
	g_assert_cmpuint(S.len, ==, 8 );
	g_assert( FLAGS & F_NON_ACGT);

	g_assert_cmpstr(subject.RS, ==, "TGCAACGT#ACGTTGCA");
	g_assert_cmpuint(subject.RSlen, ==, 8*2+1);
	g_assert( subject.gc == 0.5);

	seq_subject_free( &subject);
	seq_free( &S);

	FLAGS = F_NONE;

	seq_init( &S, "@ACGT_!0TGCA        ", "name");
	seq_subject_init( &subject, &S);

	g_assert_cmpstr(S.S, ==, "ACGT!TGCA");
	g_assert_cmpuint(S.len, ==, 9 );
	g_assert( FLAGS & F_NON_ACGT);

	g_assert_cmpstr(subject.RS, ==, "TGCA;ACGT#ACGT!TGCA");
	g_assert_cmpuint(subject.RSlen, ==, 9*2+1);

	seq_subject_free( &subject);
	seq_free( &S);

	FLAGS = F_NONE;

}

int main(int argc, char *argv[])
{
	g_test_init( &argc, &argv, NULL);
	g_test_add_func("/seq/basic", test_seq_basic);
	g_test_add_func("/seq/full", test_seq_full);
	g_test_add_func("/seq/non acgt", test_seq_nonacgt);

	return g_test_run();
}

