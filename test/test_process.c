#include "global.h"
#include "process.h"
#include <glib.h>
#include <math.h>

int FLAGS = 0;
int THREADS = 1;
long unsigned int BOOTSTRAP = 0;
double ANCHOR_P_VALUE = 0.025;
gsl_rng *RNG = NULL;
int MODEL = M_JC;

double shustring_cum_prob(size_t x, double g, size_t l);
size_t min_anchor_length(double p, double g, size_t l);

void test_shustring_cum_prob() {
	int len = 100000;
	double gc = 0.5;
	double p_value = 0.025;

	size_t threshold = min_anchor_length(p_value, gc, len);

	g_assert_cmpfloat(1 - p_value, <,
					  shustring_cum_prob(threshold + 1, gc / 2, len));
	g_assert_cmpfloat(1 - p_value, <=,
					  shustring_cum_prob(threshold, gc / 2, len));
	g_assert_cmpfloat(1 - p_value, >,
					  shustring_cum_prob(threshold - 1, gc / 2, len));
}

int main(int argc, char *argv[]) {
	g_test_init(&argc, &argv, NULL);
	g_test_add_func("/process/shustring_cum_prob", test_shustring_cum_prob);

	return g_test_run();
}
