#!/bin/sh -f

# This scripts test the accuracy of andi with random inputs. For that
# it uses the small program test_random to generate pairs of sequences
# with a given distance. By default, test_random creates a new set of
# sequences each time it is called. Thus, this test has a small, but
# non-zero probability of failing. That is a problem with Debian's
# reproducible builds effort. So this script acts as a wrapper around
# this issue.
#
# Simply calling this script via
#     % ./test/test_random.sh
# checks a new test-case every time. But with the right parameter
#     % RANDOM_SEED=1729 ./test/test_random.sh
# one specific set of sequences is validated.

./src/andi --help > /dev/null || exit 1

LENGTH=100000

# If RANDOM_SEED is set, use its value. Otherwise 0 is used to signal
# to test_random that a new set of sequences shall be generated.
SEED=${RANDOM_SEED:-0}

for dist in 0.0 0.001 0.01 0.02 0.05 0.1 0.2 0.3
do
	for n in $(seq 10)
	do
		if test $SEED -ne 0; then
			SEED=$((SEED + 1))
		fi

		res=$(./test/test_fasta -s $SEED -l $LENGTH -d $dist |
			tee ./test/test_random.fasta |
			./src/andi -r -t 1 |
			tail -n 1 |
			awk -v dist=$dist '{print $2, dist}' |
			awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1-$2) <= 0.02 && abs($1-$2) <= 0.02 * $2}')
		if test $res -ne 1; then
			echo "The last test computed a distance deviating more than two percent from its intended value."
			echo "See test_random.fasta for the used sequences."
			echo "RANDOM_SEED=$RANDOM_SEED"
			head -n 1 ./test/test_random.fasta
			exit 1;
		fi
	done
done

rm ./test/test_random.fasta
