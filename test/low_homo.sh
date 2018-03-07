#!/bin/bash -f

SEED=${RANDOM_SEED:-0}
SEED2=0
SEED3=0
if test $SEED -ne 0; then
	SEED=$((SEED + 1))
	SEED2=$((SEED + 2))
	SEED3=$((SEED + 3))
fi

./test/test_fasta -s $SEED -l 100000 > a_low.fa
./test/test_fasta -s $SEED2 -l 100000 > b_low.fa
./test/test_fasta -s $SEED3 -l 100 > both_low.fa

cat both_low.fa a_low.fa | awk -vRS='>' '{if($1 == "S0")print ">"$0 > "S0_low.fa"}'
cat both_low.fa b_low.fa | awk -vRS='>' '{if($1 == "S1")print ">"$0 > "S1_low.fa"}'

# this is expected to trigger the low homology warning
./src/andi -j S0_low.fa S1_low.fa 2>&1 | grep 'homology'
EXIT_VAL=$?

if [[ EXIT_VAL -ge 1 ]]; then
	echo "Triggering low homology failed" >&2
	grep '^>' a_low.fa b_low.fa both_low.fa
fi

rm -f a_low.fa b_low.fa both_low.fa S0_low.fa S1_low.fa
exit $EXIT_VAL
