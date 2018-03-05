#!/bin/bash -f

SEED=${RANDOM_SEED:-0}
SEED2=0
SEED3=0
if test $SEED -ne 0; then
	SEED=$((SEED + 1))
	SEED2=$((SEED + 2))
	SEED3=$((SEED + 3))
fi

./test/test_fasta -s $SEED -l 100000 > a.fa
./test/test_fasta -s $SEED2 -l 100000 > b.fa
./test/test_fasta -s $SEED3 -l 100 > both.fa

cat both.fa a.fa | awk -vRS='>' '{if($1 == "S0")print ">"$0 > "S0.fa"}'
cat both.fa b.fa | awk -vRS='>' '{if($1 == "S1")print ">"$0 > "S1.fa"}'

# this is expected to trigger the low homology warning
./src/andi -j S0.fa S1.fa 2>&1 | grep 'homology'
EXIT_VAL=$?

if [[ EXIT_VAL -ge 1 ]]; then
	echo "Triggering low homology failed" >&2
	grep '^>' a.fa b.fa both.fa
fi

rm -f a.fa b.fa both.fa S0.fa S1.fa
exit $EXIT_VAL
