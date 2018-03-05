#!/bin/bash -f

SEED=${RANDOM_SEED:-0}
SEED2=0
if test $SEED -ne 0; then
	SEED=$((SEED + 1))
	SEED2=$((SEED + 2))
fi


./test/test_fasta -s $SEED -l 10000 > a.fa
./test/test_fasta -s $SEED2 -l 10000 > b.fa

# this is expected to trigger the nan warning
./src/andi -j a.fa b.fa 2>&1 | grep 'nan'
EXIT_VAL=$?


if [[ EXIT_VAL -ge 1 ]]; then
	echo "Triggering nan failed" >&2
	grep '^>' a.fa b.fa both.fa
fi

rm -f a.fa b.fa
exit $EXIT_VAL
