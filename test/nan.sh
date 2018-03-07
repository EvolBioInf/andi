#!/bin/bash -f

SEED=${RANDOM_SEED:-0}
SEED2=0
if test $SEED -ne 0; then
	SEED=$((SEED + 1))
	SEED2=$((SEED + 2))
fi


./test/test_fasta -s $SEED -l 10000 > a_nan.fa
./test/test_fasta -s $SEED2 -l 10000 > b_nan.fa

# this is expected to trigger the nan warning
./src/andi -j a_nan.fa b_nan.fa 2>&1 | grep 'nan'
EXIT_VAL=$?


if [[ EXIT_VAL -ge 1 ]]; then
	echo "Triggering nan failed" >&2
	grep '^>' a_nan.fa b_nan.fa
fi

rm -f a_nan.fa b_nan.fa
exit $EXIT_VAL
