#!/bin/sh -f

test_fasta -l 10000 > a.fa
test_fasta -l 10000 > b.fa

# this is expected to trigger the nan warning
./src/andi -j a.fa b.fa
EXIT_VAL=$?

rm -f a.fa b.fa
exit $EXIT_VAL
