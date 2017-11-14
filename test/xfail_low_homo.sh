#!/bin/sh -f

test_fasta -l 100000 > a.fa
test_fasta -l 100000 > b.fa
test_fasta -l 100 > both.fa

cat a.fa both.fa | awk -vRS='>' '{if($1 == "S0")print ">"$0 > "S0.fa"}'
cat b.fa both.fa | awk -vRS='>' '{if($1 == "S1")print ">"$0 > "S1.fa"}'

# this is expected to trigger the low homology warning
./src/andi -j S0.fa S1.fa
EXIT_VAL=$?

rm -f a.fa b.fa both.fa S0.fa S1.fa
exit $EXIT_VAL
