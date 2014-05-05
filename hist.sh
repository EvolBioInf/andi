#!/bin/bash

echo "" > hist.log
./np -s hist $1
tail -n +2 hist.log | awk '{if(/[0-9]/){n[$1]++;s[$1]+=$2}}END{for(l in n){print l, n[l], s[l]/n[l]}}' | sort -n
