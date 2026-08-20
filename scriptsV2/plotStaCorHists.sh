#!/bin/bash

# plot the full station corrections for inspection
# requires outputs from previous plot scripts

ls recdata
test -f recdata || exit

ls tmpx
test -f tmpx || exit
grep RES tmpx > mydata.txt

test -f synStaCors.dat && echo "also plotting synStaCors.dat"
# histByGroup written by Claude:
histByGroup.sh mydata.txt recdata

