#!/bin/bash

# Run TargetScan's base script (targetscan_60.pl) on all seed files.
#
# Authors: Fabian Schmich (fabian.schmich@bsse.ethz.ch) and Mason Victors (mason.victors@recursionpharma.com)
#

SEEDS_DIR="$1"
TS_UTRSEQ="$2"
OUT_DIR="$3"
TSCAN=$(which targetscan_70.pl)

for s in $SEEDS_DIR/*.txt
do
	bname=$(basename $s)
	perl $TSCAN $s $TS_UTRSEQ $OUT_DIR/tscan.$bname
done
