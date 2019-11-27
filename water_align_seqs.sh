#!/bin/bash

aseq=$1
bseq=$2
print "Make sure that the the first input file is a single sequence.\n";
echo "My incoming sequences are $aseq and $bseq";
##Assuming that the file is always within a directory (otherwise the sed command does not work)

#aseq_n="$(echo $aseq| sed 's/.*\/\([a-z]*\).*/\1/')"
#bseq_n="$(echo $bseq| sed 's/.*\/\([a-z]*\).*/\1/')"

aseq_n="$(echo $aseq| sed 's/\([a-z]*\_[a-z]*\).*/\1/')"
bseq_n="$(echo $bseq| sed 's/\([a-z]*\_[a-z]*\).*/\1/')"

echo "My modified sequence nameis ${aseq_n}."

water -asequence $aseq -bsequence $bseq -gapopen 10 -gapextend 0.5 -outfile ${aseq_n}_vs_${bseq_n}.water
