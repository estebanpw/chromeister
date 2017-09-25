#!/usr/bin/env bash
G1=$1
G2=$2
KMER=$3
DIM=$4


FILE1=$(basename $G1)
FILE2=$(basename $G2)

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


if [ $# -lt 4 ]; then
        echo "***ERROR*** Use: $0 <G1> <G2> <KMER SIZE> <DIMENSION> <--filter/nothing>"
        exit -1
fi



(time $BINDIR/CHROMEISTER -query $G1 -db $G2 -kmer $KMER -out $FILE1-$FILE2.mat -dimension $DIM $5) &> $FILE1-$FILE2.log

$BINDIR/plot.R $FILE1-$FILE2.mat


#rm $FILE1-$FILE2.mat
