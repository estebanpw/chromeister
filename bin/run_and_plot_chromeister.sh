#!/usr/bin/env bash
G1=$1
G2=$2
DIM=$3


FILE1=$(basename $G1)
FILE2=$(basename $G2)

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


if [ $# -lt 3 ]; then
        echo "***ERROR*** Use: $0 <G1> <G2> <DIMENSION>"
        exit -1
fi



(time $BINDIR/CHROMEISTER -query $G1 -db $G2 -kmer 32 -out $FILE1-$FILE2.mat -dimension $DIM) &> $FILE1-$FILE2.log
Rscript $BINDIR/compute_score.R $FILE1-$FILE2.mat $DIM

#$BINDIR/plot.R $FILE1-$FILE2.mat


#rm $FILE1-$FILE2.mat
