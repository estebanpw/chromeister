#!/usr/bin/env bash
G1=$1
G2=$2
KMER=$3
DIM=$4
DIFF=$5

if [ $# -lt 5 ]; then
        echo "***ERROR*** Use: $0 <G1> <G2> <KMER> <DIMENSION> <DIFF> [optional: grid]"
        exit -1
fi


FILE1=$(basename $G1)
FILE2=$(basename $G2)

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ $6 == "grid" ]]; then

	(time $BINDIR/CHROMEISTER -query $G1 -db $G2 -kmer $KMER -out $FILE1-$FILE2.mat -dimension $DIM -diffuse $DIFF) &> $FILE1-$FILE2.log
	(Rscript $BINDIR/compute_score.R $FILE1-$FILE2.mat $DIM) &> $FILE1-$FILE2.scr.txt

else

	(time $BINDIR/CHROMEISTER -query $G1 -db $G2 -kmer $KMER -out $FILE1-$FILE2.mat -dimension $DIM -diffuse $DIFF) &> $FILE1-$FILE2.log
	(Rscript $BINDIR/compute_score-nogrid.R $FILE1-$FILE2.mat $DIM) &> $FILE1-$FILE2.scr.txt

fi

source $BINDIR/../chromeisterenv/bin/activate
python3 $BINDIR/detect_events.py $FILE1-$FILE2.mat.raw.txt
deactivate


#rm $FILE1-$FILE2.mat
