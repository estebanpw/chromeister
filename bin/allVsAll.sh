#!/usr/bin/env bash
DIR=$1
EXT=$2
DIM=$3
KMER=$4
DIFF=$5

array=()
x=0

if [ $# -lt 5 ]; then
	echo "***ERROR*** Use: $0 genomesDirectory extension dim kmer diff [optional: grid]"
	exit -1
fi

indexnameA=$(basename "$DIR")

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for elem in $(ls -d $DIR/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	array[$x]=$elem
	x=`expr $x + 1`
	#echo "X: $elem"
done

for ((i=0 ; i < ${#array[@]} ; i++))
do
	for ((j=i ; j < ${#array[@]} ; j++))
	do
		if [ $i != $j ]; then
				seqX=${array[$i]}
				seqY=${array[$j]}
				echo "----------${seqX}-${seqY}-----------"

				
				if [[ ! -f ${seqX}.$EXT-${seqY}.$EXT.mat ]]; then
					echo "Running $BINDIR/run_and_plot_chromeister.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT $KMER $DIM $DIFF $6"
					$BINDIR/run_and_plot_chromeister.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT $KMER $DIM $DIFF $6
					#Rscript $BINDIR/compute_score.R $seqX.$EXT-$seqY.$EXT.mat $DIM > $seqX.$EXT-$seqY.$EXT.scr.txt
				fi
			
		fi
	done
done

# generate index
if [[ ! -f index.csv.temp ]] && [ ! -f index-$indexnameA.csv  ]; then
	
	echo "Launching... $BINDIR/index_chromeister_solo.sh index-$indexnameA.csv $DIR $DIR"
	$BINDIR/index_chromeister_solo.sh . index-$indexnameA.csv $DIR $DIR
fi




