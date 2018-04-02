#!/usr/bin/env bash
DIR=$1
DIR2=$2
EXT=$3
DIM=$4
KMER=$5
TH=$6



if [ $# != 6 ]; then
	echo "***ERROR*** Use: $0 genomesDirectory1 genomesDirectory2 extension dim kmer th-phylo"
	exit -1
fi

indexnameA=$(basename "$DIR")
indexnameB=$(basename "$DIR2")

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir $indexnameA-$indexnameB
cd $indexnameA-$indexnameB

array=()
x=0
array2=()

for elem in $(ls -d $DIR/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	array[$x]=$elem
	x=`expr $x + 1`
	#echo "X: $elem"
done

x=0

for elem in $(ls -d $DIR2/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	array2[$x]=$elem
	x=`expr $x + 1`
	#echo "X: $elem"
done

for ((i=0 ; i < ${#array[@]} ; i++))
do
	for ((j=0 ; j < ${#array2[@]} ; j++))
	do
				seqX=${array[$i]}
				seqY=${array2[$j]}
				echo "----------${seqX}-${seqY}-----------"

				
				#echo "$BINDIR/run_and_plot_chromeister.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT 30 10000"
				if [[ ! -f ${seqX}.$EXT-${seqY}.$EXT.mat ]]; then
					
					$BINDIR/run_and_plot_chromeister.sh $DIR/${seqX}.$EXT $DIR2/${seqY}.$EXT $KMER $DIM
					Rscript $BINDIR/compute_score.R $seqX.$EXT-$seqY.$EXT.mat > $seqX.$EXT-$seqY.$EXT.scr.txt
				fi
			
	done
done

# generate index
if [[ ! -f index.csv.temp ]] && [ ! -f index-$indexnameA-$indexnameB.csv  ]; then
	
	echo "Launching... $BINDIR/index_chromeister.sh index-$indexnameA-$indexnameB.csv $DIR $DIR2"
	$BINDIR/index_chromeister.sh . index-$indexnameA-$indexnameB.csv $DIR $DIR2
	echo "Launching... $BINDIR/generate-one-score.sh index-$indexnameA-$indexnameB.csv"
	$BINDIR/generate-one-score.sh index-$indexnameA-$indexnameB.csv $TH
fi



