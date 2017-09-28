#!/usr/bin/env bash
DIR=$1
EXT=$2
DIM=$3
KMER=$4

array=()
x=0

if [ $# != 4 ]; then
	echo "***ERROR*** Use: $0 genomesDirectory extension dim kmer"
	exit -1
fi

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

				
				#echo "$BINDIR/run_and_plot_chromeister.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT 30 10000"
				if [[ ! -f ${seqX}.$EXT-${seqY}.$EXT.mat.png ]]; then
					$BINDIR/run_and_plot_chromeister.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT $KMER $DIM
				fi
			
		fi
	done
done