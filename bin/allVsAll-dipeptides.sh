#!/usr/bin/env bash
DIR=$1
EXT=$2

array=()
x=0

if [ $# != 2 ]; then
	echo "***ERROR*** Use: $0 genomesDirectory extension"
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
	

	seqX=${array[$i]}
	printf ':,%s,nothing\n' "$seqX"

	$BINDIR/dipeptides -db $DIR/${seqX}.$EXT

		
	
done




