#!/usr/bin/env bash
# Use in folder prior to all folders!

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for i in *; do

	
	$BINDIR/generate-one-score.sh $i/index-$i.csv 0.2
        value=$(awk '{sum+=$2} END { print sum/NR}' $i/index-$i.csv.inter)
        echo "$i-$value"


done

