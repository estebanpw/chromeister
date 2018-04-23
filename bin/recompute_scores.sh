#!/usr/bin/env bash

if [ $# -ne 1 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <#threads>"
   echo ""
   exit -1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
THREADS=$1


for p in */;
do

	echo "Entering $p"
	cd $p



	thepaths=()
	n=0
	# Grab the routes
	for i in *.mat ;
	do

	        thepaths[$n]=$i
        	n=`expr $n + 1`

	done

	i=0
	aux=0
	for ((i=0; i < $n ; i+=$THREADS))
	do
		aux=$i
		for ((j=0; j<$THREADS; j++))
		do
			if [ "$aux" -lt "$n" ]; then
				echo "Recomputing ${thepaths[$aux]}"
				goodpath=${thepaths[$aux]%.mat}
				Rscript $BINDIR/compute_score.R ${thepaths[$aux]} > ${goodpath}.scr.txt &
				aux=`expr $aux + 1`
			fi
		done

		for job in `jobs -p`
		do
		    #echo $job
		    wait $job
		done
		
	done


	echo "Exiting $p"
	cd ..


done


