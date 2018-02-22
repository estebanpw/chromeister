#!/usr/bin/env bash
CSV=$1

if [ $# -ne 1 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <index.csv>"
   echo ""
   exit -1
fi


BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# get first genome in list (they are sorted)
currgenome=$(tail -n +2 "$CSV" | head -1 | awk -F "," '{print $6}')

# fill array of chromosomes similarity

array=()
currsum=0
othergencounter=0
# for problems with chromo X and Y
highest=1 
# For all lines

cat $CSV | tail -n +2 > $1.temp

while IFS= read -r i
do

	othergenome=$(echo "$i" | awk -F "," '{print $6}')
	
	
	if [ "$othergenome" != "$currgenome" ]; then
	
		#echo "I have curr gen $currgenome and other gen $othergenome but switching. Currsum= $currsum"
		# Put value in array
		array[$highest]=$(awk -v a="$currsum" -v b="$othergencounter" 'BEGIN {print a/b}')
		#echo "array (should: $currsum) has on $highest -> ${array[${highest}]}"
		highest=`expr $highest + 1`
		currsum=0
		currgenome=$othergenome
		othergencounter=0
	else
		getvalue=$(echo "$i" | awk -F "," '{print $8}')
		currsum=$(awk -v a="$currsum" -v b="$getvalue" 'BEGIN {print a=a+(1-b); exit}')
		othergencounter=`expr $othergencounter + 1`

	fi

done < "$1.temp"

highest=`expr $highest - 1`
rm $1.temp

tsum=0
rm $1.csv.inter
for ((i = 1; i <= highest; i++)); do
	echo "${array[${i}]}" >> $1.csv.inter
	#val=${array[${i}]}
	#tsum=$(awk -v a="$tsum" -v b="$val" '{print a=a+b}')
done

#awk -v a="$tsum" b="$highest" '{print a/b}'

#sumfirst=$(awk -F "," 'BEGIN{suma=0}{suma = suma + $8}END{print suma}' "$CSV")
#echo "$sumfirst"



