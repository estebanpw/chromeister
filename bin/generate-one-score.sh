#!/usr/bin/env bash
CSV=$1
TH=$2

if [ $# -ne 2 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <index.csv> <threshold>"
   echo ""
   exit -1
fi


BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# get first genome in list (they are sorted)
currgenome=$(tail -n +2 "$CSV" | head -1 | awk -F "," '{print $6}')

# fill array of chromosomes similarity

array=()
arraytosort=()
names=()
currsum=0
othergencounter=0
# for problems with chromo X and Y
highest=1 
lastprint=$currgenome
# For all lines

cat $CSV | tail -n +2 > $1.temp

while IFS= read -r i
do

	othergenome=$(echo "$i" | awk -F "," '{print $6}')
	
	
	if [ "$othergenome" != "$currgenome" ]; then
	
		# Sort the array with temporal values
		sorted=($(printf '%s\n' "${arraytosort[@]}"|sort))
		# accumulate sum until threshold is reached
		usedValues=1
		usedValuesNext=2
		first=$sorted[$usedValues]
		next=$sorted[$usedValuesNext]
		finalvalue=$first
		divisor=1
		currdiff=$(awk -v a="$first" -v b="$next" -v c="$TH" 'BEGIN {if(b-a > c) print 1; else print 2; }')
		while [ $currdiff -eq 1 -a $usedValuesNext -lt $highest ];
		do
			finalvalue=$(awk -v a="$finalvalue" -v b="$next" 'BEGIN {print a+b}')
			usedValues=`expr $usedValues + 1`
			usedValuesNext=`expr $usedValuesNext + 1`
			first=$sorted[$usedValues]
			next=$sorted[$usedValuesNext]
			divisor=$(awk -v a="$divisor" 'BEGIN {print a+0.1}')
			
		done

		# Once we are done, finalvalue has to be divided
		finalvalue=$(awk -v a="$finalvalue" -v b="$divisor" 'BEGIN {print a/b}')
		
		# array holds the results
		#array[$highest]=$(awk -v a="$currsum" -v b="$othergencounter" 'BEGIN {print a/b}')
		array[$highest]=$finalvalue
		
		highest=`expr $highest + 1`
		currsum=0
		names+=($currgenome)
		lastprint=$currgenome
		currgenome=$othergenome
		othergencounter=0
		unset arraytosort
	else
		getvalue=$(echo "$i" | awk -F "," '{print $8}')
		# Copy value to array 
		arraytosort[$othergencounter]=$getvalue
		#currsum=$(awk -v a="$currsum" -v b="$getvalue" 'BEGIN {print a=a+(1-b); exit}')
		othergencounter=`expr $othergencounter + 1`

	fi

done < "$1.temp"

# do the last!!!
#if [ "$lastprint" == "$currgenome" ]; then
# Sort the array with temporal values
sorted=($(printf '%s\n' "${arraytosort[@]}"|sort))
# accumulate sum until threshold is reached
usedValues=1
usedValuesNext=2
first=$sorted[$usedValues]
next=$sorted[$usedValuesNext]
finalvalue=$first
divisor=1
currdiff=$(awk -v a="$first" -v b="$next" -v c="$TH" 'BEGIN {if(b-a > c) print 1; else print 2; }')
while [ $currdiff -eq 1 -a $usedValuesNext -lt $highest ];
do
	finalvalue=$(awk -v a="$finalvalue" -v b="$next" 'BEGIN {print a+b}')
	usedValues=`expr $usedValues + 1`
	usedValuesNext=`expr $usedValuesNext + 1`
	first=$sorted[$usedValues]
	next=$sorted[$usedValuesNext]
	divisor=$(awk -v a="$divisor" 'BEGIN {print a+0.1}')
	
done

# Once we are done, finalvalue has to be divided
finalvalue=$(awk -v a="$finalvalue" -v b="$divisor" 'BEGIN {print a/b}')
array[$highest]=$finalvalue

highest=`expr $highest + 1`
currsum=0
names+=($currgenome)	
currgenome=$othergenome
othergencounter=0
#fi


highest=`expr $highest - 1`
rm $1.temp

tsum=0
rm $1.csv.inter
aux=0
for ((i = 1; i <= highest; i++)); do
	echo "${names[${aux}]} ${array[${i}]}" >> $1.inter
	aux=`expr $aux + 1`
	#val=${array[${i}]}
	#tsum=$(awk -v a="$tsum" -v b="$val" '{print a=a+b}')
done

#awk -v a="$tsum" b="$highest" '{print a/b}'

#sumfirst=$(awk -F "," 'BEGIN{suma=0}{suma = suma + $8}END{print suma}' "$CSV")
#echo "$sumfirst"



