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
homologies=()
condition=0
othergencounter=0
# for problems with chromo X and Y
highest=1 
# For all lines

cat $CSV | tail -n +2 > $1.temp

while IFS= read -r i
do

	othergenome=$(echo "$i" | awk -F "," '{print $6}')
	if [ "$condition" -eq 0 ]; then
		currgenome=$othergenome
		condition=1
	fi
	
	if [ "$othergenome" != "$currgenome" ]; then
	
		# Sort the array with temporal values
		#printf '%s\n' "${arraytosort[@]}"
		#echo "name is $currgenome"
		
		sorted=($(printf '%s\n' "${arraytosort[@]}"|sort))
		
		#echo "For chroomo $currgenome we have "
		#echo $(printf '%s,' "${sorted[@]}")
		# accumulate sum until threshold is reached
		usedValues=1
		usedValuesNext=2
		first=${sorted[0]}
		next=${sorted[${usedValues}]}
		nextofnext=${sorted[${usedValuesNext}]}
		finalvalue=$first
		divisor=0
		currdiff=$(LC_NUMERIC=POSIX awk -v a="$next" -v b="$nextofnext" 'BEGIN {print b-a }')
		TH=$(printf '%4.6f' $TH)
		#echo "$(LC_NUMERIC=POSIX awk -v a="$currdiff" -v b="$TH" 'BEGIN { printf("comp %f > %f = %d",a,b,a>b)} ')"
		condition=$(LC_NUMERIC=POSIX awk -v a="$currdiff" -v b="$TH" 'BEGIN { printf("%d",a>b)} ')
		#echo "first $first next $next result $currdiff condition $condition divisor $divisor th $TH finalvalue $finalvalue"
		while [ $condition -eq 1 -a $usedValuesNext -lt ${#sorted[@]} ];
		do
			usedValues=`expr $usedValues + 1`
			usedValuesNext=`expr $usedValuesNext + 1`
			finalvalue=$(LC_NUMERIC=POSIX awk -v a="$finalvalue" -v b="$next" 'BEGIN {print (a+b)}')
			next=${sorted[${usedValues}]}
			nextofnext=${sorted[${usedValuesNext}]}
			
			currdiff=$(LC_NUMERIC=POSIX awk -v a="$nextofnext" -v b="$next" 'BEGIN {printf("%f", b-a) }')
			condition=$(LC_NUMERIC=POSIX awk -v a="$currdiff" -v b="$TH" 'BEGIN { printf("%d", a>b)} ')
			divisor=$(LC_NUMERIC=POSIX awk -v a="$divisor" 'BEGIN {print a+0.1}')
			
		done

		#echo "so this is what we got $finalvalue, when divided using $divisor"
		
		# array holds the results
		#array[$highest]=$(awk -v a="$currsum" -v b="$othergencounter" 'BEGIN {print a/b}')
		finalvalue=$(LC_NUMERIC=POSIX awk -v a="$finalvalue" -v b="$usedValues" -v c="$divisor" 'BEGIN {printf("%f", a/(b-c))}')
		array[$highest]=$finalvalue
		homologies[$highest]=$usedValues
		
		highest=`expr $highest + 1`
		condition=0
		names+=($currgenome)
		othergencounter=0
		unset arraytosort
		
		
		
		
		getvalue=$(echo "$i" | awk -F "," '{print $8}')
		# Copy value to array 
		arraytosort[$othergencounter]=$getvalue
		#currsum=$(awk -v a="$currsum" -v b="$getvalue" 'BEGIN {print a=a+(1-b); exit}')
		othergencounter=`expr $othergencounter + 1`
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

usedValues=1
usedValuesNext=2
first=${sorted[0]}
next=${sorted[${usedValues}]}
nextofnext=${sorted[${usedValuesNext}]}
finalvalue=$first
divisor=0
currdiff=$(LC_NUMERIC=POSIX awk -v a="$next" -v b="$nextofnext" 'BEGIN {print b-a }')
TH=$(printf '%4.6f' $TH)
#echo "$(LC_NUMERIC=POSIX awk -v a="$currdiff" -v b="$TH" 'BEGIN { printf("comp %f > %f = %d",a,b,a>b)} ')"
condition=$(LC_NUMERIC=POSIX awk -v a="$currdiff" -v b="$TH" 'BEGIN { printf("%d",a>b)} ')
#echo "first $first next $next result $currdiff condition $condition divisor $divisor th $TH finalvalue $finalvalue"
while [ $condition -eq 1 -a $usedValuesNext -lt ${#sorted[@]} ];
do
        usedValues=`expr $usedValues + 1`
        usedValuesNext=`expr $usedValuesNext + 1`
        finalvalue=$(LC_NUMERIC=POSIX awk -v a="$finalvalue" -v b="$next" 'BEGIN {print (a+b)}')
        next=${sorted[${usedValues}]}
        nextofnext=${sorted[${usedValuesNext}]}

        currdiff=$(LC_NUMERIC=POSIX awk -v a="$nextofnext" -v b="$next" 'BEGIN {printf("%f", b-a) }')
        condition=$(LC_NUMERIC=POSIX awk -v a="$currdiff" -v b="$TH" 'BEGIN { printf("%d", a>b)} ')
        divisor=$(LC_NUMERIC=POSIX awk -v a="$divisor" 'BEGIN {print a+0.1}')

done

#echo "so this is what we got $finalvalue, when divided using $divisor"

# array holds the results
#array[$highest]=$(awk -v a="$currsum" -v b="$othergencounter" 'BEGIN {print a/b}')
finalvalue=$(LC_NUMERIC=POSIX awk -v a="$finalvalue" -v b="$usedValues" -v c="$divisor" 'BEGIN {printf("%f", a/(b-c))}')
array[$highest]=$finalvalue
homologies[$highest]=$usedValues


highest=`expr $highest + 1`
currsum=0
names+=($currgenome)	
currgenome=$othergenome
othergencounter=0
#fi


highest=`expr $highest - 1`
rm $1.temp

tsum=0
echo "deleteme" > $1.inter
rm $1.inter
aux=0
for ((i = 1; i <= highest; i++)); do
	echo "${names[${aux}]} ${array[${i}]} ${homologies[${i}]}" >> $1.inter
	#echo "${names[${aux}]} ${array[${i}]} ${homologies[${i}]}"
	aux=`expr $aux + 1`
	#val=${array[${i}]}
	#tsum=$(awk -v a="$tsum" -v b="$val" '{print a=a+b}')
done

#awk -v a="$tsum" b="$highest" '{print a/b}'

#sumfirst=$(awk -F "," 'BEGIN{suma=0}{suma = suma + $8}END{print suma}' "$CSV")
#echo "$sumfirst"



