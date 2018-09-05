#!/usr/bin/env bash
DIR=$1
FASTAS1=$3
FASTAS2=$4
OUT=$2

echo "Computing index CSV..." > index.csv.temp

while [ "$(find . -size 0 | wc -l)" -ne 0 ]; do
        sleep 10s
done


EXT="mat"
EXTSCORE="scr.txt"
EXTGENERAL=".fasta"

# data is like
# HOMSA.Chr.10.fasta



if [ $# != 4 ]; then
	echo "***ERROR*** Use: $0 <directory> <out> <fastas_directory_1> <fastas_directory_2>"
	exit -1
fi


rm $OUT

for elem in $(ls -d $DIR/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do

	IFS='-', read -a splits <<< "$elem"   # yields MUSMU.Chr.8.fasta and  MUSMU.Chr.Y.fasta
	IFS='.', read -a getnum <<< "$elem"   # yields MUSMU  Chr  8  fasta  and MUSMU Chr Y fasta

	scorepath=$(basename $elem .mat).$EXTSCORE

	sed -i "/X.*/d" $DIR/$scorepath
	sed -i "s/\[1\]//g" $DIR/$scorepath

	score=$( head -1  $DIR/$scorepath)
	len1=$( head -1 $DIR/$elem.mat.events.txt | awk -F="," '{print $1}')
	len2=$( head -1 $DIR/$elem.mat.events.txt | awk -F="," '{print $2}')

	

	file1=${splits[0]}
	file2=$(basename ${splits[1]} .mat)

	ID1=$(head -1 $FASTAS1/$file1)
	ID2=$(head -1 $FASTAS2/$file2)

	
	#scorepath="$(basename $elem .mat).$EXTSCORE"
	#score="$(head -1 $DIR/$scorepath)"

	counter=0
	numX=${getnum[2]}
	numY=${getnum[5]}




	echo "$(basename ${splits[0]} $EXTGENERAL),$(basename ${splits[1]} ${EXTGENERAL}),$ID1,$ID2,$elem.$EXT.filt.png,$numX,$numY,$score,$len1 $len2" >> $OUT

done

sort -k5,5n -k6,6n -o $OUT $OUT 


sed -i '1iSpX, SpY, IDX, IDY, IMG, CHNumberX, CHNumberY, Score, LengthX, LengthY' $OUT

rm index.csv.temp

