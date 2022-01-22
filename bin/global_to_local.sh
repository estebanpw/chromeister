FILECSV=$1
EVENTS=$2

tail -n +3 $FILECSV | awk -F "," 'BEGIN{nSeqs=1;} /\#/{exit(0)} /[0-9]+/{lengths[nSeqs]=$3; nSeqs++;} END{for(j=1; j<nSeqs; j++) print j","lengths[j]; } ' > $FILECSV.ytmp


cat $FILECSV | awk -F "," 'BEGIN{nSeqs=1; ACTIVE=0;} /\#/{ACTIVE=1} /[0-9]+/{if(ACTIVE==1){  lengths[nSeqs]=$3; nSeqs++;}} END{for(j=1; j<nSeqs; j++) print j","lengths[j]; } ' > $FILECSV.xtmp

paste -d , $FILECSV.ytmp $FILECSV.xtmp > $FILECSV.yxtmp

rm $FILECSV.xtmp $FILECSV.ytmp

tail -n +2 $EVENTS > $EVENTS.tmp

#xStart yStart xEnd yEnd
awk -F "," 'BEGIN{sy=1; sx=1; print "xStart,yStart,xEnd,yEnd,strand,approximate length,description";} NR==FNR{ly[sy]=$2;  lx[sx]=$4; sy++; sx++;} 
NR!=FNR{
xStart=$1; yStart=$2; xEnd=$3; yEnd=$4;

lacX=0;  i=1; prevLacX=0; while(xStart > lacX && i<sx){ prevLacX=lacX; lacX+=lx[i]; i++;} newXstart=xStart-prevLacX;
lacX=0;  i=1; prevLacX=0; while(yStart > lacX && i<sy){ prevLacX=lacX; lacX+=ly[i]; i++;} newYstart=yStart-prevLacX;
lacX=0;  i=1; prevLacX=0; while(xEnd   > lacX && i<sx){ prevLacX=lacX; lacX+=lx[i]; i++;} newXend  =xEnd  -prevLacX;
lacX=0;  i=1; prevLacX=0; while(yEnd   > lacX && i<sy){ prevLacX=lacX; lacX+=ly[i]; i++;} newYend  =yEnd  -prevLacX;



print newXstart","newYstart","newXend","newYend","$5","$6","$7;


}' $FILECSV.yxtmp $EVENTS.tmp

rm $EVENTS.tmp $FILECSV.yxtmp
