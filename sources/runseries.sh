#!/bin/bash

echo "ENTER SERIES NUMBER FOR ANALYSIS: "
read SERIES_NO

DIR_NAME="S_0$SERIES_NO"
LOG_NAME="logfile_S0$SERIES_NO.txt"

cd ../build/
./data        $SERIES_NO | tee -a $LOG_NAME
./attenuation $SERIES_NO | tee -a $LOG_NAME
./timeres     $SERIES_NO | tee -a $LOG_NAME
./tconst      $SERIES_NO | tee -a $LOG_NAME 
./lightout    $SERIES_NO | tee -a $LOG_NAME

if
   cd ../results/$DIR_NAME ; then
   echo "RESULTS DIRECTORY EXISTS!"
else
   mkdir ../results/$DIR_NAME
   cd ../results/$DIR_NAME
   echo "CREATED AND ENTERED RESULTS DIRECTORY!"
fi

mv ../../build/$LOG_NAME ./
mv ../data_series$SERIES_NO.root ./
mv ../attenuation_series$SERIES_NO.root ./
mv ../timingres_series$SERIES_NO.root ./
mv ../timeconst_series$SERIES_NO.root ./
mv ../lighoutput_series_$SERIES_NO.root ./

DATE=$(date +%d-%m-%Y" "%H:%M:%S)
echo "" >> $LOG_NAME
echo "ANALYSIS DONE AT: $DATE" >> $LOG_NAME