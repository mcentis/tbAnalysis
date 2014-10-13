#!/bin/bash

#usage . ./subsituteFiles directory dierctoryWithNewFiles
#directory is the directory to look for files, with the structure specified below, and dierctoryWithNewFiles is the directory that contains the new files (usually when the reconstruction is runned again)
#the directory structure should be like base/sensor/bias where the bias folder contains the root files and a .cfg
#the base directory should be given as argument
#the root files are expected to have the .root extension

unsDir=$2 #directory with unsorted files
#jobCount=0 #number of daughter jobs

for iSens in $(ls -d $1/*/)
do
    #echo $iSens
    for iBias in $(ls -d $iSens/*/)
    do
#echo $iBias
	for iRun in $(ls $iBias/*.root)
	do
	    runNum=${iRun##*/} # extract run number (what comes after the last /)
	    runNum=${runNum%%-*} # extract run number (what comes before the first -)
	    echo $runNum

	    newFile=$unsDir/$runNum-alibava-tracking-3.root # line with fixed name!!!!!!!!!
	    echo $newFile
	    if [ -e $newFile ]
	    then
	    	rm $iRun
		mv $newFile $iBias
	    fi
	done
    done
done

