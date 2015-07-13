#!/bin/bash

#usage . ./launchAnalysis directory maxProc
#directory is the directory to look for files, with the structure specified below, and maxProc is the max number of processes to execute at the same time (should be equal to te number of cpus)
#read the root files in a directory structure like base/sensor/bias where the bias folder contains the root files and a .cfg
#the base directory should be given as argument
#the root files are expected to have the .root extension

maxjobs=$2
jobCount=0 #number of daughter jobs

for iSens in $(ls -d $1/*/)
do
    #echo $iSens
    for iBias in $(ls -d $iSens/*/)
    do
#echo $iBias
	for iRun in $(ls -d $iBias/*.root)
	do
	    iConf=$(ls -d $iBias/*.cfg)
	    #xterm -e "cd $PWD; pwd; ./executables/chargeDistr $iRun $iConf" &
	    echo $iRun
	    ./executables/chargeDistr $iRun $iConf & # to not have new windows opening
	    let jobCount+=1 #increment job coutner
#echo "======================================================================bla $iRun $iConf"
	    if [ $jobCount -eq $maxjobs ] #if the number of submitted jobs is equal to the maximum, wait for the last of them to finish
	    then
# echo waiting
# echo $jobCount
		#wait $! #wait for the last submitted job to finish
		wait #wait for the submitted jobs to finish
		jobCount=0
	    fi
	done
    done
done

wait # wait for all jobs to finish
