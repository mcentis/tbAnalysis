#!/bin/bash

#usage . ./modifyCfgOneSensor.sh directory find replace

#uses sed to modify the config files in a directory structure like sensor/bias
#the base directory should be given as argument
#the config files are expected to have the .cfg extension


for iBias in $(ls -d $1/*)
do
    #echo $iBias
    iConf=$(ls -d $iBias/*.cfg)
    if [ ! -z $iConf ] # exclude the case of files not found
    then
	#echo "======================================== $iConf"
	sed -i "s|$2|$3|g" $iConf
    fi
done
