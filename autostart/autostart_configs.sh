#!/bin/bash
mv ../config.m ../tmp_config.m
Confs=$(ls ./configs/)
for i in $Confs
do
    echo
    echo ========================================================
    echo Config number $i : | sed -e 's/config//g' -e 's/\.m//g'
    cp ./configs/$i "../config.m"
    echo   
    echo Calculating matrices 
    echo 
    ../numbers
    echo
    echo Main simulation:
    echo
    ../main
done
mv ../tmp_config.m ../config.m
