#!/bin/bash

libDir=lib
binDir=bin
common=../library

networkSize=$1
meanDegree=$2
SIII=$3
IR=$4
randomSeed=$5

name=N${networkSize}M${meanDegree}${SI_II}_${IR}-${randomSeed}

function debugBuild {
	g++ -std=c++17 -Wall -g -fsanitize=leak -I ${common} -I ${libDir}\
	    main-SIR.cpp\
        -o ${binDir}/${name}
}

function build {
	g++ -std=c++17 -O3 -flto -march=native -I ${common} -I ${libDir} -o ${binDir}/${name} \
		main-SIR.cpp
}

./bin/${name}.out ${networkSize} ${meanDegree} ${SIII} ${IR} ${randomSeed}
rm bin/${name}.out




