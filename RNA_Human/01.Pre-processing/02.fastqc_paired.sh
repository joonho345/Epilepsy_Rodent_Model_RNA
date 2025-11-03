#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2

source ${ConfigFile}

#####
DataPath=${FastqPath}
OutPath=${FastqcPath}
Data1=${FASTQ1}
Data2=${FASTQ2}

#####
make_dir ${OutPath}

#####
fastqc -o ${OutPath} ${DataPath}${Data1}
fastqc -o ${OutPath} ${DataPath}${Data2}

