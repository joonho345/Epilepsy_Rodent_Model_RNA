#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2

source ${ConfigFile}

#####
DataPath=${FastqPath}
OutPath=${FastqcPath}
Data0=${FASTQ0}

#####
make_dir ${OutPath}

#####
fastqc -o ${OutPath} ${DataPath}${Data0}

