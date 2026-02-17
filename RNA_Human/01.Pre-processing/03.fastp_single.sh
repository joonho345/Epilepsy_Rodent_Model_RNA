#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${FastqPath_unzip}
OutPath=${FastpPath}

#####
make_dir ${OutPath}

#####
fastp \
	--in1 ${DataPath}${FASTQ0_unzip} \
	--out1 ${OutPath}${FASTP0} \
	--html ${OutPath}${ID}.report.html \
	--json ${OutPath}${ID}.report.json \
	--length_required 25 \
	--low_complexity_filter \
	--overrepresentation_analysis -P 20
