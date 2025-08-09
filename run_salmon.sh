#!/bin/bash

FASTQ_DIR=$1
RESULT_DIR=$2

SALMON_INDEX="salmon_index"

samples=$(ls $FASTQ_DIR/*_1.fastq.gz | xargs -n1 basename | sed 's/_1.fastq.gz//' | sort | uniq)

for sample in $samples; do

    mkdir -p $RESULT_DIR/"$sample"

    fq1="$FASTQ_DIR/${sample}_1.fastq.gz"
    fq2="$FASTQ_DIR/${sample}_2.fastq.gz"

    salmon quant \
        -i $SALMON_INDEX \
        -l A \
        -1 $fq1 \
        -2 $fq2 \
        -o $RESULT_DIR/"$sample" \
        --meta \
        --minScoreFraction 0.25 \
        --consensusSlack 0.2 \
        --maxRecoverReadOcc 1000

done
