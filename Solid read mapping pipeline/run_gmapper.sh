#!/bin/bash

# Define the pairs of SRR suffixes you want to process
declare -a pairs=("9596" "9798" "9900" "0102")

# Iterate over each pair
for pair in "${pairs[@]}"
do
    # Extract the first and second numbers for the pair
    first="${pair:0:2}"
    second="${pair:2:2}"

    # Determine if we need to use SRR12041 for the second part of the pair
    if [[ "$pair" == "9900" ]]; then
        # For the 9900 pair, use SRR12040 for 99 and SRR12041 for 00
        gmapper-cs -L referenceData/rnafasta/gencode.v46.transcripts-cs \
        -1 fastq/lastGremoved/processed/SRR12040${first}-p.fastq.gz \
        -2 fastq/lastGremoved/processed/SRR12041${second}-p.fastq.gz \
        -N 16 -p opp-in -Q \
        --strata --single-best-mapping --all-contigs \
        --min-avg-qv 20 \
        > samfiles/intermediate/${first}${second}_map.mini.sam 2> samfiles/intermediate/${first}${second}_map.mini.log

    elif [[ "$pair" == "0102" ]]; then
        # Use SRR12041 for both parts of the 0102 pair
        gmapper-cs -L referenceData/rnafasta/gencode.v46.transcripts-cs \
        -1 fastq/lastGremoved/processed/SRR12041${first}-p.fastq.gz \
        -2 fastq/lastGremoved/processed/SRR12041${second}-p.fastq.gz \
        -N 16 -p opp-in -Q \
        --strata --single-best-mapping --all-contigs \
        --min-avg-qv 20 \
        > samfiles/intermediate/${first}${second}_map.mini.sam 2> samfiles/intermediate/${first}${second}_map.mini.log
    else
        # Use SRR12040 for both parts in other cases
        gmapper-cs -L referenceData/rnafasta/gencode.v46.transcripts-cs \
        -1 fastq/lastGremoved/processed/SRR12040${first}-p.fastq.gz \
        -2 fastq/lastGremoved/processed/SRR12040${second}-p.fastq.gz \
        -N 16 -p opp-in -Q \
        --strata --single-best-mapping --all-contigs \
        --min-avg-qv 20 \
        > samfiles/intermediate/${first}${second}_map.mini.sam 2> samfiles/intermediate/${first}${second}_map.mini.log
    fi

    echo "Processing SRR pair ${first}${second} completed."
done
