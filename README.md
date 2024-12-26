SOLiD reads mapping pipeline
================

## Intro

This document describes how to align reads for Solid color space reads.
Basics about difference between letter space reads (ATCG) generated by
illumina platforms and color spaced reads (0101) can be googled. It
should be noted that solid reads seem to be not popular in recent years
and the mapping platforms are mostly discontinued. Though attractive, it
is not proper to transform the color spaced reads into letter spaced
reads for downstream analysis (reasons can be found in this link and
more online:
<https://davetang.org/muse/2012/07/04/from-sra-to-fastq-for-solid-data/>)

This file uses sra files from this source:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1356833>

### Download sra files with sratoolkit

``` bash
for sra_id in SRR1204091 SRR1204092 SRR1204093 SRR1204094 SRR1204095 SRR1204096 SRR1204097 SRR1204099 SRR1204098 SRR1204100 SRR1204101 SRR1204102 SRR1204103 SRR1204104 SRR1204105 SRR1204106 SRR1204107 SRR1204108 SRR1204109 SRR1204110 SRR1204111 SRR1204112 SRR1204113 SRR1204114 SRR1204115 SRR1204116 SRR1204117 SRR1204118 SRR1204119 SRR1204120 SRR1204121 SRR1204122 SRR1204123 SRR1204124 SRR1204125 SRR1204126 SRR1204127 SRR1204128 SRR1204129 SRR1204130 SRR1204131 SRR1204132; do
    prefetch $sra_id
done

```

### Fastqdump into fastq files

Theoretically the solid reads should not use fastqdump, but abi-dump
instead. But I checked the results, the reads are the same, just there
are mysterious C within the reads post fastqdump for the solid color
space reads. So I removed all the Cs in the fastq files with C script
(use ‘preprocess_fastq.c’).

``` bash
for sraFile in sraDir/sra/SRR*; do
    echo "Extracting fastq from ${sraFile}"
    fastq-dump \
        --origfmt \
        --gzip \
        --outdir sraDir/fastq \
        ${sraFile}
done
```

Alternatively, it can be down with: abidump and solid2fastq pearl
scripts. The final fastq results are the same.

Anyway, fastq files were generated from the .qual and .csfasta files of
the solid reads for subsequent reads alignment.

### Reads mapping with Shrimp2.

Refer to shrimp2 documentation for more details about the command:
<https://github.com/compbio-UofT/shrimp/blob/master/README> The command
I used:

``` bash
time gmapper-cs -L referenceData/rnafasta/gencode.v46.transcripts-cs \
-1 fastq/lastGremoved/processed/SRR1204091-p.fastq.gz \
  -2 fastq/lastGremoved/processed/SRR1204092-p.fastq.gz \
-N 16 -p opp-in -Q  \
--strata --single-best-mapping --all-contigs \
--min-avg-qv 20 \
>samfiles/9192_map.mini.sam 2>samfiles/9192_map.mini.log
```

Note, the paired ends reads were mapped to total reference genome, not
whole genome (hg19), as whole genome needs to be split into 2 files due
to RAM limits, and it takes roughly 48hours for one sample. When mapping
into reference genome, comparable read count results were generated
(check reports from ‘compare_transcriptome_mapped.md’ ), and it only
take half of the total time (\<24hours)

The projection of the reference genome was generated to allow efficient
load of the reference genome before aligning:

``` bash
#for example:
#cs for color space mode
utils/project-db.py --shrimp-mode cs referencegenome.fasta

```

### Samtools to sort and filter the paired reads only

``` bash

#filter paired end reads only, and then sort
samtools view -@ 12 -b -f 2 9394_map.mini.sam | \
samtools sort -@ 12 -o 9394_sorted.bam && \
samtools sort -n -@ 12 -o 9394_paired_sorted_mini.bam 9394_sorted.bam

```

### Export read count

Since the reads were mapped to transcriptome reference, read count can
be directly counted and exported:

``` bash
#export to read count
samtools view -@ 4 9394_paired_sorted_mini.bam | awk '{print $3}' | cut -d'|' -f2,6 | sort | uniq -c > DEG/9192_mini_gene_id_name_counts.txt
```
