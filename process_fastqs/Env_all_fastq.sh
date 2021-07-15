#!/bin/bash

FQ_DIR="../fastq/"
ls ${FQ_DIR}*.fastq.gz | sed "s/\.\.\/fastq\///g" | sed "s/\.fastq\.gz//g">list.txt #get a list of all fastq's in directory
mkdir sam
mkdir bam

while read -r ONE; do
    read -r TWO

    bowtie2 -x ../seq/pNL4_3_rev_cs -1 ${FQ_DIR}/${ONE}.fastq.gz -2 ${FQ_DIR}/${TWO}.fastq.gz --fast-local -S ${ONE}.sam --rdg 100,3 --rfg 100,3
    java countDMScodons $ONE.sam 8474 8560
    samtools view -S -b $ONE.sam > bam/$ONE.bam
    rm $ONE.sam
    echo "Processed $ONE"

done < list.txt
