#!/bin/bash

ls *.fastq >list.txt #get a list of all fastq's in directory

while read -r ONE; do
    read -r TWO
    FNAME=$(echo "$ONE" | cut -d. -f1)
    
    bowtie2 -x ~/hive/jferna10/EnvLibraries/Mapping/rev_cs_bowtie2_index -1 $ONE -2 $TWO --fast-local -S $FNAME.sam --rdg 100,3 --rfg 100,3
    java countDMScodons $FNAME.sam 8474 8560
    echo "Processed $FNAME"

    #  	bowtie2 -x ~/hive/jferna10/EnvLibraries/Mapping/rev_cs_bowtie2_index -1 LLP2_256_a_S4_L001_R1_001.fastq.gz -2 LLP2_256_a_S4_L001_R2_001.fastq.gz --fast-local -S tmp.sam --rdg 100,3 --rfg 100,3
done < list.txt



