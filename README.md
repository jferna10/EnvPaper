---
title: Code Related to Functional and Structural Segregation of Overlapping Helices in HIV-1
---

This github repo contains code related to the submitted paper "Functional and Structural Segregation of Overlapping Helices in HIV-1". The files deposited here are intended to make the analyses - as they were done at the time of writing the paper transparent. However due to things like files being renamed (e.g. GEO nomenclature led to renaming the fastqs that you can download at GSE179046 from the ones that were actually pulled off the machine), the compute environment, etc you can't just run this code and get the figures.  But it should be pretty close, and you shouldn't hesitate to contact us if you notice any issues.

Here is the overview of what you will find in this repository:


## Stats
Basic QC metrics for MiSeq run for the Env Deep Mutational Scanning Data. These are run-wide stats like demultipliexing stats.


## Reports
Basic QC metrics for each of the fastqs for the Env Deep Mutational Scanning Data.

## process_fastqs
Code used to generate codon and amino acid counts.

The bulk of the work is done by countDMS, a simple java program which attempts to count codons from each fastq entry after aligning to a reference (pNL4_3_rev_cs.fa which was built into a bowtie2 index).

The BAM/SAMs are not provided here but can be sent if needed.


## aa_tab
Amino acid counts generated from countDMS with simple number relabeling to make the coordinates readable.

## codon_tab
Codon counts generated from countDMS. These may be useful if you care about specific codons.





## Fernandes_GEO_seq_template.xlsx

This file should provide metadata mapping naming changes between GEO and filenames in this repo.
