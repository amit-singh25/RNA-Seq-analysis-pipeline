#!/bin/bash
#PBS -N cuff_diff
#PBS -j oe 

##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   

name=$1

output=~/cuff_dif_out
gtf=~/genome/igenome/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-current/Genes

cuffdiff -p 20 --min-reps-for-js-test 2 --dispersion-method per-condition --output-dir ${output}/cuffdiff 
 --library-type fr-firststrand --use-sample-sheet ${gtf}/genes.gtf sample_sheet_All.txt 
