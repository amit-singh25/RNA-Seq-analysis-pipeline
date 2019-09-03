#!/bin/bash
#PBS -N htseq
#PBS -j oe 

##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   

name=$1



gtf=~/genome/igenome/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-current/Genes
out=~/cuff_diff_out

cuffmerge -g ${gtf}/genes.gtf <assembly_GTF_list.txt> -o ${out}
