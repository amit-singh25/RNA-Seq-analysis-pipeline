#!/bin/bash
#PBS -N s_merge
##PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=25G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -o ./s_merge.out
#PBS -e ./s_merge.err
#PBS -V 

#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
put=~/stringte
dir=~/alignment
gtf=~/genome/igenome/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-current/Genes
genome=~/genome/igenome/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta
stringtie --merge -p 16 -G ${gtf}/genes.gtf -o ${put}/stringtie_merged.gtf ${put}/mergelist.txt
