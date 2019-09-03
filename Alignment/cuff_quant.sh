#!/bin/bash
#PBS -N cuff_quant
#PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=20G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=16   
#PBS -o cuff_quant.stdout
#PBS -e cuff_quant.stderr

name=$1
put=~/stringte
dir=~/alignment
gtf=~/genome/igenome/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-current/Genes
genome=~/genome/igenome/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta

cuffquant -p 10 -b ${genome}/genome.fa -o ${put}/ballgown/${name}/ --library-type=fr-unstranded ${put}/stringtie_merged.gtf ${dir}/${name}/accepted_hits.bam
echo ${name}
