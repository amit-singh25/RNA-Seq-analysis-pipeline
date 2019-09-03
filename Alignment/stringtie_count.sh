#!/bin/bash
#PBS -N s_count
##PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=25G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -o ./s_count.out
#PBS -e ./s_count.err
#PBS -V 

#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
bam=~/alignment
out=~/stringte/ballgown
gtf=~/stringte
stringtie -e -B -p 16 -G ${gtf}/stringtie_merged.gtf -o ${out}/${name}.gtf ${bam}/${name}/accepted_hits.bam
