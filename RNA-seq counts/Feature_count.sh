#!/bin/bash
#PBS -N f_count 
##PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=25G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -o ./f_count_$PBS_JOBID.out
#PBS -e ./f_count_$PBS_JOBID.err
#PBS -V 

#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
gtf=~/genome/subread_genome/mouse
out=~/subread
data=~/subread/count
featureCounts -p -s 2 -T 16 -a ${gtf}/Mus_musculus.GRCm38.90.gtf -t exon -g gene_id -o ${data}/${name}_f_counts.txt ${out}/${name}_sort_accepted_hits.bam

###clean the data that would take directly to the DEseq2

cut -f1,7-15 ${data}/${name}_f_counts.txt |sed 1d >${data}/${name}_counts_clean.txt
#cut -f 1,7-15 ${data}/${name}_f_counts.txt > ${data}/${name}_counts_clean.txt
#cat -n text.txt | sed '11d' 
echo ${name}
