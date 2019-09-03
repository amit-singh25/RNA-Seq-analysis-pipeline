#!/bin/bash
#PBS -N hisat
#PBS -j oe 
##PBS -l file=32GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8	
#PBS -o ./hisat_$PBS_JOBID.out
#PBS -e ./hisat_$PBS_JOBID.err
echo $PBS_JOBID
echo $PBS_JOBNAME
cd $PBS_O_WORKDIR

#echo "=========================================================="
#echo "Starting on : $(date)"
#echo "Running on node : $(hostname)"
#echo "Current directory : $(pwd)"
#echo "Current job ID : $JOB_ID"
#echo "Current job name : $JOB_NAME"
#echo "Task index number : $SGE_TASK_ID"
#echo "=========================================================="


name=$1
gtf=~/genome/hisat_genome/mouse
data=~/raw_data
out=~/alignment/hisat_r_alignemnt
#####creat a hisat index or download from igenome portal od hisat webpage

hisat2-build mouse_90.fa ${gtf}/mouse_90_hisat_index
#extract_splice_sites.py
##get exon 
#hisat2_extract_exons.py  ${gtf}/mouse.gtf > exon.txt
##get splice site
#hisat2_extract_splice_sites.py ${gtf}/mouse.gtf  > splicesites.txt


#hisat2 -p 8 -x ${gtf}/mouse_90_hisat_index --dta -1 ${data}/${name}_R1_001.fastq.gz -2 ${data}/${name}_R2_001.fastq.gz | samtools view -Sb - > ${out}/${name}.bam
hisat2-align-l -p 8 -x ${gtf}/mouse_90_hisat_index --dta -1 ${data}/${name}_R1_001.fastq.gz -2 ${data}/${name}_R2_001.fastq.gz | samtools view -Sb - > ${out}/${name}.bam
#hisat2 -p 8 -x ${gtf}/chaetomium_hisat_index --max-intronlen 2000 --dta --rna-strandness F -U ${data}/${name}.txt.gz -S ${out}/${name}.sam
samtools sort -@ 8 -o ${out}/${name}_sort.bam ${out}/${name}.bam

