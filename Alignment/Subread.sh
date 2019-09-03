#!/bin/bash
#PBS -N subread 
##PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8  
#PBS -e ./subread_$PBS_JOBID.err           # stderr file
#PBS -o ./subread_$PBS_JOBID.out

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
data=~/mouse_data/publish_data/data
out=~/subread
index=~/genome/subread_genome/mouse
###generate subread index
subread-buildindex -o ${index}/mouse_index_90 ${index}/mouse_genome_90.fa 

subread-align -T 16 -t 0 -i ${index}/mouse_index_92 -r ${data}/${name}_R1_001.fastq.gz -R ${data}/${name}_R2_001.fastq.gz -o ${out}/${name}_accepted_hits.bam
samtools sort ${out}/${name}_accepted_hits.bam ${out}/${name}_sort_accepted_hits.bam
#samtools index ${out}/${name}_sort_accepted_hits.bam

echo $name
