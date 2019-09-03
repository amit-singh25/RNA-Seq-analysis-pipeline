#!/bin/bash
#PBS -N run_star 
#PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8	


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

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ~/genome/star_genome/mouse \
--genomeFastaFiles ~/genome/star_genome/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile ~/genome/star_genome/mouse/Mus_musculus.GRCm38.90.gtf \


#echo $name
#qstat -j $JOB_ID

#echo "=========================================================="
#echo "Finished on : $(date)"
#echo "=========================================================="
