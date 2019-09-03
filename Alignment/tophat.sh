#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8	
#PBS -o tophat.stdout
#PBS -e tophat.stderr
#PBS -o ./tophat_$PBS_JOBID.out
#PBS -e ./tophat_$PBS_JOBID.err



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

index=~/genome/igenome/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index
gtf=~/genome/igenome/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-current/Genes
input~/mouse_data
out=~/alignment

tophat2 -o ${out}/${name} -p 8 -g 2 \
-G ${gtf}/genes.gtf  \
--transcriptome-index= ${index}/genome ${input}/${name}_R1.fastq.gz ${input}/${name}_R2.fastq.gz

echo $name



#stat -j $JOB_ID
#echo "=========================================================="
#echo "Finished on : $(date)"
#echo "=========================================================="


