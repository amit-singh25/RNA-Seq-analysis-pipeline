
#!/bin/bash
#PBS -N stingtie 
##PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=25G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -o ./stingtie_$PBS_JOBID.out
#PBS -e ./stingtie_$PBS_JOBID.err
#PBS -V 

#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
put=~/stringte
dir=~/alignment
gtf=~/igenome/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-current/Genes
genome=~/igenome/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta

stringtie -p 16 -o /home/bq_asingh/stringte/${name}.gtf \
-G ${gtf}/genes.gtf \
-l ${name} ${dir}/${name}/accepted_hits.bam \

samtools sort -@ 8 -o ${dir}/${name}/accepted_sort_hits.bam ${dir}/${name}/accepted_hits.bam 
