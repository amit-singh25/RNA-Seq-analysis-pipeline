#!/bin/bash
#PBS -N Star_alignment
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./Star_alignment_$PBS_JOBID.err           # stderr file
#PBS -o ./Star_alignment_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID



name=$1

gtf=~/genome/star_genome/mouse
data=~/raw_data
out=~/alignment
STAR --runThreadN 8 --genomeDir ${gtf} --sjdbGTFfile ${gtf}/Mus_musculus.GRCm38.90.gtf \
--readFilesCommand zcat \
--readFilesIn ${data}/${name}_1.fastq.gz \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${out}/${name} \
