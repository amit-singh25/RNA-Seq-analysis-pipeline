#!/bin/bash
#PBS -N htseq
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./htseq_$PBS_JOBID.err           # stderr file
#PBS -o ./htseq_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
gtf=~/genome
out=~/alignment
data=~/count

htseq-count -f bam \
-r pos \
--type=exon --idattr=gene_id \
--stranded=reverse \
--mode=intersection-nonempty \
${input}/${name}_sort_accepted_hits.bam \
${gtf}/Mus_musculus.GRCm38.90.gtf >${output}/${name}.htseq_count.txt


#samtools sort -n ${input}/${name}/accepted_hits.bam ${input}/${name}/accepted_hits.sorted.bam
#samtools view ${input}/${name}/accepted_hits.sorted.bam | htseq-count --mode=intersection-nonempty \
#--stranded=reverse --type=exon --idattr=gene_id - ${index}/genes.gtf \
#> ${out}/${name}.htseq_count.txt

echo ${name}
