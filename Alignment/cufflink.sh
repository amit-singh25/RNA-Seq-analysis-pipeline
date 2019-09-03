
#!/bin/bash
#PBS -N cufflink 
#PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   

name=$1
out=~/cuff_out_out
dir=~/alignemnt
gtf=~/genome/igenome/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-current/Genes
genome=~/genome/igenome/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta

cufflinks -p 8 -o $(out)/cufflink/${name} \
 -G ${gtf}/genes.gtf -b${genome}/genome.fa \
 -u --library-type fr-unstranded \
${dir}/${name}/accepted_hits.bam \

echo ${name}
