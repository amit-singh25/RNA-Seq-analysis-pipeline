#!/bin/bash
#PBS -N kallisto 
##PBS -j oe 
##PBS -l file=100GB
#PBS-l mem=25G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -o ./kallisto_$PBS_JOBID.out
#PBS -e ./kallisto_$PBS_JOBID.err
#PBS -V 

#echo $PBS_JOBNAME
#echo $PBS_JOBID
#cd $PBS_O_WORKDIR

#echo "=========================================================="
#echo "Starting on : $(date)"
#echo "Running on node : $(hostname)"
#echo "Current directory : $(pwd)"
#echo "Current job ID : $JOB_ID"
#echo "Current job name : $JOB_NAME"
#echo "Task index number : $SGE_TASK_ID"
#echo "=========================================================="

name=$1
#genome=~/genome/kallisto_genome/mouse_noncode
#genome=~/genome/kallisto_genome/mouse_noncode/gencode
genome=~/genome/kallisto_genome/mouse_noncode/gencode
out=~/mouse_data/publish_data/gencode_out
data=~/mouse_data/publish_data

#kallisto index -i ${genome}/mouse_index_ensembl ${genome}/Mus_musculus.GRCm38.ncrna.fa.gz
kallisto quant -i ${genome}/mouse_index_gencod --plaintext -t 16 ${data}/${name}_1.fastq.gz ${data}/${name}_2.fastq.gz -o ${out}/${name}

echo ${name}




##wget ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz
##wget ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
###kalisto pipe line 
##cat Mus_musculus.GRCm38.cdna.all.fa.gz Mus_musculus.GRCm38.ncrna.fa.gz > Mus_musculus.GRCm38.rna.fa.gz
##kallisto index -i mouse_index_kalisto Mus_musculus.GRCm38.rna.fa.gz
####quantification for paired end sequnce
#kallisto quant -i mouse_index_kalisto -o count -b 100 read1.fastq.gz read2.fastq.gz

####make all tsv file to one file here we fisrt cut first colum thed cut all samples that are in 5th cloumn tpm  
#paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30,35,40,45 > transcript_tpms_all_samples.tsv

#ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
#cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
#mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
#rm -f header.tsv
#####
##we also need the name of gene also symbol so here belwo is the command line 
#grep "^>" mouse_cdna_ncrna_91_single_cell.fa | sed -e "s/>//g;s/ /\t/g" | cut -f 1,4,7 >tex.txt 



