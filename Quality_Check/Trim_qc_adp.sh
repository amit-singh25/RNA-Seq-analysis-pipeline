#!/bin/bash
#PBS -N Trim_qc_adp
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./Trim_qc_adp_$PBS_JOBID.err           # stderr file
#PBS -o ./Trim_qc_adp_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
data=~/raw_data

flexbar -r ${name}_R1.fastq.gz -p ${name}_R2.fastq.gz -t ${data}/${name}_trimmed.fastq  -n 20 -z GZ -m 30 -u 0  -q TAIL -qt 28 -a ${data}/adaptor.fa -qf sanger -j


###find adapator from the fasta_file or else take a illumina adaptor sequence
##bbmerge.sh in1=r1.fq in2=r2.fq outa=adapters.fa strict
###
https://github.com/optimuscoprime/autoadapt
