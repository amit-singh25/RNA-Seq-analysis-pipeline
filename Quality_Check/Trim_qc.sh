#!/bin/bash
#PBS -N fastqc
#PBS -j oe 
##PBS -l file=16GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   

name=$1

data=~/mouse_data
out=~/fastqc
fastqc ${data}${name}.fastq.gz > ${out}${name}
echo ${name}

