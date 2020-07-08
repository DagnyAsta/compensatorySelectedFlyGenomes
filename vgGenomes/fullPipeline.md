# Genomic analysis of compensated genomes from the *vg1* mutation

The scripts are writen as they will be executed on the server *Drosophila*. Therefore all location of data and output are according to that. 

First of all I did fastQC (FastQC v0.11.7) on the fasta files using the script *fastqc.sh*

Example:
```
#This is a bash script for running FastQC on the genomes from the compensationary experiment

#vestigial 2013 data
input=/data3/Fly/compensation_genomes/vestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/fastqcAnalysis/vestigial1_2013/

fastqc -o ${output} --noextract ${input}*.gz 

```

Next I did clean the reads with cutadapt (version 1.16) using the script *cutadapt.sh*:
Example:
```
#cutadapt script to clean the sequenses 

#vg samples from 2013
program=/usr/bin/
input=/data3/Fly/compensation_genomes/vestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/cutadaptTrimm/cutadaptTrimmedVestigial1_2013/

for file in ${input}*R1_001\.fastq\.gz
do
file1=${file:48}
echo ${file1}
file2=${file1/"R1_001.fastq.gz"/"R2_001.fastq.gz"}
echo ${file2}
echo ${file:48:(-16)}
${program}cutadapt -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 15 -o ${output}trimmed_${file1} -p {output}trimmed_${file2} ${input}${file1} ${input}${file2} > ${output}/summaries/summary_${file:48:(-16)}.txt
done





  ```
  #this is a script for building an indexed drosophila genome from flybase to align on, both with bowtie2 version 2.3.4.1
  
  input=/home/dagny/fly/flyRefForCompensation/
  output=/home/dagny/fly/flyRefForCompensation/bowtie2_index/
  program=/programs/bowtie2/
  
  ${program}bowtie2-build ${input}flyBase.dmel-all-chromosome-r6.15.fasta flyBase.dmel-all-chromosome-r6.15
  
  ```
