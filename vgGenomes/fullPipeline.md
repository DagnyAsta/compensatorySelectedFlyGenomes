# Genomic analysis of compensated genomes from the *vg1* mutation

This pipeline is adapted and built on Paul Knoops run through of his pooled data (https://github.com/PaulKnoops/episodicSequenceData)

The scripts are writen as they will be executed on the server *Drosophila*. Therefore all location of data and output are according to that. 

##QC and aligning

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

Indexing the reference genome was the next step. I desided to use the mapper bwa (bwa-0.7.17) and fetched the genome on flybase, 6.15 version. For indexing I used the script *bwaindexing.sh*

Example:
```
program=/programs/bwa-0.7.17/
referance=/home/dagny/fly/flyRefForCompensation/
output=/home/dagny/fly/flyRefForCompensation/bwa_index/

${program}bwa index -p ${output}bwaIndexed_flyBase.dmel-all-chromosome-r6.15 ${referance}flyBase.dmel-all-chromosome-r6.15.fasta
```

After indexing the trimmed fasta files could be aligned to the reference genome, using bwa. I used the script *bwalign.sh*

Example:

```
#for vestigial 2013
program=/programs/bwa-0.7.17/
input=/data3/Fly/compensationGenomeAnalysis/cutadaptTrimm/cutadaptTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/bwAlign/bwalignTrimmedVestigial1_2013/
referance=/home/dagny/fly/flyRefForCompensation/bwa_index/

for file in ${input}*R1_001\.fastq\.gz
do
file1=${file:${#input}}
echo ${file1}
file2=${file1/"R1_001.fastq.gz"/"R2_001.fastq.gz"}
echo ${file2}
echo ${file1:0:-16}
${program}bwa mem -t 4 -M ${referance}bwaIndexed_flyBase.dmel-all-chromosome-r6.15 ${input}${file1} ${input}${file2} > ${output}bwaAligned_${file1:0:-16}_pe.sam 
done
```

Some quality control paralell to converting the files from sam format to a bam and sorting them was the next step after the mapping. I used samtools (samtools-1.7) to do that in the script *samtoolsConvert.sh*

Example:

```
#vestigial 2013
program=/programs/samtools-1.7/
input=/data3/Fly/compensationGenomeAnalysis/bwAling/bwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/samtoolSamToBam/samtooledBwalignTrimmedVestigial1_2013/

for file in ${input}*.sam
do
fileSam=${file:${#input}}
fileBam=qFiltSamSort_${fileSam/".sam"/".bam"}
echo ${fileSam}
echo ${fileBam}
samtools view -b -S -q 20 -F 0x0004 ${input}${fileSam} | samtools sort -O bam -o ${output}${fileBam}
done
```

Ideally you remove duplicates (possable PCR amplification). First I indexed the reference genome with picard using the script *picardIndex.sh*

Example:

```
#! /bin/bash

program=/home/dagny/picardTools/picard.jar
index_dir=/home/dagny/fly/flyRefForCompensation/
ref_genome=/home/dagny/fly/flyRefForCompensation/flyBase.dmel-all-chromosome-r6.15.fasta

java -jar ${program} CreateSequenceDictionary R=${ref_genome} O=${index_dir}flyBase.dmel-all-chromosome-r6.15.dict

```


Using picard in the script *picardSortDupRemov.sh* I sorted the bam files as requered with picard tools and removed the possable dublicates in the reads.

Example:

```
#vestigial2013
program=/home/dagny/picardTools/
input=/data3/Fly/compensationGenomeAnalysis/samtooledMergedBwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/picardedSamtooledMergedBwalignTrimmedVestigial1_2013/
tmp=/home/dagny/tmp/

for file in ${input}*.bam
do
fileA=${file:${#input}}
fileB=${fileA:0:-14}
fileB=${fileB:39}
echo ${fileA}
echo ${fileB}
java -Xmx2g -Djava.io.tmpdir=${tmp} -jar ${program}picard.jar SortSam I=${input}${fileA} O=${output}picSort_${fileA} VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${tmp}
java -Xmx2g -jar ${program}picard.jar MarkDuplicates I=${output}picSort_${fileA} O=${output}dupremo_picSort_${fileA}  M=${output}/dupstat/${fileB}dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
rm ${output}picSort_${fileA}
done

```

To eliminate the noise in the data I decided to realign the reads around indels. I used a script from Katharine Pelletier, a PhD student at Ian Dworkin's lab at that time, to lean on.

```
#script modified from Katharine Pelletier


#!/bin/bash


#Path to input directory
#input=/data3/Fly/compensationGenomeAnalysis/picardReadGrouped/picReadGrouMergePicardedSamtooledBwalignTrimmed20140815_DNASeq_PEand20140826_DNASeq_PE/

#Path to output directory
#output=/data3/Fly/compensationGenomeAnalysis/gatkIndelRealign/

#Variable for reference genome (non-zipped)
#refGenome=/home/dagny/fly/flyRefForCompensation/flyBase.dmel-all-chromosome-r6.15.fasta

#Path to GATK
#program=/home/dagny/GATK/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/

#files=(${input}*_FVW_R1_*.bam)

#for file in ${files[@]}
#do

#fileName=${file:${#input}}
#echo ${fileName}
##echo ${fileName:0:-3}

#java -Xmx8g -jar ${program}GenomeAnalysisTK.jar -I ${input}${fileName} -R ${refGenome} -T RealignerTargetCreator -o ${output}${fileName:0:-3}intervals

#java -Xmx8g -jar ${program}GenomeAnalysisTK.jar -I ${input}${fileName} -R ${refGenome} -T IndelRealigner -targetIntervals ${output}${fileName:0:-3}intervals -o ${output}indelRealigned_${fileName}

#done

```


After removing duplicate I merged the bam files, in the case of the vestigal data from 2013 it was only two sequence samples from 2 lanes. However in the net and rho data set this is more complicated and needs to be taken with extra caution, not to merge wrong files. I used the script *samtoolsMerge.sh*

```
#to merge the two technical replicates from the two lanes

#vestigial 2013
program=/programs/samtools-1.7/
input=/data3/Fly/compensationGenomeAnalysis/samtooledBwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/samtooledMergedBwalignTrimmedVestigial1_2013/

for file in ${input}*L005_pe.bam
do
fileL005=${file:${#input}}
fileL006=${fileL005/"_L005_"/"_L006_"}
fileMerged=merged_${fileL005:0:-12}_pe.bam
echo ${fileL005}
echo ${fileL006}
echo ${fileMerged}
samtools merge ${output}${fileMerged} ${input}${fileL005} ${input}${fileL006} 
done

```

##Fst calculation and more
 







