# Genomic analysis of the compensated genomes 

This pipeline is adapted and built on Paul Knoops run through of his pooled data (https://github.com/PaulKnoops/episodicSequenceData)

The scripts are writen as they will be executed on the server *Drosophila*. Therefore all location of data and output are according to that. 

### Quality control and rimming using FastQC and cutadapt

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
### Aligning with bwa followed by a quality control with samtools and picard tools

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

In case I did an extra step of quality control after the merging with the script *samtoolsQC.sh*

```
#second QC if the dup removeal affected the quality scale

program=/programs/samtools-1.7/
input=/data3/Fly/compensationGenomeAnalysis/picardDupRem/picardedSamtooledMergedBwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/finalBam/vestigial/
for file in ${input}*.bam

do

fileBam=${file:${#input}}
echo ${fileBam}
samtools view -q 20 -F 0x0004 -b ${input}${fileBam} > ${output}secQfilt_${fileBam}
rm ${input}${fileBam}

done
```

YOU NEED TO ADD THE INDEXING OF THE BAM FILES


### Pi and Fst calculation with Popoolation including a step of indel realignment with GATK
 
In order to prepare the data for downstream analysis with popoolation I used samtools to pile the bam files up using the script *samtoolsMpileup.sh*

Example: 

```
#to pileup the bamfiles in order to make syncronized files for downstream popoolation analysis

####vestigial 2013
program=/programs/samtools-1.7/
input=/data3/Fly/compensationGenomeAnalysis/gatkIndelRealign/indelRealignPicReadGrouPicardedSamtooledMergedBwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/samtoolPileUp/pileUpIndelRealignPicardedSamtooledMergedBwalignTrimmedVestigial1_2013/
ref=/home/dagny/fly/flyRefForCompensation/bwa_index/
#creat one pileup file for one bam file - per population
for file in ${input}*.bam
do
fileA=${file:${#input}}
fileB=${fileA:0:-14}
fileB=vestigial2013${fileB:85}
echo ${fileA}
echo ${fileB}
${program}samtools mpileup -B -f ${ref}flyBase.dmel-all-chromosome-r6.15.fasta ${input}${fileA} > ${output}chr2L_${fileB}.mpileup
done

```
note: in the same script there are multiple options of pileup combination, pileing all togeather, only CAS and NASC *etc.*

Downstream I used Popoolation2 (version popoolation2_1201) and Popoolation (version popoolation_1.2.2) to calculate the allel frequency difference. 

For somewhat of a QC I did calculate the pi for the population with the popoolation script *Variance-sliding.pl* executed within the bash script *poopolationTajimaPiCalculation.sh*

Example:


```
#! /bin/bash    

#this script is to execute the perl script from popoolations which calculates tajima pi 

program=/home/dagny/popoolationTools/popoolation_1.2.2/
input=/data3/Fly/compensationGenomeAnalysis/samtoolPileUp/pileUpPicardedSamtooledMergedBwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/poopolationPiCalc/vestigial

for file in vestigial2013CAS_R1.mpileup vestigial2013CAS_R2.mpileup vestigial2013CAS_R3.mpileup vestigial2013FVW_WT.mpileup vestigial2013NASC_R1.mpileup vestigial2013NASC_R2.mpileup vestigial2013NASC_R3.mpileup vestigial2013VG_BASE.mpileup

do

fileBase=${file:0:-8}
echo ${file}
echo ${fileBase}

perl ${program}Variance-sliding.pl --input ${input}${file} --output ${output}${fileBase}.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output ${output}${fileBase}.snps --min-covered-fraction 0.5

done
```

In order to calculate the allel difference between populations I piled the bam files all togeather in one large mpileup file with samtools using the sam script above with another option (*samtoolsMpileup.sh*)

Example: 

```
#to pileup the bamfiles in order to make syncronized files for downstream popoolation analysis

####vestigial 2013
program=/programs/samtools-1.7/
input=/data3/Fly/compensationGenomeAnalysis/gatkIndelRealign/indelRealignPicReadGrouPicardedSamtooledMergedBwalignTrimmedVestigial1_2013/
output=/data3/Fly/compensationGenomeAnalysis/samtoolPileUp/pileUpIndelRealignPicardedSamtooledMergedBwalignTrimmedVestigial1_2013/
ref=/home/dagny/fly/flyRefForCompensation/bwa_index/


create a pileupfile for all the populations
${program}samtools mpileup -B -f ${ref}flyBase.dmel-all-chromosome-r6.15.fasta ${input}indelRealigned*.bam > ${output}vestigial2013AllPop.mpileup

```

Then I syncroniced the mpileup file using the popoolation2 script *mpileup2sync.jar* excecuted in the bash script *poopolationSyncFiles.sh*

Example:

```
#script that will use the jar file from poopolation to sync the piled up bam files from the vg data

#vestigial
#program=/home/dagny/popoolationTools/popoolation2_1201/
#input=/data3/Fly/compensationGenomeAnalysis/samtoolPileUp/pileUpIndelRealignPicardedSamtooledMergedBwalignTrimmedVestigial1_2013/
#output=/data3/Fly/compensationGenomeAnalysis/poopolationSync/

#java -ea -Xmx4g -jar ${program}mpileup2sync.jar --input ${input}vestigial2013CAS123andNASC123RG2BiologicalRepMerged.mpileup --output ${output}vestigial2013CAS123andNASC123RGBiologicalRepMerged.sync --fastq-type sanger --min-qual 20 --threads 2
```

With the syncroniced file I could calculate fst values for windows or bases. This I did with the perl popoolation2 script *fst-sliding.pl* executed in the bash script *poopolationFstCalculation.sh*

Example:

```
#this script is to calculate fst from sync files (from mpileupfiles, synced with popoolation2), with popoolation2


#vestigial calculation
#program=/home/dagny/popoolationTools/popoolation2_1201/
#input=/data3/Fly/compensationGenomeAnalysis/poopolationSync/
#output=/data3/Fly/compensationGenomeAnalysis/poopolationFstCalc/vestigial2013AllPop/

#perl ${program}fst-sliding.pl --input ${input}*.sync --output ${output}vestigial2013AllPop_50to200Cov1window.fst --min-count 6 --min-coverage 50 -max-coverage 200 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 120
```

Note: Within the same bash script there are multiple options of parameter that I did try. We did decide on that window sice of 1000 and step sice of 500 would be reasonable. The coverage was harder to establish as much of the X chromosome did get wiped out if the coverage was set to low. I did try to set the max coverage to see if that was causing the issue but it did not to much of a difference, therefore I figured the min coverage was the determing factor of this problem. I ended up seting the min coverage to 30 and max coverage around 200. The min coverage fraction I naturally set to 1 when calculating base pare wise but for windows of sice 1000 to 0.6 (this is the proportion within the frame you are calculating that needs to meet the min coverage requirement).

### Preparing and ploting Fst values in R

When preparing the fst files and plotting the values across the genome I used the R script *fstPlot.R*. First ofcourse the fst files need to be transfered from the server drosophila to your computer, this I do with a simple bash script you see below:

```
scp -r -P 24680 dagny@nf-ux148.rhi.hi.is:/data3/Fly/compensationGenomeAnalysis/poopolationFstCalc/vestigial2013AllPopBiologicalRepMergedRepMasked /Users/Dagny/Dropbox/CompensationDNAseq/data/poopolationFstCalc/afterIndelRealign/
```

Note: The R script I will not give an example of as it is long to give an example of. There is an example of a manhattan plot below. 



Takk og bless.







