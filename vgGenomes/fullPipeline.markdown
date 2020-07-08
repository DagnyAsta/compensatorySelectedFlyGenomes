# Genomic analysis of compensated genomes from the *vg1* mutation

```sh
#this is a script for building an indexed drosophila genome from flybase to align on, both with bowtie2 version 2.3.4.1
input=/home/dagny/fly/flyRefForCompensation/
output=/home/dagny/fly/flyRefForCompensation/bowtie2_index/
program=/programs/bowtie2/
${program}bowtie2-build ${input}flyBase.dmel-all-chromosome-r6.15.fasta flyBase.dmel-all-chromosome-r6.15
```
