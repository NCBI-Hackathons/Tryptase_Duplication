#!/bin/bash

SAM_FOLDER=$1
REF=$2

for i in $SAM_FOLDER/*.sam; do

 echo $i
 samtools view -bS $i | samtools sort - > $i.sorted.bam
 samtools index $i.sorted.bam

done

samtools merge $SAM_FOLDER/merged.bam $SAM_FOLDER/*.bam

samtools mpileup -a -d 1000000 -f $REF $SAM_FOLDER/merged.bam > $SAM_FOLDER/merged.pileup.txt
