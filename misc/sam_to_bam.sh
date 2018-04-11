#!/bin/bash

SAM_FOLDER=$1

for i in $SAM_FOLDER/*.sam; do

 echo $i
 samtools view -bS $i | samtools sort - > $i.sorted.bam
 samtools index $i.sorted.bam

done
