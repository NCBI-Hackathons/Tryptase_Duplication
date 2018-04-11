#!/bin/bash

for i in *.bam; do

  echo $i
  samtools view -h $i Alpha_GEX_64k_HEX Alpha_GEX_79k_dup_FAM > $i.filtered.sam
  samtools view -bS $i.filtered.sam > $i.filtered.bam
  samtools index $i.filtered.bam
  echo "$i finished"


done
