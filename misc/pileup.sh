#!/bin/bash

FOLDER=$1
REF=$2


for FILE in $FOLDER/*.sam.sorted.bam; do

  echo $FILE
  samtools mpileup -a -f $REF $FILE > $FILE.pileup.txt


done
