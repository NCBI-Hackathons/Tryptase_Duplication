#!/bin/bash

ACC_FILE=$1

for i in $( cat $ACC_FILE ); do

  magicblast -db aabspliced -sra $i -no_unaligned -num_threads 4 -reftype transcriptome -out results/$i.aml_aabspliced

done
