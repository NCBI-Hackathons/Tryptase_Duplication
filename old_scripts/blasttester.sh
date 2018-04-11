#!/usr/bin/bash

 magicblast -db $1 -sra ERR1024254 -no_unaligned -num_threads 4 -splice F -reftype transcriptome -out ERR1024254.transcriptome.$1 &
