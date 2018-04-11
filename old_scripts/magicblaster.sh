#!/bin/bash
#aabspliced
## makeblastdb -in aa-bBACclone.fasta -dbtype nucl -parse_seqids -out aabBACclone
# magicblast -db aabBACclone -sra SRR1260867 -no_unaligned -num_threads 4 -score 40 -out SRR1260867.aab.test.foo
## If you are using RNA then add -splice F
## awk '{if ($NF == "NM:i:0") print $0}' ERR1024254.transcriptome.aabspliced | sort -nk 4

magicblast -db aabspliced -sra $1 -no_unaligned -num_threads 4 -reftype transcriptome -out results/$1.aml_aabspliced
