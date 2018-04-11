#!/bin/bash

for i in * ; do

   LENGTH=${#i}

   END=${i:LENGTH-3:LENGTH}

   if [ $END != "sam" ] && [ $END != ".sh" ]

   then

     echo $i
     mv $i $i.sam

   fi

done
