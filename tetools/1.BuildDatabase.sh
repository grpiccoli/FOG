#!/bin/bash

LB=/vuw/tetools/Libraries

APP=tetools

s="singularity run -B $scra:/vuw -W /vuw/tetools --env LC_ALL=C --env LIBDIR=$LB --env BLASTDB=/vuw/tetools/databases $scra/$APP/tetools_latest.sif"

asm_dir=input/asm

asms=($(ls $scra/$asm_dir/*.fa))

for asm in "${asms[@]}"
do
 base=$(basename $asm .fa)
 name=${base%%.*}
 ls $scra/$asm_dir/${name}* > lst
 cmd="$s BuildDatabase -name $name -batch lst
 echo $cmd
 eval "$cmd"

done
