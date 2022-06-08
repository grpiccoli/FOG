#!/bin/bash
#SBATCH --job-name=Rmodeler
### each process only uses between 4 to 6 processes even though it says 32
#SBATCH -c 32
#SBATCH -o ./%x-%A-%a.out
#SBATCH -e ./%x-%A-%a.err
### 
#SBATCH --time=3-0:00
#SBATCH --mem=50G

LB=/vuw/tetools/Libraries
APP=tetools
s="singularity run -B $scra:/vuw -W /vuw/tetools --env LC_ALL=C --env LIBDIR=$LB --env BLASTDB=/vuw/tetools/databases $scra/tetools/tetools_grpiccoli.sif"

dbs=($(ls $scra/tetools/databases/*.translation | xargs -I{} basename {} .translation))

sOUT=/vuw/tetools/out
OUT=$scra/tetools/out
num=${1}
#for db in ${dbs[@]}
#do
 db=${dbs[$num]}
 mkdir -p $OUT/${db}
 $s RepeatModeler -dir $sOUT/${db} -database /vuw/tetools/databases/$db -LTRStruct -pa 32 #2>&1 > $db.out #&
# pids[$i]=$!
# ((i++))
#done

#for pid in ${pids[@]}
#do
# wait $pid
#done
