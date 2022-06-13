#!/bin/bash
#SBATCH --job-name=Rmasker
### each process only uses between 4 to 6 processes even though it says 32
#SBATCH -c 72
#SBATCH -o ./%x-%A-%a.out
#SBATCH -e ./%x-%A-%a.err
### 
#SBATCH --time=3-0:00
#SBATCH --mem=50G

cpus=$(($SLURM_CPUS_PER_TASK)/4)
LB=/vuw/tetools/Libraries
APP=tetools
s="singularity run -B $scra:/vuw -W /vuw/$APP --env LC_ALL=C --env LIBDIR=$LB --env BLASTDB=/vuw/$APP/databases $scra/$APP/tetools_grpiccoli.sif"

dbs=($(ls $scra/$APP/databases/*.translation | xargs -I{} basename {} .translation))

sOUT=/vuw/$APP/out
OUT=$scra/$APP/out
num=${1}
db=${dbs[$num]}
mkdir -p $OUT/${db}/Rmasker
cd $OUT/${db}/Rmasker
ln -s $OUT/${db}/consensi.fa.classified $OUT/${db}/Rmasker/consensi.fa.classified
$s RepeatMasker -libdir $LIBDIR -lib consensi.fa.classified -gff -pa $cpu -a -inv -xsmall -x -poly -source -html -ace -gff -e -u ${db}*.fa
