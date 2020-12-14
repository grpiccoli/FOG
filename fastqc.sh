#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=12
#SBATCH -o /nfs/scratch/rodriggu/fastqc.out
#SBATCH -e /nfs/scratch/rodriggu/fastqc.err
#SBATCH --mail-type=BEGIN,END.FAIL
#SBATCH --mail-user=guillermo.rodriguezpiccoli@vuw.ac.nz
#SBATCH --time=6-12:00
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --partition=parallel
module load singularity
out_d=/nfs/scratch/rodriggu/fastqc
mkdir -p $out_d
ls -R isoseq/ | grep fastq.gz | xargs -P4 -i time singularity run /nfs/scratch/rodriggu/FOG_fastqc.sif {} -t 12 --noextract -o $out_d
