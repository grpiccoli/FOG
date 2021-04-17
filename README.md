# FOG
New Zealand Flat Oyster Genomics

[![singularity shub](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5003)

## [PIPELINE](https://kiwiepic.github.io/FOG/pages/flowchart.html)

![flowchart](https://kiwiepic.github.io/FOG/pages/flowchart.png)

## ARCHITECTURE

## QUICK START

curl -L -o ecoli.fastq https://sra-pub-src-1.s3.amazonaws.com/SRR10971019/m54316_180808_005743.fastq.1

### Hifiasm

nextflow run kiwiepic/FOG/hifiasm.nf -resume -profile local --ref rapoi/ecoli.fastq --genomeSize 1g

### Flye

nextflow run kiwiepic/FOG/flye.nf -resume -profile local --ref rapoi/ecoli.fastq --genomeSize 1g

### Hi-Canu

nextflow run kiwiepic/FOG/hicanu.nf -resume -profile local --ref rapoi/ecoli.fastq --genomeSize 1g

nextflow run kiwiepic/FOG/purge_dups.nf -resume -profile local --ref output/1-assembly/hicanu/*.fa

### Peregrine

