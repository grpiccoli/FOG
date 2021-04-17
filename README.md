# FOG
New Zealand Flat Oyster Genomics

[![singularity shub](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5003)

## [PIPELINE](https://kiwiepic.github.io/FOG/pages/flowchart.html)

![flowchart](https://kiwiepic.github.io/FOG/pages/flowchart.png)

## ARCHITECTURE

## QUICK START

```
curl -L -o ecoli.fastq https://sra-pub-src-1.s3.amazonaws.com/SRR10971019/m54316_180808_005743.fastq.1
```

### ASSEMBLY

```
nextflow run kiwiepic/FOG/<assembler>.nf -resume -profile local --ref ecoli.fastq --genomeSize 4.8m
```

assemblers: hifiasm | flye | pbipa | hicanu | peregrine

### Purge Dups

It is recommendable to run purge dups after assemblers hicanu and peregrine

```
nextflow run kiwiepic/FOG/purge_dups.nf -resume -profile local --ref ecoli.fastq
```