!/usr/bin/env nextflow

raw_reads = params.rawReads
out_dir = file(params.outDir)

out_dir.mkdir()

process fastqc {
	tag { "${params.projectName}.rFQC.${sample}" }
	cpus { 2 }
	publishDir "${out_dir}/qc/raw/${sample}", mode: 'copy', overwrite: false

	input:
		set sample, file(in_fastq) from read_pair
}
