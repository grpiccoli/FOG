#!/usr/bin/env nextflow

raw_isoseq = params.isoseq
out_dir = file(params.outdir)

out_dir.mkdir()

isoseq_reads = Channel.fromPath("$raw_isoseq/**/*.fastq.gz", type: 'file').buffer(size:1)

process fastqc_isoseq {
	label 'fastqc'
	tag "$iso_sample"

	input:
	file iso_sample from isoseq_reads

	output:
	file "*_fastqc.{zip,html}" into iso_fastqc

	script:
	"""
	fastqc $iso_sample -t ${task.cpus} --noextract
	"""
}

process multiqc_isoseq {
    label 'multiqc'
    tag "$iso_fastqc"

    input:
	file ('iso_fastqc/*') from iso_fastqc.collect().ifEmpty([])

	output:
	file "multiqc_report.html" into iso_multiqc
	file "multiqc_data"

	script:
    """
    multiqc .
    """
}

workflow.onComplete {
	println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
 )
}
