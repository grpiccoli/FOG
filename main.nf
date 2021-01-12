#!/usr/bin/env nextflow
out_dir = file(params.outdir)

out_dir.mkdir()

bamiso = Channel.fromPath("$params.isoseq/**/*.bam", type: 'file').buffer(size:1)
bamref = Channel.fromPath("$params.ref/*.bam", type: 'file').buffer(size:1)
bamvar = Channel.fromPath("$params.hifi/**/*.bam", type: 'file').buffer(size:1)
hicref = Channel.fromPath("$params.hic/**/*.fq.gz", type: 'file').buffer(size:2)
bams = bamiso.mix(bamref,bamvar)

process bam2fastx {
	tag "bam2fastq.$x"

	input:
	file x from bams

	output:
	file "*.fastq.gz" into fastq

	when:
	params.run == 'all' || params.run == 'bam2fastx'

	script:
	"""
	bam2fastq -o $x $x
	"""
}

process fastqc {
    tag "fastqc.$x"

    input:
    file x from fastq.mix(hicref)

    output:
    file "*_fastqc.{zip,html}" into fastqc

	when:
	params.run == 'all' || params.run == 'fastqc'

    script:
    """
    fastqc $x -t ${task.cpus} --noextract
    """
}

process multiqc_ref {
    tag "multiqc_iso.$x"

    input:
	ref = fastqc.filter(~/.*${params.refname}.*/).collect().ifEmpty([])
	hic = fastqc.filter(~/.*${params.refname}.*/).collect().ifEmpty([])
	var = fastqc.filter(~/.*${params.refname}.*/).collect().ifEmpty([])
	iso = fastqc.filter(~/.*${params.refname}.*/).collect().ifEmpty([])
    file ('*') from ref.mix(hic,var,iso)

    output:
    file "multiqc_report.html" into multiqc

	when:
    params.run == 'all' || params.run == 'multiqc'

    script:
    """
    multiqc .
    """
}

process canu {
	tag "canu.$x"

	input:
	file x from fastqref

	output:
	file "*_fastqc.{zip,html}" into iso_fastqc

	when:
    params.run == 'all' || params.run == 'canu'

	script:
	"""
	canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
	"""
}

process peregrine {
	tag "peregrine.$x"

    input:
    file x from fastqref

    output:
    file "*_fastqc.{zip,html}" into iso_fastqc

	when:
    params.run == 'all' || params.run == 'peregrine'

    script:
    """
    yes yes | python3 /data/korens/devel/Peregrine/bin/pg_run.py asm chm13.list 24 24 24 24 24 24 24 24 24 --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output ./
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
