#!/usr/bin/env nextflow
Channel.fromPath("$params.refbam", type: 'file')
.buffer(size:1)
.set{
    refbam;
}

out_i="$params.outdir/indexes"
out_ref="$params.outdir/fastas"

process pbbam {
	tag "pbbam.$x"
    container "$params.bio/pbbam:1.6.0--h5b7e6e0_0"
    publishDir out_i

	input:
	file x from refbam

	output:
	tuple file("*.bam.pass"), file("*.pbi") into ref_pbi

	script:
	"""
    pbindex $x
    ln -s $x ${x}.pass
	"""
}

process bam2fastq {
	tag "bam2fastq.$bam"
    container "$params.bio/bam2fastx:1.3.0--he1c1bb9_8"
    publishDir out_ref
    cache 'lenient'

	input:
    tuple file(bam), file(index) from ref_pbi

	output:
	file "*.fastq.gz" into ref_hifiasm, ref_flye, ref_nextdenovo, ref_canu

	script:
	"""
    name=$bam
    name=\${name%.*}
    mv $bam \$name
    name=\${name%.*}
    bam2fastq -o \$name \${name}.bam
	"""
}