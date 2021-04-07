#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_asm="$params.outdir/2_assemblers"

process canu {
	tag "canu.$x"
    container "$params.bio/canu:2.1.1--he1b5a44_0"
    publishDir out_asm

	input:
	file x from ref

	output:
	file "*.contigs.fasta" into canu
    file "*.report" canu_report
    file "*.fasta.gz" reads
    file "*.unassembled.fasta" unassembled
    file "*.layout.*" layouts

	script:
	"""
    memory=`echo "${task.memory}" | sed 's/[^0-9]//g'`
    memory="\${memory}G"
    canu -assemble -p asm -d fog genomeSize=${params.genomeSize} \
    -pacbio-hifi $x useGrid=false \
    maxThreads=$task.cpus maxMemory=\$memory \
    merylThreads=$task.cpus merylMemory=\$memory merylConcurrency=1 \
    hapThreads=$task.cpus hapMemory=\$memory hapConcurrency=1 \
    batThreads=$task.cpus batMemory=\$memory batConcurrency=1 \
    gridOptionsJobName="$x"
	"""
}