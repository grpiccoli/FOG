#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_i="$params.outdir/indexes"
out_asm="$params.outdir/2_assemblers"

process nextdenovo {
	tag "nextdenovo.$x"
    container "$params.fog:nextdenovo-2.4.0"
	
	input:
    file x from ref

    output:
    file "ref.asm" into nextdenovo, nextdenovo_stats

    script:
    """
	ls $x > input.fofn
	wget https://raw.githubusercontent.com/Nextomics/NextDenovo/master/doc/run.cfg
	nextDenovo run.cfg
    """
}

process nextpolish {
	tag "nextpolish.$x"

    input:
    file x from nextdenovo

    output:
    file "*fasta" into nextpolish

    script:
    """
    canu -assemble -p asm -d asm genomeSize=$params.genomeSize -pacbio-hifi $x
    """
}