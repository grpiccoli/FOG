#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_asm="$params.outdir/seq_quality"

process fastk {
	tag "pb_assembly.$x"
    publishDir out_asm
    container "grpiccoli/fastk:latest"

    input:
    file r from ref

    output:
    file "*fasta" into longqc

    script:
    """
    mem=`echo "$task.memory" | sed 's/[^0-9]*//g'`
    mem=`expr \$mem / $task.cpus`
    longQC.py sampleqc \
    -p $task.cpus \
    -m \$mem \
    -x pb-rs2 \
    -o longqc $r \
    --index $params.index
    """
}