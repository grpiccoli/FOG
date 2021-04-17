#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_asm="$params.outdir/2_assemblers"

process peregrine {
	tag "peregrine.$x"
    container 'cschin/peregrine:1.6.3'
    publishDir out_asm

    input:
    file x from ref

    output:
    file "*fasta" into peregrine

    script:
    """
    echo "$x" > seq_dataset.lst
    . /opt/conda/etc/profile.d/conda.sh
    conda activate peregrine
    yes yes | pg_run.py asm seq_dataset.lst \
    ${task.cpus} ${task.cpus} \
    ${task.cpus} ${task.cpus} \
    ${task.cpus} ${task.cpus} \
    ${task.cpus} ${task.cpus} \
    ${task.cpus} \
    --with-consensus \
    --output peregrine
    """
}