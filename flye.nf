#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_i="$params.outdir/indexes"
out_asm="$params.outdir/2_assemblers"

process flye {
	tag "flye.$x"
    container "$params.bio/flye:2.8.2--py36h5202f60_0"
    publishDir "$out_asm"

    input:
    file x from ref

    output:
    file "flye/*fasta" into pre_flye
    file "flye/*{txt,log,gfa}"

    script:
    """
    if [[ -d "$out_asm/flye" ]];
    then
        ln -s $out_asm/flye/* .
    else
	    flye --pacbio-hifi $x -o flye -t $task.cpus -g $params.genomeSize -i 10
    fi
    """
}

process post_flye{
	tag "flye"
    publishDir "$out_asm/flye"

    input:
    file x from pre_flye.collect()

    output:
    file "flye.fasta" into flye

    script:
    """
    mv assembly.fasta flye.fasta
    """
}