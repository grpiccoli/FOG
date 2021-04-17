process racon {
	tag "racon.$x"

    input:
    file x from purge_dups

    output:
    file "*fasta" into racon

    when:
    params.all

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}