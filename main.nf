#!/usr/bin/env nextflow
out_dir = file(params.outdir)

out_dir.mkdir()

bamiso = Channel.fromPath("$params.isoseq/**/*.bam", type: 'file').buffer(size:1)
bamref = Channel.fromPath("$params.ref/*.bam", type: 'file').buffer(size:1)
bamccs = Channel.fromPath("$params.hifi/**/*.bam", type: 'file').buffer(size:1)
hicref = Channel.fromPath("$params.hic/**/*.fq.gz", type: 'file').buffer(size:2)

process bam2fastx {
	tag "bam2fastq.$x"

	input:
	file x from bamiso.mix(bamref,bamccs)

	output:
	file "*.fastq.gz" into fastq
	file "*$params.refname*" into ref_canu, ref_peregrine, ref_hifiasm, ref_flye, ref_pbipa, ref_nextdenovo, ref_pb_assembly

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

process multiqc {
    tag "multiqc.$x"

	input:
    tuple x, file('*') from fastqc.map { 
	if (it =~/.*ref.*/){  
		return ['ref', it]  
	}else if(it =~/.*hic.*/){ 
		return ['hic', it]  
	}else if(it =~/.*ccs.*/){ 
		return ['css', it]  
	}else if(it =~/.*iso.*/){ 
		return ['iso', it]  
	}  
	} 
	.groupTuple()

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
	file x from ref_canu

	output:
	file "*fasta" into canu

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
    file x from ref_peregrine

    output:
    file "*fasta" into peregrine

	when:
    params.run == 'all' || params.run == 'peregrine'

    script:
    """
    yes yes | python3 /data/korens/devel/Peregrine/bin/pg_run.py asm chm13.list 24 24 24 24 24 24 24 24 24 --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output ./
    """
}

process hifiasm {
	tag "hifiasm.$x"

	input:
	file x from ref_hifiasm

	output:
	file "ref.asm" into hifiasm

	when:
	params.run == 'all' || params.run == 'hifiasm'

	script:
	"""
	hifiasm -o ref.asm -t${task.cpus} $x
	"""
}

process pbipa {
	tag "pbipa.$x"

    input:
    file x from ref_pbipa

    output:
    file "ref.asm" into pbipa

    when:
    params.run == 'all' || params.run == 'pbipa'

    script:
    """
    ipa local --nthreads ${task.cpus} --njobs 1 -i $x
    """
}

process flye {
	tag "flye.$x"

    input:
    file x from ref_flye

    output:
    file "ref.asm" into flye

    when:
    params.run == 'all' || params.run == 'flye'

    script:
    """
	flye --pacbio-hifi $x -o large -t ${task.cpus} -g 1g -i 20
	flye --pacbio-hifi $x -o short -t ${task.cpus} -g 0.6g -i 20
    """
}

process nextdonovo {
	tag "nextdenovo.$x"
	
	input:
    file x from ref_nextdenovo

    output:
    file "ref.asm" into nextdenovo

    when:
    params.run == 'all' || params.run == 'nextdenovo'

    script:
    """
	ls $x > input.fofn
	cp doc/run.cfg ./
	nextDenovo run.cfg
    """
}

process pb_assembly {
	tag "pb_assembly.$x"

    input:
    file x from ref_pb_assembly

    output:
    file "*fasta" into pb_assembly

    when:
    params.run == 'all' || params.run == 'pb_assembly'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process mummer {
	tag "mummer.$x"

    input:
    file x from canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)

    output:
    file "*fasta" into mummer

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    nucmer --maxmatch --nosimplify delta-filter -i 98 -l 10000
    nucmer --maxmatch --noextend --nosimplify -l 500 -c 1000 delta-filter -i 99.9 -l 10000
    """
}

process quast {
    tag "mummer.$x"

    input:
    file x from canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)

    output:
    file "*fasta" into mummer

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process tetools {
    tag "mummer.$x"

    input:
    file x from canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)

    output:
    file "*fasta" into mummer

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process minimap2 {
    tag "mummer.$x"

    input:
    file x from canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)

    output:
    file "*fasta" into mummer

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process purge_dups {
	tag "purge_dups.$x"

    input:
    file x from canu

    output:
    file "*fasta" into purge_dups

    when:
    params.run == 'all' || params.run == 'purge_dups'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}



process merqury {
    tag "mummer.$x"

    input:
    file x from canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)

    output:
    file "*fasta" into mummer

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
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
