#!/usr/bin/env nextflow
out_dir = file(params.outdir)

out_dir.mkdir()

//isoseq
Channel.fromPath("$params.isoseq/**/*.bam", type: 'file')
.buffer(size:1).set{bamiso}
//hifi reference
Channel.fromPath("$params.ref/*.bam", type: 'file')
.buffer(size:1).set{bamref}
//hifi variants
Channel.fromPath("$params.hifi/**/*.bam", type: 'file')
.buffer(size:1).set{bamccs}
//hi-c reference
Channel.fromPath("$params.hic/**/*.fq.gz", type: 'file')
.buffer(size:2).set{hicref}

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

canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)
.set{i_mummer}

process mummer {
	tag "mummer.$x"

    input:
    file x from i_mummer

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

process transposonpsi {
    tag "transposonpsi.$x"

    input:
    file x from mummer

    output:
    file "*fasta" into transposonpsi

    when:
    params.run == 'all' || params.run == 'transposonpsi'

    script:
    """
    transposonPSI.pl $x nuc
    """
}

process tetools {
    tag "tetools.$x"

    input:
    file x from transposonpsi

    output:
    file "*fasta" into tetools

    when:
    params.run == 'all' || params.run == 'tetools'

    script:
    """
    BuildDatabase -name fog -engine ncbi $x
    RepeatModeler -database fog -engine ncbi -pa ${task.cpus} -LTRStruct
    RepeatMasker -lib fog-families.fa $x -pa ${task.cpus}
    """
}

process minimap2 {
    tag "mummer.$x"

    input:
    file x from tetools

    output:
    file "*fasta" into minimap2

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
    file x from minimap2

    output:
    file "*fasta" into purge_dups

    when:
    params.run == 'all' || params.run == 'purge_dups'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process racon {
	tag "racon.$x"

    input:
    file x from purge_dups

    output:
    file "*fasta" into racon

    when:
    params.run == 'all' || params.run == 'racon'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process nextpolish {
	tag "nextpolish.$x"

    input:
    file x from racon

    output:
    file "*fasta" into i_allhic, i_marginphase, i_falconphase, i_hirise

    when:
    params.run == 'all' || params.run == 'nextpolish'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process allhic {
	tag "nextpolish.$x"

    input:
    file x from i_allhic

    output:
    file "*fasta" into allhic

    when:
    params.run == 'all' || params.run == 'allhic'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process marginphase {
	tag "nextpolish.$x"

    input:
    file x from i_marginphase

    output:
    file "*fasta" into marginphase

    when:
    params.run == 'all' || params.run == 'marginphase'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process falconphase {
	tag "nextpolish.$x"

    input:
    file x from i_falconphase

    output:
    file "*fasta" into falconphase

    when:
    params.run == 'all' || params.run == 'falconphase'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process hirise {
	tag "nextpolish.$x"

    input:
    file x from i_hirise

    output:
    file "*fasta" into hirise

    when:
    params.run == 'all' || params.run == 'hirise'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

allhic.mix(marginphase, falconphase, hirise).into{ i_deepv; i_haplotypo }

process deepv {
	tag "deepv.$x"

    input:
    file x from i_deepv

    output:
    file "*fasta" into deepv

    when:
    params.run == 'all' || params.run == 'hirise'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process haplotypo {
	tag "haplotypo.$x"

    input:
    file x from i_haplotypo

    output:
    file "*fasta" into haplotypo

    when:
    params.run == 'all' || params.run == 'hirise'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

deepv.mix(haplotypo).set{ i_whatshap }

process whatshap {
	tag "whatshap.$x"

    input:
    file x from i_whatshap

    output:
    file "*fasta" into whatshap

    when:
    params.run == 'all' || params.run == 'whatshap'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process rephase {
    tag "rephase.$x"

    input:
    file x from whatshap

    output:
    file "*fasta" into rephase

    when:
    params.run == 'all' || params.run == 'rephase'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process purge_haplotigs {
    tag "purge_haplotigs.$x"

    input:
    file x from rephase

    output:
    file "*fasta" into i_quast, i_mercury, i_genomeqc, i_assembly_stats, final_ref

    when:
    params.run == 'all' || params.run == 'purge_haplotigs'

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process quast {
    tag "mummer.$x"

    input:
    file x from i_quast

    output:
    file "*fasta" into quast

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process merqury {
    tag "mummer.$x"

    input:
    file x from i_mercury

    output:
    file "*fasta" into merqury

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process genomeqc {
    tag "mummer.$x"

    input:
    file x from i_genomeqc

    output:
    file "*fasta" into genomeqc

    when:
    params.run == 'all' || params.run == 'mummer'

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process assembly_stats {
    tag "mummer.$x"

    input:
    file x from i_assembly_stats

    output:
    file "*fasta" into assembly_stats

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