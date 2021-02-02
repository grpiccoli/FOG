#!/usr/bin/env nextflow
out_dir = file(params.outdir)
out_dir.mkdir()

//isoseq
Channel.fromPath("$params.input/rnavar_*.bam", type: 'file')
.buffer(size:1).into{qrna;bamrna}
//hifi
Channel.fromPath("$params.input/dnavar_*.bam", type: 'file')
.buffer(size:1).into{qdna;bamdna}
//hifi reference
Channel.fromPath("$params.input/ref_*.bam", type: 'file')
.buffer(size:1).into{qref;bamref}
//hi-c reference
Channel.fromPath("$params.input/hic_*.fq.gz", type: 'file')
.buffer(size:2).into{qhic;hicref}
//Genomics_03\Tarakihi\TARdn2\Hifi\4_filter_conta_preassembly\refseq
contaminants=file("$params.input/contaminants.fasta")
repbase=file("$params.input/RepeatMasker*.fasta")
out_dcn="$out_dir/1_decontamination"
out_asm="$out_dir/2_assemblers"
out_rpt="$out_dir/3_repeat_analysis"
out_pol="$out_dir/4_polishing"
out_phs="$out_dir/5_phasing"
out_rph="$out_dir/6_rephasing"

process fastqc {
    tag "fastqc.$x"
    //error
    container "$params.bio/fastqc:0.11.9--0"
    //container "$params.fog/fastqc:0.11.9"

    input:
    file x from qdna.mix(qref,qrna,qhic)

    output:
    file "*_fastqc.zip" into fastqc

    script:
    """
    if [[ "$x" == *"bam"* ]];
    then
        format=bam
    else
        format=fastq
    fi
    fastqc -t $task.cpus -f \$format $x
    """
}

process multiqc {
    tag "multiqc.$x"
    container "$params.bio/multiqc:1.9--py_1"
    publishDir params.input

	input:
    file "fastqc/*" from fastqc.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc

    script:
    """
    multiqc .
    """
}

process pbmm2_icontaminants {
	tag "pbmm2_decontamination.$x"
    publishDir params.input

	input:
    file c from contaminants

	output:
    file "*ccs.mmi" into i_ccs
    file "*isoseq.mmi" into i_isoseq

	script:
	"""
    pbmm2 index --preset CCS $c ${c}_ccs.mmi
    pbmm2 index --preset ISOSEQ $c ${c}_isoseq.mmi
	"""
}

process pbmm2_mcontaminats {
	tag "pbmm2_mcontaminants.$x"
    publishDir out_dcn

	input:
	file x from bamref.mix(bamdna,bamrna)
    file iccs from i_ccs.collect()
    file iisoseq from i_isoseq.collect()

	output:
    file "*cont.bam" into pbcont

	script:
	"""
    memory=`echo "${task.memory}" | sed 's/[^0-9]//g'`G
    preset=`echo ref_4W_A_hifi.bam | cut -d"_" -f4 | awk -F '.' '{ print toupper(\$1) }'`
    if [[ preset == "isoseq" ]];
    then
        index=$iisoseq
    else
        index=$iccs
    fi
	pbmm2 align \$index $x ${x}.cont.bam --sort -j ${task.cpus} -J ${task.cpus} -m \$memory
	"""
}

process samtools_decon {
	tag "samtools_decon.$x"
    publishDir out_dcn

	input:
	file x from pbcont

	output:
    file "*decon.bam" into decon
    file "*.con" into con
    file "*.cnt" into cnt

	script:
	"""
    bam=`echo "$x" | sed 's/\\.cont\\.bam\$//'`
    ln -s "$workflow.launchDir/$params.input/\$bam" .
    samtools view $x | cut -f 1 > ${x}.con
    samtools view -h \$bam | grep -vf ${x}.con | samtools view -bS -o \${bam%.*}.decon.bam -
    samtools view -c \$bam > ${x}.cnt
	"""
}

process gnuplot_decon {
	tag "r_graphdecon.$x"
    publishDir out_dcn

	input:
	file x from con
    file n from cnt

	output:
    file "*pdf"
    file "*dat"
    file "*plg"

	script:
	"""
    data=dataframe.dat
    n=`cat $n`
    cat $x | cut -d"_" -f1 | sort | uniq -c | sort -k1nr,2n | awk -F' ' -v cnt=0 -v num=\$n '{print cnt"\\t"\$2"\\t"\$1/num;cnt++;}' > \$data
    tee -a tmp.plg <<EOT
    set terminal pdfcairo enhanced color notransparent
    set output '${x}.pdf'
    set key noautotitle
    set xlabel "Taxonomic Order"
    set ylabel "Read Count"
    set title "Read Contamination by Taxonomic Order"
    set key box top left spacing 1.5
    plot '\$data' using 1:3:xtic(2) with boxes
    EOT
    gnuplot tmp.plg
	"""
}

decon.branch {
    ref: it =~/ref/
    dnavar: it =~/dnavar/
    rnavar: it =~/dnavar/
}.set{decon}

decon.ref.into{ref_pbbam; ref_peregrine; ref_hifiasm; ref_flye; ref_pbipa; ref_nextdenovo; ref_pb_assembly}

process pbbam {
	tag "pbbam.$x"
    container "$params.bio/pbbam:1.6.0--h5b7e6e0_0"
    publishDir params.input

	input:
	file x from ref_pbbam

	output:
	file "*.pbi" into pbi

	script:
	"""
    pbindex $x
	"""
}

process bam2fastx {
	tag "bam2fastq.$x"
    container "$params.bio/bam2fastx:1.3.0--he1c1bb9_8"
    publishDir params.input

	input:
    file x from pbi

	output:
	file "*.fastq.gz" into ref_canu

	script:
	"""
    bam="$x"
    bam="\${bam%.*}"
    ln -s "$out_dcn/\$bam" .
    bam2fastq -o \${bam%.*} \$bam
	"""
}

process canu {
	tag "canu.$x"
    container "${params.bio}/canu:2.1.1--he1b5a44_0"
    publishDir out_asm

	input:
	file x from ref_canu

	output:
	file "*.contigs.fasta" into canu, quast_canu, genomeqc_canu
    file "*.report" canu_report
    file "*.fasta.gz" reads
    file "*.unassembled.fasta" unassembled
    file "*.layout.*" layouts

	script:
	"""
    gigs=`echo "${task.memory}" | sed 's/[^0-9]//g'`
    if [[ \$gigs > 23 ]];
    then
        memory="\$gigs\\G"
	    canu -assemble -p asm -d fog genomeSize=1g -pacbio-hifi $x useGrid=false maxThreads=${task.cpus} maxMemory=\$memory
    else
        touch ${x}.contigs.layout.tigInfo
        touch ${x}.contigs.layout.readToTig
        touch ${x}.contigs.layout
        touch ${x}.unassembled.fasta
        touch ${x}.contigs.fasta
        touch ${x}.trimmedReads.fasta.gz
        touch ${x}.correctedReads.fasta.gz
    fi
	"""
}

process peregrine {
	tag "peregrine.$x"
    container 'docker://cschin/peregrine:1.6.3'
    publishDir out_asm

    input:
    file x from ref_peregrine

    output:
    file "*fasta" into peregrine, quast_peregrine, genomeqc_peregrine

	when:
    params.all

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
	file "ref.asm" into hifiasm, quast_hifiasm, genomeqc_hifiasm

	when:
    params.all

	script:
	"""
	hifiasm -o ref.asm -t${task.cpus} $x
    awk '/^S/{print ">"\$2;print \$3}' ref.asm.p_ctg.gfa > ref.asm.fasta
	"""
}

process pbipa {
	tag "pbipa.$x"

    input:
    file x from ref_pbipa

    output:
    file "ref.asm" into pbipa, quast_pbipa, genomeqc_pbipa

    when:
    params.all || params.pbipa

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
    file "ref.asm" into flye, quast_flye, genomeqc_flye

    when:
    params.all || params.flye

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
    file "ref.asm" into nextdenovo, quast_nextdenovo, genomeqc_nextdenovo

    when:
    params.all || params.nextdenovo

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
    file "*fasta" into pb_assembly, quast_pb, genomeqc_pb

    when:
    params.all || params.pb_assembly

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
    params.all || params.mummer

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
    params.all || params.transposonpsi

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
    params.all || params.tetools

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
    params.all || params.mummer

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
    file "*fasta" into purge_dups, quast_dups, genomeqc_dups

    when:
    params.all || params.purge_dups

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
    file "*fasta" into racon, quast_racon, genomeqc_racon

    when:
    params.all || params.racon

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
    file "*fasta" into i_allhic, i_marginphase, i_falconphase, i_hirise, i_dipasm, quast_nextpolish, genomeqc_nextpolish

    when:
    params.all || params.nextpolish

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
    file "*fasta" into allhic, mercury_allhic

    when:
    params.all || params.allhic

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
    file "*fasta" into marginphase, mercury_marginphase

    when:
    params.all || params.marginphase

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
    file "*fasta" into falconphase, mercury_falconphase

    when:
    params.all || params.falconphase

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
    file "*fasta" into hirise, mercury_hirise

    when:
    params.all || params.hirise

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process dipasm {
	tag "nextpolish.$x"

    input:
    file x from i_dipasm

    output:
    file "*fasta" into dipasm, mercury_dipasm

    when:
    params.all || params.hirise

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

allhic.mix(marginphase, falconphase, hirise, dipasm)
.into{ i_haplotypo; i_purge_haplotigs }

process haplotypo {
	tag "haplotypo.$x"

    input:
    file x from i_haplotypo

    output:
    file "*fasta" into haplotypo

    when:
    params.all || params.hirise

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process purge_haplotigs {
    tag "purge_haplotigs.$x"

    input:
    file x from i_purge_haplotigs

    output:
    file "*fasta" into final_ref, quast_purge_haplotigs, mercury_purge_haplotigs, genomeqc_purge_haplotigs, assembly_stats_purge_haplotigs

    when:
    params.all || params.purge_haplotigs

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

process quast {
    tag "mummer.$x"

    input:
    file x from quast_canu.mix(quast_peregrine,quast_hifiasm,quast_pbipa,quast_flye,quast_nextdenovo,quast_pb,quast_dups,quast_racon,quast_nextpolish).collect()

    output:
    file "*fasta" into quast
    file "other" into squat

    when:
    params.all || params.mummer

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process icarus {
    tag "icarus.$x"

    input:
    file x from quast

    output:
    file "*fasta" into icarus

    when:
    params.all || params.icarus

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process genomeqc {
    tag "genomeqc.$x"

    input:
    file x from genomeqc_canu.mix(genomeqc_peregrine,genomeqc_hifiasm,genomeqc_pbipa,genomeqc_flye,genomeqc_nextdenovo,genomeqc_pb,genomeqc_dups,genomeqc_racon,genomeqc_nextpolish).collect()

    output:
    file "*fasta" into genomeqc

    when:
    params.all || params.mummer

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process merqury {
    tag "mummer.$x"

    input:
    file x from mercury_allhic.mix(mercury_marginphase,mercury_falconphase,mercury_hirise,mercury_dipasm).collect()

    output:
    file "*fasta" into merqury

    when:
    params.all || params.mummer

    script:
    """
    Quast.py --large --skip-unaligned-mis-contigs    
    """
}

process assembly_stats {
    tag "mummer.$x"

    input:
    file x from assembly_stats_purge_haplotigs

    output:
    file "*fasta" into assembly_stats

    when:
    params.all || params.mummer

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