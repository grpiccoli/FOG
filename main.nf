#!/usr/bin/env nextflow

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
out_dcn="$params.outdir/1_decontamination"
out_asm="$params.outdir/2_assemblers"
out_rpt="$params.outdir/3_repeat_analysis"
out_pol="$params.outdir/4_polishing"
out_phs="$params.outdir/5_phasing"
out_rph="$params.outdir/6_rephasing"
out_qul="$params.outdir/seq_quality/"
out_i="$params.outdir/indexes/"

//0.quality

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
    memory=`echo "$task.memory" | sed 's/[^0-9]//g' | awk '{print int(\$1/$task.cpus)"G"}'`
    export JAVA_TOOL_OPTIONS="-Xmx\$memory -Xms\$memory -Xss1G"
    fastqc -t $task.cpus -f \$format $x
    """
}

process multiqc {
    tag "multiqc.$x"
    container "$params.bio/multiqc:1.9--py_1"
    publishDir out_qul

	input:
    file "fastqc/*" from fastqc.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc

    script:
    """
    multiqc .
    """
}

//1.deconatamination

process pbmm2_icontaminants {
	tag "pbmm2_decontamination.$x"
    publishDir out_i

	input:
    file c from contaminants

	output:
    file "*ccs.mmi" into i_ccs
    file "*isoseq.mmi" into i_isoseq

	script:
	"""
    name="$c"
    name=\${name%.*}
    pbmm2 index --preset CCS $c \${name}_ccs.mmi
    pbmm2 index --preset ISOSEQ $c \${name}_isoseq.mmi
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
    tuple file("*.bam.pass"), file("*_cont.bam") into pbcont

	script:
	"""
    cpus=`echo "$task.cpus" | awk '{r=int(\$0/2);if(r==0){print 1}else{print r}}'`
    memory=`echo "$task.memory" | sed 's/[^0-9]//g' | awk -v cpus=\$cpus '{print int(\$0/cpus)}'`G
    preset=`echo ref_4W_A_hifi.bam | cut -d"_" -f4 \
    | awk -F '.' '{ print toupper(\$1) }'`
    if [[ preset == "isoseq" ]];
    then
        index=$iisoseq
    else
        index=$iccs
    fi
    name="$x"
    name=\${name%.*}
	pbmm2 align \$index $x \${name}_cont.bam --sort \
    -j \$cpus -J \$cpus -m \$memory
    ln -s $x ${x}.pass
	"""
}

process samtools_decon {
	tag "samtools_decon.$orig"
    publishDir out_dcn
    cache 'lenient'

	input:
	tuple file(orig), file(cont) from pbcont

	output:
    file "*_decon.bam" into decon
    file "*.con" into con
    file "*.cnt" into cnt

	script:
	"""
    name="$orig"
    name=\${name%.*}
    mv $orig \$name
    name=\${name%.*}
    samtools view $cont | cut -f 1 > \$name.con
    samtools view -h \$name.bam | grep -vf \$name.con \
    | samtools view -bS -o \${name}_decon.bam -
    samtools view -c \$name.bam > \$name.cnt
	"""
}

process gnuplot_decon {
	tag "r_graphdecon.$x"
    publishDir "$out_dcn/report"

	input:
	file x from con
    file n from cnt

	output:
    file "*pdf"
    file "*dat"
    file "*plg"

	script:
	"""
    name="$x"
    name=\${name%.*}
    data=\${name}.dat
    n=`cat $n`
    cat $x | cut -d"_" -f1 | sort | uniq -c | sort -k1nr,2n \
    | awk -F' ' -v cnt=0 -v num=\$n '{print cnt"\\t"\$2"\\t"\$1/num;cnt++;}' \
    > \$data
    tee -a \${name}.plg <<EOT
    set terminal pdfcairo enhanced color notransparent
    set output '\${name}_cont.pdf'
    set key noautotitle
    set xlabel "Taxonomic Order"
    set ylabel "Read Count"
    set title "Read Contamination by Taxonomic Order"
    set key box top left spacing 1.5
    plot '\$data' using 1:3:xtic(2) with boxes
    EOT
    gnuplot \${name}.plg
	"""
}

decon.branch {
    ref: it =~/ref/
    dnavar: it =~/dnavar/
    rnavar: it =~/dnavar/
}.set{decon}

decon.ref.into{
    ref_pbbam;
    ref_peregrine;
    ref_pbipa;
    ref_pb_assembly
}

//bam 2 fastq or fasta

process pbbam {
	tag "pbbam.$x"
    container "$params.bio/pbbam:1.6.0--h5b7e6e0_0"
    publishDir out_i

	input:
	file x from ref_pbbam

	output:
	tuple file("*.bam.pass"), file("*.pbi") into ref_pbi

	script:
	"""
    name="$x"
    pbindex $x
    ln -s $x ${x}.pass
	"""
}

process bam2fastx {
	tag "bam2fastq.$bam"
    container "$params.bio/bam2fastx:1.3.0--he1c1bb9_8"
    publishDir out_dcn
    cache 'lenient'

	input:
    tuple file(bam), file(index) from ref_pbi

	output:
	file "*.fastq.gz" into ref_hifiasm, ref_flye, ref_nextdenovo, ref_canu

	script:
	"""
    name="$bam"
    name=\${name%.*}
    mv $bam \$name
    name=\${name%.*}
    bam2fastq -o \$name \${name}.bam
	"""
}

//2.assemblers

process hifiasm {
	tag "hifiasm.$x"
    container "$params.bio/hifiasm:0.13--h8b12597_0"
    publishDir "$out_asm/hifiasm"

	input:
	file x from ref_hifiasm

	output:
	file "*.fasta" into hifiasm
    file "*asm*"

    //errorStrategy { task.exitStatus=0 ? 'ignore' : 'terminate' }

	script:
	"""
    gigs=`echo "$task.memory" | sed 's/[^0-9]//g'`
    if [[ \$gigs > 23 ]];
    then
        hifiasm -o hifiasm.asm -t $task.cpus $x 2> hifiasm.asm.log
        awk '/^S/{print ">"\$2;print \$3}' hifiasm.asm.p_ctg.gfa > hifiasm.fasta
    else
        touch hifiasm.asm.log
        touch hifiasm.fasta
    fi
	"""
}

process flye {
	tag "flye.$x"
    container "$params.bio/flye:2.8.2--py36h5202f60_0"
    publishDir "$out_asm/flye"

    input:
    file x from ref_flye

    output:
    file "*fasta" into flye
    file "*{txt,log}"

    script:
    """
	flye --pacbio-hifi $x -o flye -t $task.cpus -g $params.genomeSize -i 10
    mv flye/assembly.fasta flye.fasta
    mv flye/*txt .
    mv flye/*log .
    """
}

process nextdonovo {
	tag "nextdenovo.$x"
	
	input:
    file x from ref_nextdenovo

    output:
    file "ref.asm" into nextdenovo

    when:
    params.all

    script:
    """
	ls $x > input.fofn
	wget https://raw.githubusercontent.com/Nextomics/NextDenovo/master/doc/run.cfg
	nextDenovo run.cfg
    """
}

process pbipa {
	tag "pbipa.$x"
    container "$params.bio/pbipa:1.3.2--hee625c5_0"

    input:
    file x from ref_pbipa

    output:
    file "ref.asm" into pbipa

	when:
    params.all

    script:
    """
    ipa local -i $x --nthreads ${task.cpus} --njobs 1
    """
}

process canu {
	tag "canu.$x"
    container "${params.bio}/canu:2.1.1--he1b5a44_0"
    publishDir out_asm

	input:
	file x from ref_canu

	output:
	file "*.contigs.fasta" into canu
    file "*.report" canu_report
    file "*.fasta.gz" reads
    file "*.unassembled.fasta" unassembled
    file "*.layout.*" layouts

	when:
    params.all

	script:
	"""
    gigs=`echo "${task.memory}" | sed 's/[^0-9]//g'`
    if [[ \$gigs > 23 ]];
    then
        memory="\${gigs}G"
	    canu -assemble -p asm -d fog genomeSize=${params.genomeSize} \
        -pacbio-hifi $x useGrid=false \
        maxThreads=$task.cpus maxMemory=\$memory \
        merylThreads=$task.cpus merylMemory=\$memory merylConcurrency=1 \
        hapThreads=$task.cpus hapMemory=\$memory hapConcurrency=1 \
        batThreads=$task.cpus batMemory=\$memory batConcurrency=1 \
        gridOptionsJobName="$x"
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
    file "*fasta" into peregrine

	when:
    params.all

    script:
    """
    yes yes | python3 /data/korens/devel/Peregrine/bin/pg_run.py asm chm13.list 24 24 24 24 24 24 24 24 24 --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output ./
    """
}

process pb_assembly {
	tag "pb_assembly.$x"

    input:
    file x from ref_pb_assembly

    output:
    file "*fasta" into pb_assembly

    when:
    params.all || params.pb_assembly

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

//3.repeat analysis

canu.mix(peregrine,hifiasm,pbipa,flye,nextdenovo,pb_assembly)
.into{
    i_assembler_mummer;
    i_assembler_assembly_stats;
    i_assembler_quast;
    i_assembler_genomeqc
}

process mummer {
	tag "mummer.$x"

    input:
    file x from i_assembler_mummer

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

//4.polishing

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

//5.phasing

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

//6.re-phasing

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
    file "*fasta" into final_ref, quast_purge_haplotigs, mercury_purge_haplotigs, genomeqc_purge_haplotigs, assembly_stats_purge_haplotigs, i_prapi

    when:
    params.all || params.purge_haplotigs

    script:
    """
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}

//assembly quality assesment

process quast {
    tag "mummer.$x"

    input:
    file x from i_assembler_quast.mix(quast_dups,quast_racon,quast_nextpolish).collect()

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
    file x from i_assembler_genomeqc.mix(genomeqc_dups,genomeqc_racon,genomeqc_nextpolish).collect()

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
    tag "assembly-stats.$x"
    container "$params.fog/assembly-stats:17.02"
    publishDir "$out_asm/assembly-stats"

    input:
    file x from i_assembler_assembly_stats.mix(assembly_stats_purge_haplotigs).collect()

    output:
    file "*"

    script:
    """
    create-stats ${x} > assembly-stats.html
    """
}

decon.rnavar.into{i_isoseq; i_pbmm2}

//isoseq
process isoseq3 {
    input:
    file x from i_isoseq

    output:
    file "*fasta" into o_isoseq

    when:
    params.all

    script:
    """
    """
}

process prapi {
    input:
    file x from o_isoseq
    file g from i_prapi

    output:
    file "*fasta" into i_tama, i_maker

    when:
    params.all

    script:
    """
    """
}

//genome annotation
process tama {
    input:
    file x from i_tama

    output:
    file "*fasta" into tama

    when:
    params.all

    script:
    """
    """
}

process maker {
    input:
    file x from i_maker

    output:
    file "*fasta" into maker

    when:
    params.all

    script:
    """
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