#!/usr/bin/env nextflow
Channel.fromPath("$params.bam", type: 'file')
.buffer(size:1)
.set{
    refbam;
}

Channel.fromPath("$params.fasta", type: 'file')
.buffer(size:1)
.set{
    reffa;
}

Channel.fromPath("$params.hic1", type: 'file')
.buffer(size:1)
.set{
    hic1;
}

Channel.fromPath("$params.hic2", type: 'file')
.buffer(size:1)
.set{
    hic2;
}

out_i="$params.outdir/indexes"
out_asm="$params.outdir/2_assemblers"

process yamls {
  tag "yamls"
  publishDir out_i

  output:
  file "pb-assembly.yaml" into pbipa_yaml

  script:
  """
  tee -a pb-assembly.yaml <<EOT
  name: pb-assembly
  channels:
    - defaults
    - conda-forge
    - bioconda
  dependencies:
    - pb-assembly=0.0.8-1
  EOT
  """
}

process pb_assembly {
	tag "pb_assembly.$x"
    publishDir out_asm
    conda "$out_i/$y"

    input:
    file f from reffa
    file b from refbam
    file c from hic1
    file d from hic2
    file y from pbipa_yaml

    output:
    file "*fasta" into pb_assembly

    script:
    """
    echo "$b" > subreads.bam.fofn
    echo "$f" > subreads.fasta.fofn
    tee -a fc_run.cfg <<EOT

    EOT
    tee -a fc_unzip.cfg <<EOT
    
    EOT
    tee -a fc_phase.cfg <<EOT
    
    EOT
    canu -assemble -p asm -d asm genomeSize=0.6g -pacbio-hifi $x
    """
}