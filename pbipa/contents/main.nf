#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_i="$params.outdir/indexes"
out_asm="$params.outdir/2_assemblers"

process yamls {
	tag "yamls"
  publishDir out_i

  output:
  file "pbipa.yaml" into pbipa_yaml

  script:
  """
  tee -a pbipa.yaml <<EOT
  name: pbipa
  channels:
    - defaults
    - conda-forge
    - bioconda
  dependencies:
    - pbipa=1.3.2
  EOT
  """
}

process pbipa {
	tag "pbipa.$x"
  publishDir out_asm
  conda "$out_i/$y"
  //container "$params.bio/pbipa:1.3.2--hee625c5_0"

  input:
  file x from ref
  file y from pbipa_yaml

  output:
  file "pbipa/assembly-results/*.fasta" into pre_pbipa

  script:
  """
  #--no-phase
  ipa local --nthreads ${task.cpus} --njobs 1 --run-dir pbipa -i $x
  """
}

process post_pbipa {
  tag "post_pbipa"
  publishDir "$out_asm/pbipa"

  input:
  file x from pre_pbipa.collect()

  output:
  file "*.fasta" into pbipa

  script:
  """
  mv final.a_ctg.fasta pbipa_alternate.fasta
  mv final.p_ctg.fasta pbipa_primary.fasta
  """
}