#!/usr/bin/env nextflow
Channel.fromPath("$params.asm", type: 'file')
.buffer(size:1)
.set{
    asm;
}

Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_asm="$params.outdir/2_assemblers"
out_i="$params.outdir/indexes"

process purge_dups_yaml {
	tag "yamls"
  publishDir out_i

  output:
  file "purge_dups.yaml" into purge_dups_yaml

  script:
  """
  tee -a purge_dups.yaml <<EOT
  name: purge_dups
  channels:
    - defaults
    - conda-forge
    - bioconda
  dependencies:
    - purge_dups=1.2.5
  EOT
  """
}

process purge_dups {
	tag "purge_dups.$x"
  publishDir out_asm
  conda "$out_i/$y"

  input:
  file x from ref
  file a from asm
  file y from purge_dups_yaml

  output:
  file "*fasta" into purge_dups, purge_dups_stats

  script:
  """
  wget https://raw.githubusercontent.com/dfguan/purge_dups/master/scripts/pd_config.py
  chmod +x pd_config.py
  echo "$x" > pb.fofn
  ./pd_config.py \
  -l purge_dups \
  -n config.purge_dups.json \
  $a \
  pb.fofn
  run_purge_dups.py config.purge_dups.json purge_dups purge_dups
  """
}