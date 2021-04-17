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

    input:
    file x from canu
    file y from purge_dups_yaml

    output:
    file "*fasta" into purge_dups, purge_dups_stats

    script:
    """
    ./scripts/pd_config.py \
    -l iHelSar1.pri -s 10x.fofn \
    -n config.iHelSar1.PB.asm1.json \
    ~/vgp/release/insects/iHelSar1/iHlSar1.PB.asm1/iHelSar1.PB.asm1.fa.gz \
    pb.fofn

    python scripts/run_purge_dups.py config.iHelSar1.json src iHelSar1
    """
}