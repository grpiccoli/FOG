#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_asm="$params.outdir/2_assemblers"
out_i="$params.outdir/indexes"

process canu_yaml {
	tag "yaml"
  publishDir out_i

  output:
  file "canu.yaml" into canu_yaml

  script:
  """
  tee -a canu.yaml <<EOT
  name: canu
  channels:
    - defaults
    - conda-forge
    - bioconda
  dependencies:
    - canu=2.1.1
  EOT
  """
}

process canu {
	tag "canu.$x"
  //container "$params.bio/canu:2.1.1--he1b5a44_0"
  publishDir out_asm
  conda "$out_i/$y"

	input:
	file x from ref
  file y from canu_yaml

	output:
  file "hicanu/*.contigs.fasta" into canu
  file "hicanu/*.report" canu_report
  file "hicanu/*.unassembled.fasta" unassembled
  file "hicanu/*.layout.*" layouts

	script:
	"""
  memory=`echo "${task.memory}" | sed 's/[^0-9]//g'`
  memory="\${memory}G"
  canu \
  -p asm -d hicanu \
  genomeSize=${params.genomeSize} \
  -pacbio-hifi $x \
  useGrid=false \
  maxThreads=$task.cpus maxMemory=\$memory \
  merylThreads=$task.cpus merylMemory=\$memory merylConcurrency=1 \
  hapThreads=$task.cpus hapMemory=\$memory hapConcurrency=1 \
  batThreads=$task.cpus batMemory=\$memory batConcurrency=1 \
  gridOptionsJobName="hicanu"
	"""
}