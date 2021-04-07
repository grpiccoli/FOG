#!/usr/bin/env nextflow
Channel.fromPath("$params.ref", type: 'file')
.buffer(size:1)
.set{
    ref;
}

out_i="$params.outdir/indexes"
out_asm="$params.outdir/2_assemblers"

process hifiasm {
	tag "hifiasm.$x"
    container "$params.bio/hifiasm:0.13--h8b12597_0"
    publishDir "$out_asm/hifiasm"
    cache 'lenient'

	input:
	file x from ref

	output:
    file "*tg.gfa" into hifiasm_o
    file "*log"

    //errorStrategy { task.exitStatus=0 ? 'ignore' : 'terminate' }

	script:
	"""
    hifiasm=$workflow.launchDir/$out_asm/hifiasm
    gigs=`echo "$task.memory" | sed 's/[^0-9]//g'`
    if [[ \$gigs > 23 ]];
    then
        hifiasm -o hifiasm.asm -t $task.cpus $x 2> hifiasm.asm.log
    else
        echo "unsufficient memory: \$gigs GB, minimum required 24 GB";
        exit 1; 
    fi
	"""
}

process post_hifiasm {
	tag "post_hifiasm.$x"
    publishDir "$out_asm/hifiasm"

	input:
	file x from hifiasm_o

	output:
	file "*.fasta" into hifiasm

    script:
    """
    IFS=" " read -r -a array <<< "$x"
    for e in "\${array[@]}"
    do
        case \$e in *r_utg*)
            awk '/^S/{print ">"\$2;print \$3}' $x > hifiasm_raw.fasta
            ;;
        *p_utg*)
            awk '/^S/{print ">"\$2;print \$3}' $x > hifiasm_processed.fasta
            ;;
        *p_ctg*)
            awk '/^S/{print ">"\$2;print \$3}' $x > hifiasm_primary.fasta
            ;;
        *a_ctg*)
            awk '/^S/{print ">"\$2;print \$3}' $x > hifiasm_alternate.fasta
            ;;
        esac
    done
    """
}
