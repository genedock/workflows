task bamTocram{
    File bam
    String refname
    command <<<
        set -ex
        cd /var/data
        if [ ${refname} == "hg19" ]; then
            ref="/rdata/genedock/hg19_broad/ucsc.hg19.fasta"
            
        fi
        if [ ${refname} == "b37" ]; then
            ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
        fi
        
        samtools view -C -T $ref ${bam} > aln.cram
        
    >>>
    output{
        File cram="aln.cram"
    }
    runtime {
        docker: "seqflow/genedock_wgs:1.0"
        memory: "30G"
        disk: "1000G"
        cpu: 16
    }
}

workflow BAMtoCRAM{
    File input_bam
    String refname="hg19"
    call bamTocram{
        input:
            bam=input_bam,
            refname=refname
    }
    output{
        bamTocram.cram
    }
}
