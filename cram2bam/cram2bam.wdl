task cramTobam{
    File cram
    command <<<
        set -ex
        cd /var/data
        samtools view -bS aln.cram >out.bam 
        
    >>>
    output{
        File bam="out.bam"
    }
    runtime {
        docker: "seqflow/genedock_wgs:1.0"
        memory: "30G"
        disk: "1000G"
        cpu: 16
    }
}

workflow CRAMtoBAM{
    File input_cram
    call cramTobam{
        input:
            cram=input_cram
    }
    output{
        cramTobam.bam
    }
}
