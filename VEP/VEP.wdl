task VEP{
    File vcf
    String refseq="Yes"
    String assembly
    command <<<
        set -ex
        cd /var/data
        refseq=""
        if [ ${refseq} == "Yes" ]; then
            refseq="--refseq"
        fi
        if [ ${assembly} == "GRCh37" ]; then
            /opt/vep/src/ensembl-vep/vep -i ${vcf} --dir_cache /rdata/genedock/VEP/Anno_Database -o annotated.vcf $refseq --force_overwrite --everything --vcf -fork 4 --assembly GRCh37 --offline -plugin dbscSNV,/rdata/genedock/VEP/Anno_Database/dbscSNV.txt.gz --plugin SPIDEX,/rdata/genedock/VEP/Anno_Database/hg19_spidex.txt.gz  --fasta /rdata/genedock/VEP/Anno_Database/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --custom /rdata/genedock/VEP/Anno_Database/clinvar_20180603.vcf.gz,clinvar,vcf,exact,0,CLNDN,CLNSIG,CLNREVSTAT
        fi
        if [ ${assembly} == "GRCh38" ]; then
            /opt/vep/src/ensembl-vep/vep -i ${vcf} --dir_cache /rdata/genedock/VEP/Anno_Database -o annotated.vcf $refseq --force_overwrite --everything --vcf -fork 4 --assembly GRCh38 --offline --fasta /rdata/genedock/VEP/Anno_Database/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --custom /rdata/genedock/VEP/Anno_Database/clinvar_20201107_GRCh38.vcf.gz,clinvar,vcf,exact,0,CLNDN,CLNSIG,CLNREVSTAT
        fi
        
    >>>
    output{
        File out="annotated.vcf"
    }
    runtime {
        docker: "seqflow/vep:101"
        memory: "8G"
        disk: "200G"
        cpu: 4
    }
}


workflow VEP_annotation{
    File input_vcf
    String refseq="Yes"
    String assembly="GRCh37"
    call VEP{
        input:
            vcf=input_vcf,
            refseq=refseq,
            assembly=assembly,
    }
    
    output{
        VEP.out
    }
}
