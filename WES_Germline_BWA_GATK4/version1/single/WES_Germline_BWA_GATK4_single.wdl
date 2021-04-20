##WES
task soapnuke_MergeFq{
  input{
        Array[File] read1
        Array[File] read2
        String sample_id
        String adapter_read1
        String adapter_read2
        Float misMatch
        Float nRate
        Float lowQual
        Float qualRate
        Int qualSys
        }
        command {
                set -ex
                SOAPnuke=/bioapp/SOAPnuke/SOAPnuke-1.5.6
                mkdir ${sample_id}_SOAPnuke
                cd ${sample_id}_SOAPnuke
                cat ${sep=' ' read1} > ${sample_id}_merge_read1.fq.gz
                cat ${sep=' ' read2} > ${sample_id}_merge_read2.fq.gz
                if [ ${qualSys} == 2 ] ;then 
                $SOAPnuke filter -f ${adapter_read1} -r ${adapter_read2} -1 ${sample_id}_merge_read1.fq.gz -2 ${sample_id}_merge_read2.fq.gz -M ${default=1 misMatch} -n ${default=0.1 nRate} -l ${default=5 lowQual} -q ${default=0.5 qualRate}  -Q ${default=2 qualSys} -G -C ${sample_id}_clean.1.fq -D  ${sample_id}_clean.2.fq   
                else
                $SOAPnuke filter -f ${adapter_read1} -r ${adapter_read2} -1 ${sample_id}_merge_read1.fq.gz -2 ${sample_id}_merge_read2.fq.gz -M ${default=1 misMatch} -n ${default=0.1 nRate} -l ${default=5 lowQual} -q ${default=0.5 qualRate}  -Q ${qualSys} -C ${sample_id}_clean.1.fq -D  ${sample_id}_clean.2.fq 
                fi
                mkdir result
                perl /genedock/ref/script/soapnuke_stat.pl Basic_Statistics_of_Sequencing_Quality.txt Statistics_of_Filtered_Reads.txt > result/${sample_id}.filter.stat.xls
                perl /genedock/ref/script/filter_stat.pl -indir  `pwd`   -output ${sample_id}_output.xls
        }
        output{
                File clean_read1="/var/data/${sample_id}_SOAPnuke/${sample_id}_clean.1.fq.gz"
                File clean_read2="/var/data/${sample_id}_SOAPnuke/${sample_id}_clean.2.fq.gz"
                File stat="/var/data/${sample_id}_SOAPnuke/result/${sample_id}.filter.stat.xls"
                File out="/var/data/${sample_id}_SOAPnuke/${sample_id}_output.xls"
        }
        runtime {
                docker: "genedockdx/qc:1.4"
                memory: "8G"
                disk: "400G"
                cpu: 4
        }
}
task align_bwa{
        File read1
        File read2
        String sample_lb
        String sample_id
        String sample_pu
        String sample_sm
        String sample_pl
        String refname
        Int thread
        command <<<
                set -x
                if [ ${refname} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                fi
                if [ ${refname} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                fi
                if [ ${refname} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                fi
                if [ ${refname} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                fi
                bwa=/bioapp/bwa-0.7.17/bwa
                samtools=/bioapp/samtools-1.10/samtools
                sample_id=${sample_id}
                sample_lb="sample_lb"
                sample_pu=${sample_pu}
                sample_sm=${sample_sm}
                sample_pl=${sample_pl}
                $bwa mem -t ${thread} -M -R "@RG\tID:$sample_id\tPL:$sample_pl\tPU:$sample_pu\tLB:$sample_lb\tSM:$sample_sm" $reference ${read1} ${read2} >${sample_id}.sam
                $samtools sort -@  ${thread} -o ${sample_id}_${sample_sm}.sort.bam -O bam ${sample_id}.sam
                $samtools index -@ 8 ${sample_id}_${sample_sm}.sort.bam
        >>>
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 8
        }
        output{
                File SortBam="${sample_id}_${sample_sm}.sort.bam"
                File SortBai="${sample_id}_${sample_sm}.sort.bam.bai"
        }
}

### 去重、BQSR、callgvcf
task WES_MarkdupGvcf_GATK4{
        File bam
        File bai
        File bed
        String refname
        String sample_id
        Int stand_call_conf
        command <<<
                set -x
                if [ ${refname} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                        dir="/rdata/genedock/hg19_broad"
                        BR="1000G_phase1.indels.hg19.sites.vcf&Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&dbsnp_138.hg19.vcf"
                fi
                if [ ${refname} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                        dir="/rdata/genedock/hg38_broad"
                        BR="Mills_and_1000G_gold_standard.indels.hg38.vcf&dbsnp_138.hg38.vcf"
                fi
                if [ ${refname} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                        dir="/rdata/genedock/b37_broad"
                        BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
                fi
                if [ ${refname} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                        dir="/rdata/genedock/b37_broad"
                        BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
                fi
                mkdir -p /var/data/tmp
                gatk=/bioapp/gatk-4.1.4.1/gatk
                javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
                db_BR=`echo $BR|sed "s#\&# --known-sites \$dir/#g" |xargs -i echo  "--known-sites  $dir/{}"`
                ##markdup 
                $gatk --java-options "$javaOption" MarkDuplicates -I ${bam}  -O ${sample_id}.mkdup.bam -M file.mat --CREATE_INDEX true


                $gatk --java-options "$javaOption" BaseRecalibrator -R $reference -I ${sample_id}.mkdup.bam $db_BR -L ${bed} -O ${sample_id}_recal_data.table && \
                $gatk --java-options "$javaOption" ApplyBQSR -R $reference -bqsr ${sample_id}_recal_data.table -I ${sample_id}.mkdup.bam -L ${bed} -O bqsr.bam && \
                $gatk --java-options "$javaOption" HaplotypeCaller -R $reference -I bqsr.bam -L ${bed} -ERC GVCF -O ${sample_id}.g.vcf.gz  && \
                $gatk --java-options "$javaOption" GenotypeGVCFs -R $reference  -V ${sample_id}.g.vcf.gz -O ${sample_id}.raw.vcf.gz -stand-call-conf ${stand_call_conf}
        >>>
        output{
                File Gvcf="${sample_id}.g.vcf.gz"
                Float gvcf_size_in_GB = size("${sample_id}.g.vcf.gz", "G")
                File RawVcf="${sample_id}.raw.vcf.gz"
                File markdup_bam="${sample_id}.mkdup.bam"
                File markdup_bai="${sample_id}.mkdup.bai"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 8
        }
}
task WES_VcfFilter_GATK4{
        File vcf
        String sample_id
        String refname
        command <<<
                set -x
                if [ ${refname} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                fi
                if [ ${refname} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                fi
                if [ ${refname} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                fi
                if [ ${refname} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                fi
                gatk=/bioapp/gatk-4.1.4.1/gatk
                mkdir -p /var/data/tmp
                javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
                ### snp filter
                tabix ${vcf}
                $gatk --java-options "$javaOption" SelectVariants -R $reference -V ${vcf} -select-type SNP -O raw.snp.vcf.gz
                $gatk --java-options "$javaOption"  VariantFiltration -R $reference -V raw.snp.vcf.gz --filter-name "filter" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0||SOR > 3.0" -O ${sample_id}.flt.snp.vcf.gz
                ### indel filter
                $gatk --java-options "$javaOption" SelectVariants -R $reference -V ${vcf}  -select-type INDEL -O raw.indel.vcf.gz
                $gatk --java-options "$javaOption" VariantFiltration -R $reference -V raw.indel.vcf.gz --filter-name "filter" --filter-expression "QD < 2.0 || FS > 200.0|| SOR > 10.0"  -O ${sample_id}.flt.indel.vcf.gz
                $gatk  --java-options "$javaOption" SortVcf -I ${sample_id}.flt.snp.vcf.gz -I ${sample_id}.flt.indel.vcf.gz -O ${sample_id}.flt.vcf.gz
        >>>
        output{
                File FltSnp="${sample_id}.flt.snp.vcf.gz"
                File FltIndel="${sample_id}.flt.indel.vcf.gz"
                File FltVCF="${sample_id}.flt.vcf.gz"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "50G"
                cpu: 4
        }
}
task bamstat{
  File bam
  File bai
  File region
  String sample_id
  command <<< 
   set -ex
   mkdir ${sample_id}
   cd ${sample_id}
   perl /bioapp/Bam_StatScripts/markdup_bam_stat_4thread.pl -ad ${bam} -r ${region} -o output
   python /bioapp/Bam_StatScripts/merge_bam_stat.py output/*.information.xls > information.xls
   python /bioapp/Bam_StatScripts/merge_cumu.py output/*.cumu.xls > cumu.xls
   python /bioapp/Bam_StatScripts/merge_depth_frequency.py output/*.depth_frequency.xls > depth_frequency.xls
   cat output/*chrall.stat > chrall.stat
   rm -rf output
   cd ..
   tar zcvf ${sample_id}_BamStat.tar.gz ${sample_id}
   awk -F '[:]' '{print $1"\t"$NF}' ${sample_id}/information.xls |cut -f 1,3  > tmp.txt
   cat tmp.txt|grep -w -E 'Total_reads_num_in_bam|Mapped_reads_num|Mapping_rate|Duplication_rate|Mapped_bases_num|Mismatch_rate|Target_region_size|BaseNum_covered_on_target|Coverage_of_target_region|Average_sequencing_depth_on_target'|grep -v cigar > ${sample_id}_bamstat.txt
   grep Fraction_of_target_covered_with_at_least tmp.txt >> ${sample_id}_bamstat.txt
  >>>
 runtime {
        docker: "genedockdx/qc:1.4"
    memory: "16G"
        disk: "200G" 
        cpu: 8
    }
 output{
   File outstat="${sample_id}_BamStat.tar.gz"
   File out_txt="${sample_id}_bamstat.txt"
  }
}
task VEP{
    File vcf
    String refseq="Yes"
    String refname
    String sample_id
    command <<<
        set -ex
        cd /var/data
        refseq=""
        if [ ${refseq} == "Yes" ]; then
            refseq="--refseq"
        fi
        if [[ ${refname} == "b37" || ${refname} == "hg19" || ${refname} == "hs37d5" ]]; then
            /opt/vep/src/ensembl-vep/vep -i ${vcf} --dir_cache /rdata/genedock/VEP/Anno_Database -o ${sample_id}_annotated.vcf $refseq --force_overwrite --everything --vcf --hgvs --symbol --numbers  --domains --canonical --af --af_1kg  --fork 4 --assembly GRCh37 --offline -plugin dbscSNV,/rdata/genedock/VEP/Anno_Database/dbscSNV.txt.gz --plugin SPIDEX,/rdata/genedock/VEP/Anno_Database/hg19_spidex.txt.gz  --fasta /rdata/genedock/VEP/Anno_Database/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --custom /rdata/genedock/VEP/Anno_Database/clinvar_20180603.vcf.gz,clinvar,vcf,exact,0,CLNDN,CLNSIG,CLNREVSTAT --pick --pick_order canonical,tsl,biotype,rank,ccds,length
        elif [ ${refname} == "hg38" ]; then
            /opt/vep/src/ensembl-vep/vep -i ${vcf} --dir_cache /rdata/genedock/VEP/Anno_Database -o ${sample_id}_annotated.vcf $refseq --force_overwrite --everything --vcf --hgvs --symbol --numbers  --domains --canonical --af --af_1kg   -fork 4 --assembly GRCh38 --offline --fasta /rdata/genedock/VEP/Anno_Database/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --custom /rdata/genedock/VEP/Anno_Database/clinvar_20201107_GRCh38.vcf.gz,clinvar,vcf,exact,0,CLNDN,CLNSIG,CLNREVSTAT --pick --pick_order canonical,tsl,biotype,rank,ccds,length
        fi
        
    >>>
    output{
        File out="${sample_id}_annotated.vcf"
    }
    runtime {
        docker: "public/vep:101"
        memory: "16G"
        disk: "200G"
        cpu: 4
    }
}
task wes_vcf_json{
  File vcf="genedockdx:/home/admin/wdl_test/wes_test/test_again/call-MarkdupGvcf/sm.raw.vcf.gz"
  String sample_name="sm"
  String ref="b37"
  command <<<
    set -ex
                if [ ${ref} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                        BR="/rdata/genedock/hg19_broad/bsnp_138.hg19.vcf"
                fi
                if [ ${ref} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                        dbsnp="/rdata/genedock/hg38_broad/dbsnp_138.hg38.vcf"
                fi
                if [ ${ref} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                        dbsnp="/rdata/genedock/b37_broad/dbsnp_138.b37.vcf"
                fi
                if [ ${ref} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                        dbsnp="/rdata/genedock/b37_broad/dbsnp_138.b37.vcf"
                fi
    date
    picard=/bioapp/picard_2.23/picard.jar
    num=`zcat ${vcf} |grep -v '^##'|head -1|awk '{for(i=1;i<=NF;i++){if($i=="'${sample_name}'"){print i}}}'`
    zless ${vcf}|cut -f 1-9,$num |grep -w -v -E '0\/0|\.\/\.' > /var/data/${sample_name}.vcf
    bgzip ${sample_name}.vcf
    tabix -p vcf ${sample_name}.vcf.gz
    mkdir /var/data/${sample_name}_vcf_qc
    cp /bioapp/QC_script/vcf_list ${sample_name}_vcf_qc/vcf_list
    cd ${sample_name}_vcf_qc
    zcat /var/data/${sample_name}.vcf.gz|grep -v '#'|awk -F '\t' '{if(length($4)==1&&length($5)==1){x=split($8,tag,";");for(i=1;i<=x;i++){if(tag[i]~"QD"){qd=i}else if(tag[i]~"SOR"){sor=i}}}{print tag[qd]"\t"tag[sor]"\tQUAL="$6"\t"$1}}'|tr "=" "\t" |awk 'BEGIN{print "chr\tQD\tSOR\tQUAL"}{print $7"\t"$2"\t"$4"\t"$6}' > QD_SOR.txt
    cut -f 2 QD_SOR.txt > QD.txt
    cut -f 3 QD_SOR.txt > SOR.txt
    cut -f 1,4 QD_SOR.txt > QUAL.txt
    zcat /var/data/${sample_name}.vcf.gz|grep -v '#'|awk -F '\t' '{if(length($4)==1&&length($5)==1){print $0}}'|cut -f 9,10|awk -F '\t' '{x=split($1,tag,":");for(i=1;i<=x;i++){if(tag[i]=="AD"){ad=i}else if(tag[i]=="DP"){dp=i}};split($2,val,":");split(val[ad],ref,",");if(val[dp]!=0&&val[dp]~"[0-9]"){print (val[dp]-ref[1])/val[dp]}}' > vcf_depth.txt
    time java -jar $picard CollectVariantCallingMetrics INPUT=/var/data/${sample_name}.vcf.gz OUTPUT=${sample_name} DBSNP=$dbsnp
    cat ${sample_name}.variant_calling_detail_metrics|head -8|tail -n 2 > picard_vcf.txt
    bcftools stats -s- /var/data/${sample_name}.vcf.gz > ${sample_name}_bcftools.txt
    grep '^TSTV'  ${sample_name}_bcftools.txt|awk '{print "biTsTvRatio\t"$5}' > list.txt
    less ${sample_name}_bcftools.txt |grep '^PSC'|awk '{print "homoVariantsNum\t"$5"\nheteroVariantsNum\t"$6"\nhomoHeteroRatio\t"$5/$6}' >> list.txt
    grep ^IDD ${sample_name}_bcftools.txt|awk 'BEGIN{print "Length\tCount"}{print $3"\t"$4}' > INDEL_length.txt
    grep '^ST'  ${sample_name}_bcftools.txt|awk 'BEGIN{print "Type\tCount"}{print $3"\t"$4}' > ST.txt
    Rscript /bioapp/QC_script/WES_vcf_json.r
    cp VariantCallingQC.json /var/data
    rm ST.txt  picard_vcf.txt vcf_depth.txt list.txt INDEL_length.txt 
    cd /var/data
    wait 
    tar -zcvf ${sample_name}_vcf_qc.tar.gz ${sample_name}_vcf_qc
    mv /var/data/VariantCallingQC.json /var/data/${sample_name}_VariantCallingQC.json
    date
>>>
  runtime {
    docker: "genedockdx/qc:1.4"
    memory: "16G"
    disk: "100G" 
    cpu: 4
    }
  output{
    File picard_out="${sample_name}_vcf_qc.tar.gz"
    File VariantCall_json="/var/data/${sample_name}_VariantCallingQC.json"
}
}
task anno_json{
  File snp_vep
  File indel_vep
  String sample_id
    command <<<
    set -ex
    echo """3_prime_UTR_variant
5_prime_UTR_variant
downstream_gene_variant
frameshift_variant
inframe_insertion
intergenic_variant
intron_variant
missense_variant
non_coding_transcript_exon_variant
non_coding_transcript_variant
regulatory_region_variant
splice_acceptor_variant
splice_donor_variant
splice_region_variant
start_lost
stop_gained
synonymous_variant
TF_binding_site_variant
upstream_gene_variant""" > vep_list
    cat vep_list|while read f;do awk 'BEGIN{sum=0}{if($8~"'$f'"){sum++}}END{print "'$f'""\t"sum}' ${snp_vep} >>  /var/data/${sample_id}_SnpAnnoQC.tsv ;done
    cat vep_list|while read f;do awk 'BEGIN{sum=0}{if($8~"'$f'"){sum++}}END{print "'$f'""\t"sum}' ${indel_vep} >>  /var/data/${sample_id}_IndelAnnoQC.tsv ;done

    >>>
    output{
        File SnpAnnoQC="/var/data/${sample_id}_SnpAnnoQC.tsv"
        File IndelAnnoQC="/var/data/${sample_id}_IndelAnnoQC.tsv"
    }
    runtime {
        docker: "genedockdx/qc:1.4"
            memory: "4G"
            disk: "100G"
            cpu: 2
    }
}
workflow wf_WES{
    input{
    Array[File] lane_fq1=["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_1.fastq.gz"]
    Array[File] lane_fq2=["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_2.fastq.gz"]
    File bed="public:/demo-data/WES-Germline_NA12878_NGv3-LAB1360/NA12878-NGv3-LAB1360-A-regions.bed"
    String sample_pu="pu"
    String sample_id="sample"
    String sample_pl="illumina"
    Int stand_call_conf=30
    String adapter_read1
    String adapter_read2
    String sample_lb
    String refname="b37"
    Float misMatch=1
    Float nRate=0.1
    Float lowQual=5
    Float qualRate=0.5
    Int qualSys=2
    Int bwa_thread=8
    String refseq="Yes"
    }
    call soapnuke_MergeFq{
      input:
        read1=lane_fq1,
        read2=lane_fq2,
        sample_id=sample_id,
        adapter_read1=adapter_read1,
        adapter_read2=adapter_read2,
        misMatch=misMatch,
        lowQual=lowQual,
        nRate=nRate,
        qualRate=qualRate,
        qualSys=qualSys
    }
    call align_bwa{
      input:
        read1=soapnuke_MergeFq.clean_read1,
        read2=soapnuke_MergeFq.clean_read2,
        sample_id=sample_id,
        sample_pu=sample_pu,
        sample_sm=sample_id,
        sample_pl=sample_pl,
        sample_lb=sample_lb,
        refname=refname,
        thread=bwa_thread,
}
    call WES_MarkdupGvcf_GATK4{
      input:
        bam=align_bwa.SortBam,
        bai=align_bwa.SortBai,
        bed=bed,
        refname=refname,
        sample_id=sample_id,
        stand_call_conf=stand_call_conf,
    }
    call WES_VcfFilter_GATK4{
      input:
        vcf=WES_MarkdupGvcf_GATK4.RawVcf,
        sample_id=sample_id,
        refname=refname,
    }
    call VEP as VEP_SNP{
      input:
        vcf=WES_VcfFilter_GATK4.FltSnp,
        refseq=refseq,
        refname=refname,
        sample_id=sample_id+"SNP",
    }
    call VEP as VEP_INDEL{
      input:
        vcf=WES_VcfFilter_GATK4.FltIndel,
        refseq=refseq,
        refname=refname,
        sample_id=sample_id+"INDEL",
    }
    call bamstat{
      input:
        bam=WES_MarkdupGvcf_GATK4.markdup_bam,
        bai=WES_MarkdupGvcf_GATK4.markdup_bai,
        region=bed,
        sample_id=sample_id,
    }
    call wes_vcf_json{
      input:
      vcf=WES_VcfFilter_GATK4.FltVCF,
      sample_name=sample_id,
      ref=refname,
    }
    call anno_json{
      input:
        snp_vep=VEP_SNP.out,
        indel_vep=VEP_INDEL.out,
        sample_id=sample_id,
    }
    output{
    File  vep_snp=VEP_SNP.out
    File  vep_indel=VEP_INDEL.out
    File  fltsnp=WES_VcfFilter_GATK4.FltSnp
    File  fltindel=WES_VcfFilter_GATK4.FltIndel
    File  fltvcf=WES_VcfFilter_GATK4.FltVCF
    File  gvcf=WES_MarkdupGvcf_GATK4.Gvcf
    Float  gvcf_size=WES_MarkdupGvcf_GATK4.gvcf_size_in_GB
    File rawvcf=WES_MarkdupGvcf_GATK4.RawVcf
    File bam=WES_MarkdupGvcf_GATK4.markdup_bam
    File bai=WES_MarkdupGvcf_GATK4.markdup_bai
    File soapnuke_clean_read1=soapnuke_MergeFq.clean_read1
    File soapnuke_clean_read2=soapnuke_MergeFq.clean_read2
    File soapnuke_stat=soapnuke_MergeFq.stat
    File soapnuke_out=soapnuke_MergeFq.out
    File bam_stat=bamstat.outstat
    File bam_stat_txt=bamstat.out_txt
    File anno_snp=anno_json.SnpAnnoQC
    File anno_indel=anno_json.IndelAnnoQC
    File vcf_json=wes_vcf_json.VariantCall_json
    File vcf_picard=wes_vcf_json.picard_out
    }
}
