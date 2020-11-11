### fq文件按照指定行数拆分
task split_fq{
    File fq
    Int? split_line
    command <<<
        set -x
        fq_name=$(basename ${fq} .gz)
        echo "fq name is $fq_name"
        pigz -dc ${fq}|/bioapp/split -d -l ${default=50000000 split_line} --filter='pigz > $FILE.gz' - $fq_name.
        rm ${fq}
    >>>
    runtime {
        docker: "seqflow/genedock_wgs:1.0"
        memory: "4G"
        disk: "400G"
        cpu: 2
    }
    output{
        Array[File] output_fq = glob("*.gz")
    }
}
### bwa mem比对到参考基因组
task align_bwa{
        Array[File] fq
        String sample_id
        String sample_lb
        String sample_pu
        String sample_sm
        String sample_pl
        String ref
        Int thread
        String read_name
        command <<<
                set -x
                if [ ${ref} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                fi
                if [ ${ref} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                fi
                if [ ${ref} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                fi
                bwa=/bioapp/bwa-0.7.17/bwa
                samtools=/bioapp/samtools-1.10/samtools
                sample_id=${sample_id}
                sample_lb=${sample_lb}
                sample_pu=${sample_pu}
                sample_sm=${sample_sm}
                sample_pl=${sample_pl}
                read_name=${read_name}
                $bwa mem -t ${thread} -M -R "@RG\tID:$sample_id\tPL:$sample_pl\tPU:$sample_pu\tLB:$sample_lb\tSM:$sample_sm" $reference ${sep=' ' fq} >$sample_id.sam
                $samtools sort -@  ${thread} -o $sample_id.$read_name.sort.bam -O bam $sample_id.sam
        >>>
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 8
        }
        output{
                File SortBam="${sample_id}.${read_name}.sort.bam"
        }
}
### 分块比对结果合并
task MergeBam_GATK4{
    Array[File] bam
    command <<<
        set -x
        mkdir -p /var/data/tmp
        gatk=/bioapp/gatk-4.1.4.1/gatk
        javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
        $gatk --java-options "$javaOption" MergeSamFiles -I ${sep=' -I ' bam}  -O merge.bam --USE_THREADING true --TMP_DIR /var/data/tmp
    >>>
        output{
                File MergeBam="merge.bam"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 4
        }
}
### 去重、BQSR、callgvcf
task WES_MarkdupGvcf_GATK4{
        File bam
        File bed
        String ref
        String sample_name
        Int stand_call_conf
        command <<<
                set -x
                if [ ${ref} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                        dir="/rdata/genedock/hg19_broad"
                        BR="1000G_phase1.indels.hg19.sites.vcf&Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&dbsnp_138.hg19.vcf"
                fi
                if [ ${ref} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                        dir="/rdata/genedock/hg38_broad"
                        BR="Mills_and_1000G_gold_standard.indels.hg38.vcf&dbsnp_138.hg38.vcf"
                fi
                if [ ${ref} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                        dir="/rdata/genedock/b37_broad"
                        BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
                fi
                mkdir -p /var/data/tmp
                gatk=/bioapp/gatk-4.1.4.1/gatk
                javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
                db_BR=`echo $BR|sed "s#\&# --known-sites \$dir/#g" |xargs -i echo  "--known-sites  $dir/{}"`
                $gatk --java-options "$javaOption" MarkDuplicates -I ${bam}  -O ${sample_name}.mkdup.bam -M ${sample_name}.mat --CREATE_INDEX true
                $gatk --java-options "$javaOption" BaseRecalibrator -R $reference -I ${sample_name}.mkdup.bam $db_BR -L ${bed} -O ${sample_name}_recal_data.table && \
                $gatk --java-options "$javaOption" ApplyBQSR -R $reference -bqsr ${sample_name}_recal_data.table -I ${sample_name}.mkdup.bam -L ${bed} -O bqsr.bam && \
                $gatk --java-options "$javaOption" HaplotypeCaller -R $reference -I bqsr.bam -L ${bed} -ERC GVCF -O ${sample_name}.g.vcf.gz  && \
                $gatk --java-options "$javaOption" GenotypeGVCFs -R $reference  -V ${sample_name}.g.vcf.gz -O ${sample_name}.raw.vcf.gz -stand-call-conf ${stand_call_conf}
        >>>
        output{
                File Gvcf="${sample_name}.g.vcf.gz"
                File RawVcf="${sample_name}.raw.vcf.gz"
                File MarkdupBam="${sample_name}.mkdup.bam"
                File MarkdupBai="${sample_name}.mkdup.bai"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 8
        }
}
### 变异过滤
task WES_VcfFilter_GATK4{
        File vcf
        String sample_name
        String ref
        command <<<
                set -x
                if [ ${ref} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                fi
                if [ ${ref} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                fi
                if [ ${ref} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                fi
                gatk=/bioapp/gatk-4.1.4.1/gatk
                mkdir -p /var/data/tmp
                javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
                ### snp filter
                tabix ${vcf}
                $gatk --java-options "$javaOption" SelectVariants -R $reference -V ${vcf} -select-type SNP -O raw.snp.vcf.gz
                $gatk --java-options "$javaOption"  VariantFiltration -R $reference -V raw.snp.vcf.gz --filter-name "filter" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0||SOR > 3.0" -O ${sample_name}.flt.snp.vcf.gz
                ### indel filter
                $gatk --java-options "$javaOption" SelectVariants -R $reference -V ${vcf}  -select-type INDEL -O raw.indel.vcf.gz
                $gatk --java-options "$javaOption" VariantFiltration -R $reference -V raw.indel.vcf.gz --filter-name "filter" --filter-expression "QD < 2.0 || FS > 200.0|| SOR > 10.0"  -O ${sample_name}.flt.indel.vcf.gz
        >>>
        output{
                File FltSnp="${sample_name}.flt.snp.vcf.gz"
                File FltIndel="${sample_name}.flt.indel.vcf.gz"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "50G"
                cpu: 4
        }
}

workflow wf_WES{
    File sample_fq1
    File sample_fq2
    File bed
    String sample_id="id"
    String sample_lb="lb"
    String sample_pu="pu"
    String sample_sm="sm"
    String sample_pl="illumina"
    String ref="hg19"
    Int bwa_thread=8
    Int stand_call_conf=30
    Boolean saw = true ### 是否将fq文件拆分去比对，是：true；否：false
    Int split_line=50000000  ### 指定拆分的行数，默认五千万行
    if(saw){
        call split_fq as split_fq1{
                input:
                    fq=sample_fq1,
                    split_line=split_line
            }
        call split_fq as split_fq2{
                input:
                    fq=sample_fq2,
                    split_line=split_line
            }
        scatter(Fq in transpose([split_fq1.output_fq,split_fq2.output_fq])){
            call align_bwa as align_saw{
                input:
                    fq=Fq,
                    sample_id=sample_id,
                    sample_lb=sample_lb,
                    sample_pu=sample_pu,
                    sample_sm=sample_sm,
                    sample_pl=sample_pl,
                    ref=ref,
                    thread=bwa_thread,
                    read_name=basename(Fq[0],".gz"),
            }
        }
        call MergeBam_GATK4 as MergeBam{
            input:
                bam=align_saw.SortBam
        }
    }
    if(!saw){
        call align_bwa as align{
            input:
                fq=[sample_fq1,sample_fq2],
                sample_id=sample_id,
                sample_lb=sample_lb,
                sample_pu=sample_pu,
                sample_sm=sample_sm,
                sample_pl=sample_pl,
                ref=ref,
                thread=bwa_thread,
                read_name=sample_sm
        }
    }
    call WES_MarkdupGvcf_GATK4 as MarkdupGvcf{
        input:
            bam=select_first([MergeBam.MergeBam,align.SortBam]),
            bed=bed,
            sample_name=sample_sm,
            ref=ref,
            stand_call_conf=stand_call_conf
    }
    call WES_VcfFilter_GATK4 as VcfFilter{
        input:
            vcf=MarkdupGvcf.RawVcf,
            sample_name=sample_sm,
            ref=ref,
    }
    output{
        File MarkdupBam=MarkdupGvcf.MarkdupBam
        File MarkdupBai=MarkdupGvcf.MarkdupBai
        File Gvcf=MarkdupGvcf.Gvcf
        File Rawvcf=MarkdupGvcf.RawVcf
        File FltSnp=VcfFilter.FltSnp
        File FltIndel=VcfFilter.FltIndel
    }
}