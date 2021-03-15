### fq文件按照指定行数拆分
task split_fq{
    Array[File] fq
    Int? compression=4
    Int? thread=4
    String sample
    command <<<
    set -ex
    date "+%G-%m-%d %H:%M:%S"
    /bioapp/fastp -i ${sep=' -I ' fq} --split_by_lines 50000000 -A -G -Q -L --thread ${thread} --compression ${compression} -o ${sample}.1.new.gz -O ${sample}.2.new.gz
    date "+%G-%m-%d %H:%M:%S"
    >>>
    runtime {
        docker: "seqflow/genedock_wgs:1.1"
        memory: "4G"
        disk: "400G"
        cpu: 2
    }
    output{
        Array[File] read1_fq = glob("*.1.new.gz")
        Array[File] read2_fq = glob("*.2.new.gz")
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
                if [ ${ref} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                fi
                bwa=/bioapp/bwa-0.7.17/bwa
                samtools=/bioapp/samtools-1.10/samtools
                sample_id=${sample_id}
                sample_lb=${sample_lb}
                sample_pu=${sample_pu}
                sample_sm=${sample_sm}
                sample_pl=${sample_pl}
                read_name=${read_name}
                $bwa mem -t ${thread} -M -R "@RG\tID:$sample_id\tPL:$sample_pl\tPU:$sample_pu\tLB:$sample_lb\tSM:$sample_sm" $reference ${sep=' ' fq} >${sample_id}.sam
                $samtools sort -@  ${thread} -o ${sample_id}.${read_name}.sort.bam -O bam ${sample_id}.sam
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
###不分块比对
task align{
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
                if [ ${ref} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                fi
                bwa=/bioapp/bwa-0.7.17/bwa
                samtools=/bioapp/samtools-1.10/samtools
                gatk=/bioapp/gatk-4.1.4.1/gatk
                javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
                sample_id=${sample_id}
                sample_lb=${sample_lb}
                sample_pu=${sample_pu}
                sample_sm=${sample_sm}
                sample_pl=${sample_pl}
                read_name=${read_name}
                $bwa mem -t ${thread} -M -R "@RG\tID:$sample_id\tPL:$sample_pl\tPU:$sample_pu\tLB:$sample_lb\tSM:$sample_sm" $reference ${sep=' ' fq} >${sample_id}.sam
                $samtools sort -@  ${thread} -o ${sample_id}.${read_name}.sort.bam -O bam ${sample_id}.sam
                $gatk --java-options "$javaOption" MarkDuplicates -I ${sample_id}.${read_name}.sort.bam  -O ${sample_sm}.mkdup.bam -M file.mat --CREATE_INDEX true
        >>>
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 8
        }
        output{
                File SortBam="${sample_sm}.mkdup.bam"
                File SortBai="${sample_sm}.mkdup.bai"
        }
}

### 分块比对结果合并
task MergeBam_GATK4{
    Array[File?] bam
    String sample_name
    command <<<
        set -x
        mkdir -p /var/data/tmp
        gatk=/bioapp/gatk-4.1.4.1/gatk
        javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
        $gatk --java-options "$javaOption" MergeSamFiles -I ${sep=' -I ' bam}  -O merge.bam --USE_THREADING true --TMP_DIR /var/data/tmp
        $gatk --java-options "$javaOption" MarkDuplicates -I merge.bam  -O ${sample_name}.mkdup.bam -M file.mat --CREATE_INDEX true
    >>>
        output{
                File MergeBam="${sample_name}.mkdup.bam"
                File MergeBai="${sample_name}.mkdup.bai"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "100G"
                cpu: 4
        }
}
task mergebam{
    Array[File] bam
    String sample_tmp_name
    command <<<
        set -x
        mkdir -p /var/data/tmp
        gatk=/bioapp/gatk-4.1.4.1/gatk
        javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
        $gatk --java-options "$javaOption" MergeSamFiles -I ${sep=' -I ' bam}  -O ${sample_tmp_name}.mkdup.bam --USE_THREADING true --TMP_DIR /var/data/tmp
    >>>
        output{
                File MergeBam="${sample_tmp_name}.mkdup.bam"
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
        File bai
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
                if [ ${ref} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
                        dir="/rdata/genedock/b37_broad"
                        BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
                fi
                mkdir -p /var/data/tmp
                gatk=/bioapp/gatk-4.1.4.1/gatk
                javaOption="-Xmx8g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:-UseGCOverheadLimit"
                db_BR=`echo $BR|sed "s#\&# --known-sites \$dir/#g" |xargs -i echo  "--known-sites  $dir/{}"`
                $gatk --java-options "$javaOption" BaseRecalibrator -R $reference -I ${bam} $db_BR -L ${bed} -O ${sample_name}_recal_data.table && \
                $gatk --java-options "$javaOption" ApplyBQSR -R $reference -bqsr ${sample_name}_recal_data.table -I ${bam} -L ${bed} -O bqsr.bam && \
                $gatk --java-options "$javaOption" HaplotypeCaller -R $reference -I bqsr.bam -L ${bed} -ERC GVCF -O ${sample_name}.g.vcf.gz  && \
                $gatk --java-options "$javaOption" GenotypeGVCFs -R $reference  -V ${sample_name}.g.vcf.gz -O ${sample_name}.raw.vcf.gz -stand-call-conf ${stand_call_conf}
        >>>
        output{
                File Gvcf="${sample_name}.g.vcf.gz"
                Float gvcf_size_in_GB = size("${sample_name}.g.vcf.gz", "G")
                File RawVcf="${sample_name}.raw.vcf.gz"
                File MarkdupBam="${sample_name}.mkdup.bam"
                Float bam_size_in_GB = size("${sample_name}.mkdup.bam","G")
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
                if [ ${ref} == "hs37d5" ];then
                        reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
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
task bamstat{
  File bam
  File bai
  File region
  command <<< 
   set -ex
   mkdir Bam_Statistics_Out
   cd Bam_Statistics_Out
   perl /bioapp/Bam_StatScripts/markdup_bam_stat_4thread.pl -ad ${bam} -r ${region} -o output
   python /bioapp/Bam_StatScripts/merge_bam_stat.py output/*.information.xls >markdup_bam_stat.xls
   python /bioapp/Bam_StatScripts/merge_cumu.py output/*.cumu.xls >cumu.xls
   python /bioapp/Bam_StatScripts/merge_depth_frequency.py output/*.depth_frequency.xls >depth_frequency.xls
   Rscript /bioapp/Bam_StatScripts/cumuPlot.R cumu.xls
   Rscript /bioapp/Bam_StatScripts/histPlot.R  depth_frequency.xls
   cd ..
   tar zcvf BamStat.tar.gz Bam_Statistics_Out
   awk '{if($1>=250){a+=$2;b+=$3}else{print $0}}END{print 250"\t"a"\t"b}' Bam_Statistics_Out/depth_frequency.xls >depth_frequency.xls
   head -251 Bam_Statistics_Out/cumu.xls >cumu.xls
   cat Bam_Statistics_Out/output/*chrall.stat >chrall.stat
   cat Bam_Statistics_Out/markdup_bam_stat.xls >information.xls
   awk -F '[:]' '{print $1"\t"$NF}' information.xls |cut -f 1,3   > json.txt
   Rscript /bioapp/QC_script/wes_gatk4_Alignment.r
  >>>
 runtime {
        docker: "genedockdx/qc:1.4"
    memory: "16G"
        disk: "200G" 
        cpu: 8
    }
 output{
   File outstat="BamStat.tar.gz"
   File AlignmentQC="/var/data/Alignment.json"
  }
}
task sample_json{
    Float fq2_size
    Float fq1_size
    Float bam_size
    Float gvcf_size
    String sample_ID
    String gender
    String workflow_name="WES_GATK4"
    File vcf
    command <<<
    set -ex 
    tabix -p vcf ${vcf}
    echo ${sample_ID}|awk '{print "sampleID\t"$0}' >> sample.txt
    echo ${workflow_name}|awk '{print "sampleWorkflow\t"$0}' >> sample.txt
    echo ${gender}|awk '{print "sampleReportedSex\t"$0}' >> sample.txt
    echo ${fq1_size}|awk '{print "sampleRead1Size\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${fq2_size}|awk '{print "sampleRead2Size\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${bam_size}|awk '{print "sampleBAMSize\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${gvcf_size}|awk '{print "sampleGVCFSize\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    bcftools stats -r X -s- ${vcf} |grep '^PSC'|awk '{print "sexHRatio\t"$6/$5}' >> sample.txt
    echo  'library(jsonlite)
    sample=read.table("/var/data/sample.txt",sep="\t",header=F,stringsAsFactors=F)
    list=c();tmp=c();for(i in 1:dim(sample)[1]){tmp=list(sample[i,2]);names(tmp)=sample[i,1];list=c(list,tmp)}
    cat(toJSON(list,pretty=T,auto_unbox=T),file="/var/data/SampleQC.json",append=F)' > script.r
    Rscript script.r
    >>>
    runtime {
        docker: "genedockdx/qc:1.4"
        memory: "4G"
        disk: "20G"
        cpu: 2
    }
    output{
      File sampleqc="/var/data/SampleQC.json"
    }
}
task fastqc{
        Array[File] read
        String name
        command {
                set -ex
                fastqc=/bioapp/FastQC_0_11_9/fastqc
                mkdir ${name}
                cat ${sep=' ' read} > read.fq.gz
                $fastqc -t 4 -o ${name} read.fq.gz
                tar zcvf ${name}.tar.gz ${name}
                unzip ${name}/*zip
                cp */fastqc_data.txt ${name}.fastqc_data.txt
        }
        output{
                File Outfastqc="${name}.tar.gz"
                File fastqc="${name}.fastqc_data.txt"
                Float fq_size_in_GB=size("read.fq.gz","G")
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "8G"
                disk: "400G"
                cpu: 4
        }
}
task getFqJson{
      File Fq1Qc
      File Fq2Qc
      command{
          set -ex
          perl /bioapp/QC_script/getFqJson.pl ${Fq1Qc} ${Fq2Qc} >rawReadsQC.json 
      }
      runtime {
                docker: "genedockdx/qc:1.4"
                memory: "4G"
                disk: "20G"
                cpu: 2
        }
        output{
                File OutFqJson="rawReadsQC.json"
        }
}
task report{
    File config
    File AlignmentQC
    File FastQC
    File SampleQC
    File VariantCallingQC
    String sample_name
    command <<<
    set -ex
    tar -zxvf ${config}
    mkdir config
    mv wes/MetricsQC.xlsx  config/
    mv ${AlignmentQC} wes/AlignmentQC.json
    mv ${FastQC} wes/RawReadQC.json
    mv ${SampleQC} wes/SampleQC.json
    mv ${VariantCallingQC} wes/VariantCallingQC.json
    cd wes
    /bioapp/knitter.sh SinglePDF/SinglePDF.Rmd pdf SinglePDF/config_params.json
    /bioapp/knitter.sh SingleHTML/SingleHTML.Rmd html SingleHTML/config_params.json
    mv SingleHTML/SingleHTML.html /var/data/${sample_name}_Report.html
    mv SinglePDF/SinglePDF.pdf /var/data/${sample_name}_Report.pdf
    >>>
    runtime {
        docker: "genedockdx/report:1.0"
        memory: "4G"
        disk: "20G"
        cpu: 2
    }
    output{
        File report_pdf = "/var/data/${sample_name}_Report.pdf"
        File report_html = "/var/data/${sample_name}_Report.html"
    }

}
task wes_vcf_json{
  #File vcf="genedockdx:/home/admin/WDL/zx/test1/call-CombFltcall-combine_filter/out.raw.vcf.gz"
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
    File VariantCall_json="/var/data/VariantCallingQC.json"
}
}


workflow wf_WES{
    Array[File] lane_fq1=["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_1.fastq.gz"]
    Array[File]lane_fq2=["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_2.fastq.gz"]
    File report_config="genedockdx:/home/admin/script/WES_report.tar.gz"
    File bed="public:/demo-data/WES-Germline_NA12878_NGv3-LAB1360/NA12878-NGv3-LAB1360-A-regions.bed"
    Array[String] lane_id=["L1"]
    String sample_lb="lb"
    String sample_pu="pu"
    String sample_sm="sm"
    String sample_pl="illumina"
    String ref="hs37d5"
    String gender="Male"
    Int? compression=4
    Int? thread=4
    Int bwa_thread=8
    Int stand_call_conf=30
    Boolean saw = false ### 是否将fq文件拆分去比对，是：true；否：false
    Int split_line=50000000  ### 指定拆分的行数，默认五千万行

    
    scatter(worker in range(length(lane_id))){
      Array[File] FQ = transpose([lane_fq1,lane_fq2])[worker]
      String lane = lane_id[worker]
    if(saw){
        call split_fq{
            input:
                fq=FQ,
                compression=compression,
                thread=thread,
                sample=sample_sm,
        }
        scatter(Fq in transpose([split_fq.read1_fq,split_fq.read2_fq])){
            call align_bwa as align_saw{
                input:
                    fq=Fq,
                    sample_id=lane,
                    sample_lb=sample_lb,
                    sample_pu=sample_pu,
                    sample_sm=sample_sm,
                    sample_pl=sample_pl,
                    ref=ref,
                    thread=bwa_thread,
                    read_name=basename(Fq[0],".gz"),
            }
        }
        call mergebam{
            input:
            bam=align_saw.SortBam,
            sample_tmp_name=basename(align_saw.SortBam[0],".bam")
        }
        }
    if(!saw){
        call align_bwa{
            input:
                fq=FQ,
                sample_id=lane,
                sample_lb=sample_lb,
                sample_pu=sample_pu,
                sample_sm=sample_sm,
                sample_pl=sample_pl,
                ref=ref,
                thread=bwa_thread,
                read_name=sample_sm
        }
    }
  }
    if(saw){
          call MergeBam_GATK4 as MergeBam{
            input:
                bam=mergebam.MergeBam,
                sample_name=sample_sm
        }
    }
    if(!saw){
          call MergeBam_GATK4 as MergeBam_2{
            input:
                bam=align_bwa.SortBam,
                sample_name=sample_sm
        }
    }
    call WES_MarkdupGvcf_GATK4 as MarkdupGvcf{
        input:
            bam=select_first([MergeBam.MergeBam,MergeBam_2.MergeBam]),
            bai=select_first([MergeBam.MergeBai,MergeBam_2.MergeBai]),
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

  call fastqc as qc1{
    input:
      read=lane_fq1,
      name="fq1qc",
  }
  call fastqc as qc2{
    input:
      read=lane_fq2,
      name="fq2qc"
  }
  call getFqJson{
      input:
        Fq1Qc=qc1.fastqc,
        Fq2Qc=qc2.fastqc,
  }
  call bamstat{
    input:
      bam=select_first([MergeBam.MergeBam,MergeBam_2.MergeBam]),
      bai=select_first([MergeBam.MergeBai,MergeBam_2.MergeBai]),
      region=bed,
  }
  call wes_vcf_json{
    input:
      vcf=MarkdupGvcf.RawVcf,
      sample_name=sample_sm,
      ref=ref
  }
  call sample_json{
    input:
        fq2_size=qc2.fq_size_in_GB,
        fq1_size=qc1.fq_size_in_GB,
        bam_size=MarkdupGvcf.bam_size_in_GB,
        gvcf_size=MarkdupGvcf.gvcf_size_in_GB,
        sample_ID=sample_sm,
        gender=gender,
        vcf=MarkdupGvcf.RawVcf,
  }
call report{
  input:
    config=report_config,
    AlignmentQC=bamstat.AlignmentQC,
    FastQC=getFqJson.OutFqJson,
    SampleQC=sample_json.sampleqc,
    VariantCallingQC=wes_vcf_json.VariantCall_json,
    sample_name=sample_sm
}

    output{
        File MarkdupBam=MarkdupGvcf.MarkdupBam
        File MarkdupBai=MarkdupGvcf.MarkdupBai
        File Gvcf=MarkdupGvcf.Gvcf
        File Rawvcf=MarkdupGvcf.RawVcf
        File FltSnp=VcfFilter.FltSnp
        File FltIndel=VcfFilter.FltIndel
        File report_pdf=report.report_pdf
        File report_html=report.report_html
    }
}
