#WES_Somatic_BWA_GATK_Mutect2_ContainRef_stat
task GD_toolkit_mapping_ContainRef{
    File read1
    File read2
    File adaptor
    Int base_quality_type=33
    Int crop=10000
    Int head_crop=0
    Int trailing=3
    Int leading=3
    Int minlength=36
    String read_id
    String read_sm
    String read_lib="lib"
    String refname
    String sliding_window="4:15"
    String illumina_clip="2:30:10"
    String sample_name
    command <<<
        set -x
        if [ ${refname} == "hg19" ] ;then
            reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
        fi
        if [ ${refname} == "b37" ] ;then
            reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
        fi
        if [ ${refname} == "hg38" ] ;then
            reference=/rdata/genedock/hg38_broad/hg38.fasta
        fi
        echo "===== step1: trim and filter the fastq files ===="
        java -Xmx16g -jar /bioapp/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 ${read1} ${read2} forward_paired.fq.gz forward_unpaired.fq.gz reverse_paired.fq.gz reverse_unpaired.fq.gz ILLUMINACLIP:${adaptor}:${illumina_clip} SLIDINGWINDOW:${sliding_window} LEADING:${leading} TRAILING:${trailing} CROP:${crop} HEADCROP:${head_crop} MINLEN:${minlength}
        bwa mem -t 8 -M -R "@RG\tID:${read_id}\tSM:${read_sm}\tPL:illumina\tLB:${read_lib}" $reference forward_paired.fq.gz reverse_paired.fq.gz >pe.sam
        samtools sort -@ 8 -o sorted.bam -O bam pe.sam 
        samtools index sorted.bam
        mkdir -p /var/data/Fastq_Data
        basename1=${sample_name}"_sample_1.fq.gz"
        basename2=${sample_name}"_sample_2.fq.gz"
        mv ${read1} /var/data/Fastq_Data/$basename1
        mv ${read2} /var/data/Fastq_Data/$basename2
        /bioapp/Fastq_StatScripts/gd_qc_stat -i /var/data/Fastq_Data/$basename1,/var/data/Fastq_Data/$basename2 -a ${adaptor},${adaptor} -N 0.1 -q ${base_quality_type} -L 15 -p 0.5 -c -k -o /var/data/Fastq_Stat
        tar zcvf Fastq_Stat.tar.gz Fastq_Stat/ --exclude=Fastq_Stat/*.fq --exclude=Fastq_Stat/*.fastq
    >>>
    output{
        File bam="sorted.bam"
        File fqstat="Fastq_Stat.tar.gz"
    }
    runtime {
        docker: "public/i1mapping:1.0"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
}
task rmdup_bqsr_ContainRef_stat{
        File inbam
        File intervals
        String ref
        String sample_name
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
        inbam=${inbam}
        PICARD=/bioapp/picard-2.23.1/picard.jar
        GATK=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
        intervals=${intervals}
        db_BR=`echo "$BR"|sed "s#\&# -knownSites \$dir/#g" |xargs -i echo  "-knownSites $dir/{}"`
        echo "===== step0: picard MarkDuplicates ====="

        java -Xmx8G -jar $PICARD MarkDuplicates I=$inbam O=dups_marked.bam METRICS_FILE=out_dups_metrics.txt REMOVE_DUPLICATES=true 
        samtools index dups_marked.bam

        echo "===== step1: base recalibration ====="
        java  -Xmx8G -jar $GATK -T BaseRecalibrator -R $reference -I dups_marked.bam $db_BR -L  $intervals -o recal_data.table

        java -Xmx8G  -jar $GATK -T PrintReads -R $reference -I dups_marked.bam -L $intervals -BQSR recal_data.table  -o  ${sample_name}.recal.bam

        echo "===== step2: stat bam result =====" 
        perl /bioapp/Bam_StatScripts/markdup_bam_stat_4thread.pl -ad  dups_marked.bam -r $intervals -o Bam_Stat -plot
        tar zcvf Bam_Stat.tar.gz Bam_Stat/ --exclude=Bam_Stat/*.target.depth
        >>>
        output{
                File bam_stat="Bam_Stat.tar.gz"
        File recalBam="${sample_name}.recal.bam"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}
task realign{
        File normal
        File tumor
        File intervals
        String ref
        command <<<
                set -x
                if [ ${ref} == "hg19" ] ;then
                        reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
                        dir="/rdata/genedock/hg19_broad"
                        knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&1000G_phase1.indels.hg19.sites.vcf"
                fi
                if [ ${ref} == "hg38" ] ;then
                        reference=/rdata/genedock/hg38_broad/hg38.fasta
                        dir="/rdata/genedock/hg38_broad"
                        knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.hg38.vcf&1000G_phase1.indels.hg38.vcf"
                fi
                if [ ${ref} == "b37" ];then
                        reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
                        dir="/rdata/genedock/b37_broad"
                        knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.b37.vcf&1000G_phase1.indels.b37.vcf"
                fi
        db_IR=`echo "$knowsites_IndelRealigner"|sed "s#\&# -known \$dir/#g" |xargs -i echo  "-known $dir/{}"`

        samtools index ${normal}
        samtools index ${tumor}
        GATK=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
        intervals=${intervals}
        java -Xmx8G -jar $GATK -T RealignerTargetCreator -R $reference -I ${tumor} -I ${normal} $db_IR -L $intervals -o target_intervals.list
        java -Xmx8G -jar $GATK -T IndelRealigner -R $reference -I ${tumor} -I ${normal} -targetIntervals target_intervals.list $db_IR --nWayOut .realign.bam
        >>>
        output{
          File normalBam="normal.recal.realign.bam"
          File tumorBam="tumor.recal.realign.bam"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}
task MuTect2{
        File normal
        File tumor
        File cosmic_data
        File dbsnp_data
        File intervals
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
        samtools index ${normal}
        samtools index ${tumor}
        GATK=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
        java -jar $GATK -T MuTect2 -R $reference -I:tumor ${tumor} -I:normal ${normal} --dbsnp ${dbsnp_data} --cosmic ${cosmic_data} -L ${intervals} -o out.vcf
    >>>
        output{
          File outvcf="out.vcf"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}

task merge_fq_bam_stat{
    File fq_stat
    File bam_stat
    String sample_name
    command <<<
        cd /var/data
        set -x
        echo "====1. merge fq stat===="
        for i in *tar.gz;do tar zxvf $i;done
        python /bioapp/Fastq_StatScripts/merge_stat.py /var/data/Fastq_Stat/*.stat >fastq_statistics.txt
        python /bioapp/Fastq_StatScripts/merge_GC.py  /var/data/Fastq_Stat/raw_*.GC >raw_sample.GC
        python /bioapp/Fastq_StatScripts/merge_GC.py  /var/data/Fastq_Stat/clean_*.GC >clean_sample.GC
        python /bioapp/Fastq_StatScripts/merge_QD.py  /var/data/Fastq_Stat/raw_*.QD >raw_sample.QD
        python /bioapp/Fastq_StatScripts/merge_QD.py  /var/data/Fastq_Stat/clean_*.QD >clean_sample.QD
        python /bioapp/Fastq_StatScripts/merge_QM.py -raw fastq_statistics.txt /var/data/Fastq_Stat/raw_*.QM >raw_sample.QM
        python /bioapp/Fastq_StatScripts/merge_QM.py  -clean fastq_statistics.txt /var/data/Fastq_Stat/clean_*.QM >clean_sample.QM
        Rscript /bioapp/Fastq_StatScripts/plot_raw-clean_GC-QD-QM.R raw_sample.GC raw_sample.QD raw_sample.QM clean_sample.GC clean_sample.QD clean_sample.QM
        Rscript /bioapp/Fastq_StatScripts/plot_raw-clean_GC-QD-QM_pdf.R raw_sample.GC raw_sample.QD raw_sample.QM clean_sample.GC clean_sample.QD clean_sample.QM
        mkdir Fastq_Statistics_Out
        cp fastq_statistics.txt *.GC *.QD *.QM *.png *.pdf Fastq_Statistics_Out/
        tar zcvf ${sample_name}.Fastq_Statistics_Out.tar.gz Fastq_Statistics_Out/
        echo "====2. merge bam stat===="
        mkdir /var/data/Bam_Statistics_Out
        cd /var/data/Bam_Statistics_Out
        python /bioapp/Bam_StatScripts/merge_bam_stat.py /var/data/Bam_Stat/*.information.xls >mardup_bam_stat.xls
        python /bioapp/Bam_StatScripts/merge_cumu.py /var/data/Bam_Stat/*.cumu.xls >cumu.xls
        python /bioapp/Bam_StatScripts/merge_depth_frequency.py /var/data/Bam_Stat/*.depth_frequency.xls >depth_frequency.xls
        Rscript /bioapp/Bam_StatScripts/cumuPlot.R cumu.xls
        Rscript /bioapp/Bam_StatScripts/histPlot.R  depth_frequency.xls
        cd /var/data
        tar zcvf ${sample_name}.Bam_Statistics_Out.tar.gz Bam_Statistics_Out/
    >>>
    runtime {
        docker: "seqflow/genedock_wgs:1.0"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
    output{
        File fastq_stat="${sample_name}.Fastq_Statistics_Out.tar.gz"
        File bam_stats="${sample_name}.Bam_Statistics_Out.tar.gz"
    }
}
task cnvkit{
  File tumor_bam
  File normal_bam
  String refname
  File target_bed
  File refFlat
  String sample_name
  command {
  set -ex
  export PATH=/root/miniconda3/bin:/root/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
   if [ ${refname} == "hg19" ] ;then
     reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
     zless ${refFlat} > /var/data/refFlat.txt
   fi
   if [ ${refname} == "hg38" ] ;then
     reference=/rdata/genedock/hg38_broad/hg38.fasta
     zless ${refFlat} > /var/data/refFlat.txt
   fi
   if [ ${refname} == "b37" ];then
     reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
     zless ${refFlat}|sed 's/chr//g' > /var/data/refFlat.txt
   fi
   mkdir ${sample_name}_cnvkit_result
   cd ${sample_name}_cnvkit_result
   /root/miniconda3/bin/cnvkit.py access $reference -o access.ref.bed
   /root/miniconda3/bin/cnvkit.py antitarget ${target_bed} -g access.ref.bed -o my_antitargets.bed
   /root/miniconda3/bin/cnvkit.py autobin ${tumor_bam} ${normal_bam} -t ${target_bed} -g access.ref.bed --annotate /var/data/refFlat.txt
   target=`ls $(pwd)/*.target.bed|grep -v ${target_bed}`
   antitarget=`ls $(pwd)/*antitarget.bed` 
   # For each sample...
   /root/miniconda3/bin/cnvkit.py coverage ${tumor_bam} $target  -p 8 -o ${sample_name}.tumor.targetcoverage.cnn
   /root/miniconda3/bin/cnvkit.py coverage ${tumor_bam} $antitarget -p 8 -o ${sample_name}.tumor.antitargetcoverage.cnn
   /root/miniconda3/bin/cnvkit.py coverage ${normal_bam} $target  -p 8 -o ${sample_name}.normal.targetcoverage.cnn
   /root/miniconda3/bin/cnvkit.py coverage ${normal_bam} $antitarget  -p 8 -o ${sample_name}.normal.antitargetcoverage.cnn
   # With all normal samples...
   /root/miniconda3/bin/cnvkit.py reference *normal.*targetcoverage.cnn --fasta $reference -o my_reference.cnn
   /root/miniconda3/bin/cnvkit.py fix ${sample_name}.tumor.targetcoverage.cnn ${sample_name}.tumor.antitargetcoverage.cnn  my_reference.cnn -o ${sample_name}.cnr
   /root/miniconda3/bin/cnvkit.py segment ${sample_name}.cnr -o ${sample_name}.cns
    /root/miniconda3/bin/cnvkit.py call ${sample_name}.cns  -m  clonal -o ${sample_name}.call.cns
    /root/miniconda3/bin/cnvkit.py scatter ${sample_name}.cnr -s ${sample_name}.call.cns -o ${sample_name}-scatter.pdf
    /root/miniconda3/bin/cnvkit.py diagram ${sample_name}.cnr -s ${sample_name}.call.cns -o ${sample_name}-diagram.pdf
    rm  *.bed
    cd /var/data
    tar -zcvf ${sample_name}_cnvkit.tar.gz ${sample_name}_cnvkit_result/
    }
      runtime {
        docker: "genedockdx/cnv_filt:1.0"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
    output{
        File cnvkit_out="${sample_name}_cnvkit.tar.gz"
        File scatter_out="/var/data/${sample_name}_cnvkit_result/${sample_name}-scatter.pdf"
        File diagram_out="/var/data/${sample_name}_cnvkit_result/${sample_name}-diagram.pdf"
    }
}
task VEP{
    File vcf
    String refseq="Yes"
    String refname
    command <<<
        set -ex
        cd /var/data
        refseq=""
        if [ ${refseq} == "Yes" ]; then
            refseq="--refseq"
        fi
        if [[ ${refname} == "b37" || ${refname} == "hg19" ]]; then
            /opt/vep/src/ensembl-vep/vep -i ${vcf} --dir_cache /rdata/genedock/VEP/Anno_Database -o annotated.vcf $refseq --force_overwrite --everything --vcf --fork 4 --assembly GRCh37 --offline -plugin dbscSNV,/rdata/genedock/VEP/Anno_Database/dbscSNV.txt.gz --plugin SPIDEX,/rdata/genedock/VEP/Anno_Database/hg19_spidex.txt.gz  --fasta /rdata/genedock/VEP/Anno_Database/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --custom /rdata/genedock/VEP/Anno_Database/clinvar_20180603.vcf.gz,clinvar,vcf,exact,0,CLNDN,CLNSIG,CLNREVSTAT 
        fi
        if [ ${refname} == "hg38" ]; then
            /opt/vep/src/ensembl-vep/vep -i ${vcf} --dir_cache /rdata/genedock/VEP/Anno_Database -o annotated.vcf $refseq --force_overwrite --everything --vcf  -fork 4 --assembly GRCh38 --offline --fasta /rdata/genedock/VEP/Anno_Database/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --custom /rdata/genedock/VEP/Anno_Database/clinvar_20201107_GRCh38.vcf.gz,clinvar,vcf,exact,0,CLNDN,CLNSIG,CLNREVSTAT
        fi
        
    >>>
    output{
        File out="annotated.vcf"
    }
    runtime {
        docker: "public/vep:101"
        memory: "8G"
        disk: "200G"
        cpu: 4
    }
}
#WES Somatic tumor and normal pipeline
task fastqc{
        File read
        String name
        command {
                set -ex
                fastqc=/bioapp/FastQC_0_11_9/fastqc
                mkdir ${name}
                $fastqc -t 4 -o ${name} ${read}
                tar zcvf ${name}.tar.gz ${name}
                unzip ${name}/*zip
                cp */fastqc_data.txt ${name}.fastqc_data.txt
        }
        output{
                File Outfastqc="${name}.tar.gz"
                File fastqc="${name}.fastqc_data.txt"
                Float fq_size_in_GB=size("${read}","G")
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
      String type
      command{
          set -ex
          perl /bioapp/QC_script/getFqJson.pl ${Fq1Qc} ${Fq2Qc} > rawReadsQC.json
          mv rawReadsQC.json /var/data/rawReadsQC_${type}.json 
      }
      runtime {
                docker: "genedockdx/qc:1.4"
                memory: "4G"
                disk: "20G"
                cpu: 2
        }
        output{
                File OutFqJson="/var/data/rawReadsQC_${type}.json"
        }
}
task sample_json_somatic{
    Float normal_fq2_size
    Float normal_fq1_size
    Float tumor_fq2_size
    Float tumor_fq1_size
    Float normal_bam_size
    Float tumor_bam_size
    String sample_ID
    String gender="NULL"
    String workflow_name="WES_Somatic_tumor_only"
    File vcf
    command <<<
    set -ex 
    zless ${vcf} > sample.vcf
    bgzip  sample.vcf
    tabix -p vcf sample.vcf.gz
    echo ${sample_ID}|awk '{print "sampleID\t"$0}' >> sample.txt
    echo ${workflow_name}|awk '{print "sampleWorkflow\t"$0}' >> sample.txt
    echo ${gender}|awk '{print "sampleReportedSex\t"$0}' >> sample.txt
    echo ${normal_fq1_size}|awk '{print "normalRead1Size\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${normal_fq2_size}|awk '{print "normalRead2Size\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${tumor_fq1_size}|awk '{print "tumorRead1Size\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${tumor_fq2_size}|awk '{print "tumorRead2Size\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${normal_bam_size}|awk '{print "normalBAMSize\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo ${tumor_bam_size}|awk '{print "tumorBAMSize\t"sprintf("%.2f",$0)" GB"}' >> sample.txt
    echo "NULL"|awk '{print "normalGVCFSize\t"$0}' >> sample.txt
    echo "NULL"|awk '{print "tumorGVCFSize\t"$0}' >> sample.txt
    chr=`zless sample.vcf.gz|grep -v '#'|awk '{if($1~"chr"){print "chrX"}else{print "X"}}'`
    num=`zless sample.vcf.gz|grep -v $chr|wc -l`
    if [ $num == 0 ];then 
    bcftools stats -r $chr -s-  sample.vcf.gz |grep '^PSC'|awk '{print "sexHRatio\t"$6/$5}' >> sample.txt
    else echo "No chrX"; fi
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
task report{
    File config
    File AlignmentQC_tumor
    File AlignmentQC_normal
    File FastQC_tumor
    File FastQC_normal
    File SampleQC
    File VariantCallingQC
    String sample_id
    command <<<
    set -ex
    tar -zxvf ${config}
    mkdir config 
    mkdir json
    mv wes/Somatic_MetricsQC.xlsx  config/
    mv ${AlignmentQC_tumor} json/AlignmentQC.tumor.json
    mv ${AlignmentQC_normal} json/AlignmentQC.normal.json
    mv ${FastQC_tumor} json/RawReadQC.tumor.json
    mv ${FastQC_normal} json/RawReadQC.normal.json
    mv ${SampleQC} json/SampleQC.json
    mv ${VariantCallingQC} json/VariantCallingQC.json
    cd wes
    /bioapp/knitter.sh SinglePDF/SinglePDF.Rmd pdf SinglePDF/config_params.json
    #/bioapp/knitter.sh SingleHTML/SingleHTML.Rmd html SingleHTML/config_params.json
    #mv SingleHTML/SingleHTML.html /var/data/${sample_id}_Report.html
    mv SinglePDF/SinglePDF.pdf /var/data/${sample_id}_Report.pdf
    >>>
    runtime {
        docker: "genedockdx/report:1.0"
        memory: "4G"
        disk: "20G"
        cpu: 2
    }
    output{
        File report_pdf = "/var/data/${sample_id}_Report.pdf"
        #File report_html = "/var/data/${sample_id}_Report.html"
    }
}
task wes_vcf_json_somatic{
  File vcf
  String sample_id
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
    num=`zless ${vcf} |grep -v '^##'|head -1|awk '{for(i=1;i<=NF;i++){if($i=="TUMOR"){print i}}}'`
    zless ${vcf}|cut -f 1-9,$num |grep -w -v -E '0\/0|\.\/\.' > /var/data/${sample_id}.vcf
    bgzip ${sample_id}.vcf
    tabix -p vcf ${sample_id}.vcf.gz
    mkdir /var/data/${sample_id}_vcf_qc
    cp /bioapp/QC_script/vcf_list ${sample_id}_vcf_qc/vcf_list
    cd ${sample_id}_vcf_qc
    zcat /var/data/${sample_id}.vcf.gz|grep -v '#'|awk -F '\t' '{if(length($4)==1&&length($5)==1){print $0}}'|cut -f 9,10|awk -F '\t' '{x=split($1,tag,":");for(i=1;i<=x;i++){if(tag[i]=="AD"){ad=i}else if(tag[i]=="DP"){dp=i}};split($2,val,":");split(val[ad],ref,",");print ref[2]/(ref[1]+ref[2])}' > vcf_depth.txt
    time java -jar $picard CollectVariantCallingMetrics INPUT=/var/data/${sample_id}.vcf.gz OUTPUT=${sample_id} DBSNP=$dbsnp
    cat ${sample_id}.variant_calling_detail_metrics|head -8|tail -n 2 > picard_vcf.txt
    bcftools stats -s- /var/data/${sample_id}.vcf.gz > ${sample_id}_bcftools.txt
    grep '^TSTV'  ${sample_id}_bcftools.txt|awk '{print "biTsTvRatio\t"$5}' > list.txt
    less ${sample_id}_bcftools.txt |grep '^PSC'|awk '{print "homoVariantsNum\t"$5"\nheteroVariantsNum\t"$6"\nhomoHeteroRatio\t"$5/$6}' >> list.txt
    grep ^IDD ${sample_id}_bcftools.txt|awk 'BEGIN{print "Length\tCount"}{print $3"\t"$4}' > INDEL_length.txt
    grep '^ST'  ${sample_id}_bcftools.txt|awk 'BEGIN{print "Type\tCount"}{print $3"\t"$4}' > ST.txt
    sed -i '9,10s/^/#/g'  /bioapp/QC_script/WES_vcf_json.r
    sed -i '59s/^/#/g'  /bioapp/QC_script/WES_vcf_json.r
    Rscript /bioapp/QC_script/WES_vcf_json.r
    cp VariantCallingQC.json /var/data
    rm ST.txt  picard_vcf.txt vcf_depth.txt list.txt INDEL_length.txt 
    cd /var/data
    wait 
    tar -zcvf ${sample_id}_vcf_qc.tar.gz ${sample_id}_vcf_qc
    date
>>>
  runtime {
    docker: "genedockdx/qc:1.4"
    memory: "16G"
    disk: "100G" 
    cpu: 4
    }
  output{
    File picard_out="${sample_id}_vcf_qc.tar.gz"
    File VariantCall_json="/var/data/VariantCallingQC.json"
}
}
task bamstat{
        File inbam
        String type
        File intervals
        String ref
        command <<<
                set -ex
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
      mkdir Bam_Stat
      cd Bam_Stat
        perl /bioapp/Bam_StatScripts/markdup_bam_stat_4thread.pl -ad ${inbam} -r ${intervals} -o output
        python /bioapp/Bam_StatScripts/merge_bam_stat.py output/*.information.xls >markdup_bam_stat.xls
        python /bioapp/Bam_StatScripts/merge_cumu.py output/*.cumu.xls >cumu.xls
        python /bioapp/Bam_StatScripts/merge_depth_frequency.py output/*.depth_frequency.xls >depth_frequency.xls
        Rscript /bioapp/Bam_StatScripts/cumuPlot.R cumu.xls
        Rscript /bioapp/Bam_StatScripts/histPlot.R  depth_frequency.xls
        cd ..
        tar zcvf Bam_Stat.tar.gz Bam_Stat/ --exclude=Bam_Stat/*.target.depth
        awk '{if($1>=250){a+=$2;b+=$3}else{print $0}}END{print 250"\t"a"\t"b}' Bam_Stat/depth_frequency.xls >depth_frequency.xls
        head -251 Bam_Stat/cumu.xls >cumu.xls
        cat Bam_Stat/output/*chrall.stat >chrall.stat
        cat Bam_Stat/markdup_bam_stat.xls >information.xls
        awk -F '[:]' '{print $1"\t"$NF}' information.xls |cut -f 1,3  > json.txt
        Rscript /bioapp/QC_script/wes_gatk4_Alignment.r
        mv /var/data/Alignment.json /var/data/AlignmentQC_${type}.json
        >>>
        output{
                File bam_stat="Bam_Stat.tar.gz"
                File AlignmentQC="/var/data/AlignmentQC_${type}.json"
                Float bam_size_in_GB = size("${inbam}","G")
        }
        runtime {
                docker: "genedockdx/qc:1.4"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}



workflow WES_Somatic_BWA_GATK_Mutect2_ContainRef_stat{
    String sample_id
    File normal_fq1="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_normal_1000000_1.fq.gz"
    File normal_fq2="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_normal_1000000_2.fq.gz"
    File tumor_fq1="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_tumor_1000000_1.fq.gz"
    File tumor_fq2="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_tumor_1000000_2.fq.gz"
    File cosmic_data="public:/database/cosmic/b37_cosmic_v54_120711.vcf"
    File dbsnp_data="public:/reference/b37_broad/dbsnp_138.b37.vcf"
    File adaptor="public:/reference/adaptor_Trimmomatic/TruSeq3-PE-2.fa"
    File intervals
    File report_config="genedockdx:/home/admin/script/somatic_wes_pair.tar.gz"
    File refFlat
    String gender="NULL"
    String refname="b37"
    call GD_toolkit_mapping_ContainRef as FqFltMapNormal{
        input:
            read1=normal_fq1,
            read2=normal_fq2,
            adaptor=adaptor,
            sample_name="normal",
            read_sm="normal",
            read_id="normal",
            refname=refname
    }
    call GD_toolkit_mapping_ContainRef as FqFltMapTumor{
        input:
            read1=tumor_fq1,
            read2=tumor_fq2,
            adaptor=adaptor,
            sample_name="tumor",
            read_sm="tumor",
            read_id="tumor",
            refname=refname
    }
    call rmdup_bqsr_ContainRef_stat as rmbqstatNormal{
    input:
        inbam=FqFltMapNormal.bam,
        intervals=intervals,
        ref=refname,
        sample_name="normal"
    }
    call rmdup_bqsr_ContainRef_stat as rmbqstatTumor{
    input:
        inbam=FqFltMapTumor.bam,
        intervals=intervals,
        ref=refname,
        sample_name="tumor"
    }
    call realign{
    input:
        normal=rmbqstatNormal.recalBam,
        tumor=rmbqstatTumor.recalBam,
        intervals=intervals,
        ref=refname
    }
    call MuTect2{
    input:
        normal=realign.normalBam,
        tumor=realign.tumorBam,
        cosmic_data=cosmic_data,
        dbsnp_data=dbsnp_data,
        intervals=intervals,
        ref=refname
    }
   
    call merge_fq_bam_stat as FqBamQCNormal{
        input:
            fq_stat=FqFltMapNormal.fqstat,
            bam_stat=rmbqstatNormal.bam_stat,
            sample_name="normal"
    }
    call merge_fq_bam_stat as FqBamQCTumor{
        input:
            fq_stat=FqFltMapTumor.fqstat,
            bam_stat=rmbqstatTumor.bam_stat,
            sample_name="tumor"
    }
    call cnvkit{
    input:
    tumor_bam=realign.tumorBam,
    normal_bam=realign.normalBam,
    refname=refname,
    target_bed=intervals,
    sample_name="sample_id",
    refFlat=refFlat,
    }
    call VEP{
      input:
        vcf=MuTect2.outvcf,
        refseq="Yes",
        refname=refname
    }
    call bamstat as normal_bamstat{
    input:
      inbam=rmbqstatNormal.recalBam,
      intervals=intervals,
      type="normal",
      ref=refname,
    }
  call bamstat as tumor_bamstat{
    input:
      inbam=rmbqstatTumor.recalBam,
      type="tumor",
      intervals=intervals,
      ref=refname,
    }
  call fastqc as tumor_qc1{
    input:
      read=tumor_fq1,
      name="tumor_fq1qc",
  }
  call fastqc as tumor_qc2{
    input:
      read=tumor_fq2,
      name="tumor_fq2qc"
  }
  call getFqJson as tumor_getFqJson{
      input:
        Fq1Qc=tumor_qc1.fastqc,
        Fq2Qc=tumor_qc2.fastqc,
        type="tumor"
  }
  call fastqc as normal_qc1{
    input:
      read=normal_fq1,
      name="normal_fq1qc",
  }
  call fastqc as normal_qc2{
    input:
      read=normal_fq2,
      name="normal_fq2qc"
  }
  call getFqJson as normal_getFqJson{
      input:
        Fq1Qc=normal_qc1.fastqc,
        Fq2Qc=normal_qc2.fastqc,
        type="normal"
  }
  call wes_vcf_json_somatic{
    input:
      vcf=MuTect2.outvcf,
      sample_id=sample_id,
      ref=refname
  }
  call sample_json_somatic{
    input:
     normal_fq2_size=normal_qc2.fq_size_in_GB,
     normal_fq1_size=normal_qc1.fq_size_in_GB,
     tumor_fq2_size=tumor_qc2.fq_size_in_GB,
     tumor_fq1_size=tumor_qc1.fq_size_in_GB,
     normal_bam_size=normal_bamstat.bam_size_in_GB,
     tumor_bam_size=tumor_bamstat.bam_size_in_GB,
     sample_ID=sample_id,
     gender="NULL",
     workflow_name="WES_Somatic",
     vcf=MuTect2.outvcf
  }
  call report{
    input:
      config=report_config,
      AlignmentQC_tumor=tumor_bamstat.AlignmentQC,
      AlignmentQC_normal=normal_bamstat.AlignmentQC,
      FastQC_tumor=tumor_getFqJson.OutFqJson,
      FastQC_normal=normal_getFqJson.OutFqJson,
      SampleQC=sample_json_somatic.sampleqc,
      VariantCallingQC=wes_vcf_json_somatic.VariantCall_json,
      sample_id=sample_id
  }
    output{
    FqBamQCNormal.fastq_stat
    FqBamQCNormal.bam_stats
    FqBamQCTumor.fastq_stat
    FqBamQCTumor.bam_stats
    rmbqstatNormal.recalBam
    rmbqstatTumor.recalBam
    MuTect2.outvcf
    cnvkit.cnvkit_out
    cnvkit.scatter_out
    cnvkit.diagram_out
    VEP.out
    report.report_pdf
    }
}
