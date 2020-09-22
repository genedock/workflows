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
    String read_id="sample"
    String read_sm="sample"
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
task Sentieon_rmdup_realign_bqsr_hc_ContainRef_stat{
    File intervals
    File sort_bam
    Int call_conf=30
    Int emit_conf=30
    Int thread=8
    String refname
    String variant_caller_algo="Genotyper"
    String sample_name
    command <<<
        export SENTIEON_LICENSE=node3.genedock.com:16527
        date "+%G-%m-%d %H:%M:%S"
        set -x
        cd /var/data
        if [ ${refname} == "hg19" ]; then
            ref="/rdata/genedock/hg19_broad/ucsc.hg19.fasta"
            dir="/rdata/genedock/hg19_broad"
            knowsites_BaseRecalibrator="dbsnp_138.hg19.vcf&Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&1000G_phase1.indels.hg19.sites.vcf"
            knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&1000G_phase1.indels.hg19.sites.vcf"
            knowsites_dbsnp="dbsnp_138.hg19.vcf"
        fi
        if [ ${refname} == "b37" ]; then
            ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
            dir="/rdata/genedock/b37_broad"
            knowsites_BaseRecalibrator="dbsnp_138.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&1000G_phase1.indels.b37.vcf"
            knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.b37.vcf&1000G_phase1.indels.b37.vcf"
            knowsites_dbsnp="dbsnp_138.b37.vcf"
        fi
        if [ ${variant_caller_algo} == "Genotyper" ]; then
            algo="--algo Genotyper --var_type BOTH "
        elif [ ${variant_caller_algo} == "Haplotyper" ]; then
            algo="--algo Haplotyper"
        fi
        db_IR=`echo $knowsites_IndelRealigner|sed "s#\&# -k \$dir/#g" |xargs -i echo  "-k $dir/{}"`
        db_BR=`echo $knowsites_BaseRecalibrator|sed "s#\&# -k \$dir/#g" |xargs -i echo  "-k $dir/{}"`
        db_snp=`echo $knowsites_dbsnp|sed "s#\&# -d \$dir/#g" |xargs -i echo  "-d $dir/{}"`
        interval=${intervals}
        #sentieon util vcfindex $dir/*.vcf
        date "+%G-%m-%d %H:%M:%S"
        echo "====1. Index sort_bam===="
        sentieon util index ${sort_bam}
        date "+%G-%m-%d %H:%M:%S"
        echo "====2. Metrics===="
        sentieon driver -r $ref -t ${thread} -i  ${sort_bam} --interval $interval --algo MeanQualityByCycle smoke_mq_metrics.txt --algo QualDistribution smoke_qd_metrics.txt --algo GCBias --summary smoke_gc_summary.txt smoke_gc_metrics.txt --algo AlignmentStat smoke_aln_metrics.txt --algo InsertSizeMetricAlgo smoke_is_metrics.txt
        sentieon plot metrics -o metrics_report.pdf gc=smoke_gc_metrics.txt qd=smoke_qd_metrics.txt mq=smoke_mq_metrics.txt isize=smoke_is_metrics.txt
        date "+%G-%m-%d %H:%M:%S"
        echo "====3. Remove Duplicate Reads===="
        sentieon driver  -t ${thread} -i ${sort_bam} --algo LocusCollector --fun score_info smoke_score.txt
        sentieon driver  -t ${thread} -i ${sort_bam} --algo Dedup --score_info smoke_score.txt --metrics smoke_dedup_metrics.txt deduped.bam
        rm ${sort_bam}
        date "+%G-%m-%d %H:%M:%S"
        echo "====4. Indel realigner===="
        sentieon driver -r $ref -t ${thread} -i deduped.bam --algo Realigner $db_IR  --interval_list $interval  realigned.bam
        date "+%G-%m-%d %H:%M:%S"
        echo "====5. BQSR===="
        sentieon driver -r $ref -t ${thread} -i realigned.bam --interval $interval  --algo QualCal $db_BR smoke_recal_data.table
        sentieon driver -r $ref -t ${thread} -i realigned.bam --interval $interval -q smoke_recal_data.table --algo QualCal $db_BR smoke_recal_data.table.post
        sentieon driver -t ${thread} --algo QualCal   --plot --before smoke_recal_data.table --after smoke_recal_data.table.post smoke_recal_data.csv
        sentieon plot QualCal -o bqsr_report.pdf smoke_recal_data.csv
        date "+%G-%m-%d %H:%M:%S"
        echo "====6b. UG Variant caller===="
        sentieon driver -r $ref -t ${thread} -i realigned.bam --interval $interval  -q smoke_recal_data.table $algo $db_snp --emit_conf=${emit_conf} --call_conf=${call_conf} output-ug.vcf
        date "+%G-%m-%d %H:%M:%S"
        echo "====7. compress UG result===="
        mv output-ug.vcf ${sample_name}.vcf
        mv deduped.bam ${sample_name}.markdup.bam
        mv bqsr_report.pdf ${sample_name}.bqsr_report.pdf
        mv metrics_report.pdf ${sample_name}.metrics_report.pdf
        date "+%G-%m-%d %H:%M:%S"
        echo "====8. stat bam result===="
        myIntervals=${intervals}
        perl /bioapp/Bam_StatScripts/markdup_bam_stat_4thread.pl -ad realigned.bam -r $myIntervals -o Bam_Stat -plot
        tar zcvf Bam_Stat.tar.gz Bam_Stat/ --exclude=Bam_Stat/*.target.depth
        date "+%G-%m-%d %H:%M:%S"
    >>>
    runtime {
        docker: "public/sentieon:201808"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
    output{
        File bam_stat="Bam_Stat.tar.gz"
        File bqsr_report="${sample_name}.bqsr_report.pdf"
        File hc_vcf="${sample_name}.vcf"
        File markdup_bam="${sample_name}.markdup.bam"
        File metrics_report="${sample_name}.metrics_report.pdf"
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
        docker: "public/sentieon:3.0"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
    output{
        File fastq_stat="${sample_name}.Fastq_Statistics_Out.tar.gz"
        File bam_stats="${sample_name}.Bam_Statistics_Out.tar.gz"
    }
}

workflow WES_Germline_BWA_Sentieon_ContainRef_stat{
    File sample_fq1
    File sample_fq2
    File adaptor="public:/reference/adaptor_Trimmomatic/TruSeq3-PE-2.fa"
    File intervals
    String sample_name="sample"
    String refname="b37"
    call GD_toolkit_mapping_ContainRef as FqFltMap{
        input:
            read1=sample_fq1,
            read2=sample_fq2,
            adaptor=adaptor,
            sample_name=sample_name,
            refname=refname
    }
    call Sentieon_rmdup_realign_bqsr_hc_ContainRef_stat as RmdupVcf{
        input:
            sort_bam=FqFltMap.bam,
            intervals=intervals,
            sample_name=sample_name,
            refname=refname
    }
    call merge_fq_bam_stat as FqBamQC{
        input:
            fq_stat=FqFltMap.fqstat,
            bam_stat=RmdupVcf.bam_stat,
            sample_name=sample_name
    }
    output{
        File MarkdupBam=RmdupVcf.markdup_bam
        File HcVcf=RmdupVcf.hc_vcf
        File BqsrReport=RmdupVcf.bqsr_report
        File MetricsReport=RmdupVcf.metrics_report
        File Fqstat=FqBamQC.fastq_stat
        File Bamstats=FqBamQC.bam_stats
    }
}