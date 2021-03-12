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
        docker: "seqflow/i1mapping:1.0"
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

		java -Xmx12G -jar $PICARD MarkDuplicates I=$inbam O=dups_marked.bam METRICS_FILE=out_dups_metrics.txt REMOVE_DUPLICATES=true 
		samtools index dups_marked.bam

		echo "===== step1: base recalibration ====="
		java  -Xmx12G -jar $GATK -T BaseRecalibrator -R $reference -I dups_marked.bam $db_BR -L  $intervals -o recal_data.table

		java -Xmx12G  -jar $GATK -T PrintReads -R $reference -I dups_marked.bam -L $intervals -BQSR recal_data.table  -o  ${sample_name}.recal.bam

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
		java -Xmx12G -jar $GATK -T RealignerTargetCreator -R $reference -I ${tumor} -I ${normal} $db_IR -L $intervals -o target_intervals.list
		java -Xmx12G -jar $GATK -T IndelRealigner -R $reference -I ${tumor} -I ${normal} -targetIntervals target_intervals.list $db_IR --nWayOut .realign.bam
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

workflow WES_Somatic_BWA_GATK_Mutect2_ContainRef_stat{
    File normal_fq1="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_normal_1000000_1.fq.gz"
    File normal_fq2="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_normal_1000000_2.fq.gz"
    File tumor_fq1="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_tumor_1000000_1.fq.gz"
    File tumor_fq2="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_smallsize/synthetic_challenge_set3_tumor_1000000_2.fq.gz"
    File cosmic_data="public:/database/cosmic/b37_cosmic_v54_120711.vcf"
    File dbsnp_data="public:/reference/b37_broad/dbsnp_138.b37.vcf"
    File adaptor="public:/reference/adaptor_Trimmomatic/TruSeq3-PE-2.fa"
    File intervals="public:/reference/b37_broad/Broad.human.exome.b37.interval_list"
    String refname="b37"
    call GD_toolkit_mapping_ContainRef as FqFltMapNormal{
        input:
            read1=normal_fq1,
            read2=normal_fq2,
            adaptor=adaptor,
            sample_name="normal",
            refname=refname
    }
    call GD_toolkit_mapping_ContainRef as FqFltMapTumor{
        input:
            read1=tumor_fq1,
            read2=tumor_fq2,
            adaptor=adaptor,
            sample_name="tumor",
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
    output{
        FqBamQCNormal.fastq_stat
        FqBamQCNormal.bam_stats
	FqBamQCTumor.fastq_stat
        FqBamQCTumor.bam_stats
	rmbqstatNormal.recalBam
	rmbqstatTumor.recalBam
	MuTect2.outvcf
    }
}
