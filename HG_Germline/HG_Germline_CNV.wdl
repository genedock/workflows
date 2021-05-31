#HG_WDL
task split_fq{
  Array[File] fq
  String fq_name
  command{
    set -ex
    date "+%G-%m-%d %H:%M:%S"
    /bioapp/fastp -i ${sep=' -I ' fq} --split_by_lines 50000000 -A -G -Q -L --thread 4 --compression 4 -o ${fq_name}.1.new.gz -O ${fq_name}.2.new.gz
    date "+%G-%m-%d %H:%M:%S"
  }
  runtime{
    docker: "genedockdx/genedock_wgs:1.1"
    memory: "8192m"
    disk: "409600m"
    cpu: 4
  }
  output{
    	Array[File] read1_fq = glob("*.1.new.gz")
    	Array[File] read2_fq = glob("*.2.new.gz")
	}
}

task align{
	Array[File] reads
	String sample_id
	String refname
	Int? thread
	String read_name
  File adaptor
  Int base_quality_type=33
    Int crop=10000
    Int head_crop=0
    Int trailing=3
    Int leading=3
    Int minlength=36
    String sliding_window="4:15"
    String illumina_clip="2:30:10"
	command <<<
		set -x
		if [ ${refname} == "hg19" ] ;then
			reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
			otherbed=/bioapp/otherbed/hg19.other.bed
		fi
		if [ ${refname} == "hg38" ] ;then
			reference=/rdata/genedock/hg38_broad/hg38.fasta
			otherbed=/bioapp/otherbed/hg38.other.bed
		fi
		if [ ${refname}  == "b37" ]; then
			reference="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
			otherbed=/bioapp/otherbed/b37.other.bed
		fi

		bwa=/bioapp/bwa-0.7.17/bwa
		samtools=/bioapp/samtools-1.10/samtools
		picard=/bioapp/picard-2.23.1/picard.jar
		sample_id=${sample_id}
		read_name=${read_name}
		mkdir -p /var/data/tmp
    java -Xmx16g -jar /bioapp/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 ${sep=' ' reads} forward_paired.fq.gz forward_unpaired.fq.gz reverse_paired.fq.gz reverse_unpaired.fq.gz ILLUMINACLIP:${adaptor}:${illumina_clip} SLIDINGWINDOW:${sliding_window} LEADING:${leading} TRAILING:${trailing} CROP:${crop} HEADCROP:${head_crop} MINLEN:${minlength}
		$bwa mem -t ${default=8 thread} -M -R "@RG\tID:$sample_id\tPL:illumina\tPU:bar\tLB:$sample_id\tSM:$sample_id" $reference forward_paired.fq.gz reverse_paired.fq.gz |$samtools view -bST $reference - > $sample_id.bam
		java -Xmx8g -jar $picard  SortSam I=$sample_id.bam O=$sample_id.$read_name.sort.bam SO=coordinate  TMP_DIR=/var/data/tmp VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

if [ ~{refname} = "hg19" ] || [ ~{refname} = "hg38" ] ;then
array="chr1|chr2\nchr3|chr4|chr5\nchr6|chr7|chr8\nchr9|chr10|chr11|chr12\nchr13|chr14|chr15|chr16\nchr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY|chrM"
elif [ ~{refname} = "b37" ] ; then
array="1|2\n3|4|5\n6|7|8\n9|10|11|12\n13|14|15|16\n17|18|19|20|21|22|X|Y|MT"
fi

index=0
echo  -e $array|while read chrom;do index=$(( $index + 1 )) 
	chr=`samtools view -H $sample_id.$read_name.sort.bam|grep "@SQ" | awk 'NR<=25{split($2,x,":");print x[2]}'|grep -E -w "$chrom"`
    if [ $? -eq 0 ];then
      echo $chr|awk -v BAM=$sample_id.$read_name.sort.bam -v INDEX=$index '{cmd="samtools view -bh "BAM" "$0" > "INDEX".bam";system(cmd)}'
    else
      samtools view -H $sample_id.$read_name.sort.bam > $index.bam
    fi
done

    cp 1.bam ${sample_id}.${read_name}.1.bam
    cp 2.bam ${sample_id}.${read_name}.2.bam
    cp 3.bam ${sample_id}.${read_name}.3.bam
    cp 4.bam ${sample_id}.${read_name}.4.bam
    cp 5.bam ${sample_id}.${read_name}.5.bam
    cp 6.bam ${sample_id}.${read_name}.6.bam
    mkdir -p /var/data/Fastq_Data
    basename1=${sample_id}"_sample_1.fq.gz"
    basename2=${sample_id}"_sample_2.fq.gz"
    mv `echo ${sep=' ' reads}|awk '{print $1}'` /var/data/Fastq_Data/$basename1
    mv `echo ${sep=' ' reads}|awk '{print $2}'` /var/data/Fastq_Data/$basename2
    /bioapp/Fastq_StatScripts/gd_qc_stat -i /var/data/Fastq_Data/$basename1,/var/data/Fastq_Data/$basename2 -a ${adaptor},${adaptor} -N 0.1 -q ${base_quality_type} -L 15 -p 0.5 -c -k -o /var/data/Fastq_Stat
    tar zcvf Fastq_Stat.tar.gz Fastq_Stat/ --exclude=Fastq_Stat/*.fq --exclude=Fastq_Stat/*.fastq
    cp Fastq_Stat.tar.gz ${sample_id}.${read_name}.Fastq_Stat.tar.gz
    >>>
	output{
		File bam0="${sample_id}.${read_name}.1.bam"
		File bam1="${sample_id}.${read_name}.2.bam"
		File bam2="${sample_id}.${read_name}.3.bam"
		File bam3="${sample_id}.${read_name}.4.bam"
		File bam4="${sample_id}.${read_name}.5.bam"
		File bam5="${sample_id}.${read_name}.6.bam"
		File bam="${sample_id}.${read_name}.sort.bam"
    	File fqstat="${sample_id}.${read_name}.Fastq_Stat.tar.gz"
	}
	runtime {
		docker: "public/genedock_wgs:1.0"
		memory: "16G"
	 	disk: "100G"
		cpu: 8
  	}
}

task GD{
  String ref="b37"
  command <<<
    set -ex
    if [ ${ref} == "hg19" ]; then
        interval="/rdata/genedock/pipeline_config/WGS_script/ref/hg19/hg19_noN.bed"
    elif [ ${ref} == "b37" ]; then
        interval="/rdata/genedock/pipeline_config/WGS_script/ref/b37/b37_noN.bed"
    fi
    if [ ~{ref} = "hg19" ] || [ ~{ref} = "hg38" ] ;then
      array="chr1|chr2\nchr3|chr4|chr5\nchr6|chr7|chr8\nchr9|chr10|chr11|chr12\nchr13|chr14|chr15|chr16\nchr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY|chrM"
    elif [ ~{ref} = "b37" ] ; then
      array="1|2\n3|4|5\n6|7|8\n9|10|11|12\n13|14|15|16\n17|18|19|20|21|22|X|Y|MT"
    else
      echo "the type parameter can only be b37 or hg19 or hg38 "
      exit 1
    fi
    index=0
    echo -e $array|while read chrom;do index=$(( $index + 1 ))
    awk -v chrom=$chrom -v FS="\t" '$1!~/^@/{if(match($1,"^("chrom")$")){print $0}}' $interval  >  $index.intervals
    ls
    done
  >>>
  output{
  File intervals0="1.intervals"
  File intervals1="2.intervals"
  File intervals2="3.intervals"
  File intervals3="4.intervals"
  File intervals4="5.intervals"
  File intervals5="6.intervals"
  }
    runtime {
        docker: "public/genedock_wgs:1.0"
        memory: "1G"
        disk: "10G"
        cpu: 1
    }
}

task DNAseq_HC{
    Array[File] bams
    Int nt
    String refname
    String sample_name
    File bed
    String chrlist
    command <<<
        export SENTIEON_LICENSE=node3.genedock.com:16527
        date "+%G-%m-%d %H:%M:%S"
        set -ex
        cd /var/data
        if [ ${refname} == "hg19" ]; then
            ref="/rdata/genedock/hg19_broad/ucsc.hg19.fasta"
            dir="/rdata/genedock/hg19_broad"
            known_dbsnp="$dir/dbsnp_138.hg19.vcf"
            known_1000GIndels="$dir/1000G_phase1.indels.hg19.sites.vcf"
            known_Mills="$dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
            
        fi
        if [ ${refname} == "b37" ]; then
            ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
            dir="/rdata/genedock/b37_broad"
            known_dbsnp="$dir/dbsnp_138.b37.vcf"
            known_1000GIndels="$dir/1000G_phase1.indels.b37.vcf"
            known_Mills="$dir/Mills_and_1000G_gold_standard.indels.b37.vcf"
        fi
        inbam=`echo ${sep=' ' bams}|sed 's/ / I=/g'`
        picard=/bioapp/picard-2.23.8/picard.jar
        java -Xmx8g -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $picard MarkDuplicates I=$inbam O=inbam.bam METRICS_FILE=samplename.bam.mat CREATE_INDEX=true TMP_DIR=/var/data/tmp
        echo "===== step1: Remove Duplicate Reads ===="
        dedup_param="--traverse_param 1000000/10000"
        sentieon driver -r $ref -t ${nt} -i inbam.bam $dedup_param --algo LocusCollector --fun score_info ${sample_name}.score.txt.gz \
            --algo MeanQualityByCycle ${sample_name}.mq_metrics.txt --algo QualDistribution ${sample_name}.qd_metrics.txt \
            --algo GCBias --summary ${sample_name}.gc_summary.txt ${sample_name}.gc_metrics.txt --algo AlignmentStat ${sample_name}.aln_metrics.txt --algo InsertSizeMetricAlgo ${sample_name}.is_metrics.txt
        sentieon driver -r $ref -t ${nt} -i inbam.bam \
                     --traverse_param 1000000/10000 --algo Dedup --rmdup \
                     --score_info ${sample_name}.score.txt.gz \
                     --metrics ${sample_name}.dedup_metrics.txt \
                     --bam_compression 1 ${sample_name}.deduped.bam
        date "+%G-%m-%d %H:%M:%S"
        
        echo "===== step2: Indel Realignment ===="
        
        sentieon driver -t ${nt} -r $ref -i ${sample_name}.deduped.bam --algo Realigner  -k $known_Mills -k $known_1000GIndels ${sample_name}.realigned.bam
        date "+%G-%m-%d %H:%M:%S"
        
        echo "===== step3: BQSR ===="

        sentieon driver -r $ref  -t ${nt} -i ${sample_name}.realigned.bam  --algo QualCal -k $known_dbsnp -k $known_1000GIndels -k $known_Mills ${sample_name}.recal_data.table
        date "+%G-%m-%d %H:%M:%S"
        
        echo "===== step4: HC Variant caller===="
        hc_option="--emit_conf=30 --call_conf=30"
        sentieon driver -r $ref -t ${nt} -i ${sample_name}.realigned.bam  -q ${sample_name}.recal_data.table --algo Haplotyper  --emit_mode gvcf ${sample_name}.${chrlist}.output-hc.g.vcf.gz && 
        sentieon driver -r $ref -t ${nt} --algo GVCFtyper $hc_option -v ${sample_name}.${chrlist}.output-hc.g.vcf.gz -d $known_dbsnp ${sample_name}.${chrlist}.output-hc.vcf.gz
        
        perl /bioapp/Bam_StatScripts/markdup_bam_stat_4thread.pl -ad inbam.bam  -r ${bed} -o Bam_Stat -plot
        tar zcvf Bam_Stat.tar.gz Bam_Stat/ --exclude=Bam_Stat/*.target.depth
        cp Bam_Stat.tar.gz ${sample_name}.${chrlist}.Bam_Stat.tar.gz
        mv inbam.bam ${sample_name}.${chrlist}.bam
    >>>
    output{
        File gvcf="${sample_name}.${chrlist}.output-hc.g.vcf.gz"
        File vcf="${sample_name}.${chrlist}.output-hc.vcf.gz"
        File bam_stat="${sample_name}.${chrlist}.Bam_Stat.tar.gz"
        File bam = "${sample_name}.${chrlist}.bam"
    }
    runtime {
        docker: "public/sentieon:201911"
        memory: "32G"
        disk: "400G"
        cpu: 8
    }
}
task combinevcf{
	Array[File] vcf
	String sample_id
	command <<<
		set -xe
		cd /var/data
		mkdir /var/data/tmp
    picard=/bioapp/picard-2.23.1/picard.jar
    java -Xmx6G -jar $picard MergeVcfs I=${sep=' I=' vcf} O=${sample_id}.vcf.gz 
  >>>
	output{
		File combinevcf="${sample_id}.vcf.gz"
	}
	runtime {
		docker: "seqflow/genedock_wgs:1.0"
    	memory: "8G"
    	disk: "100G"
    	cpu: 4
    }
}
task merge_fq_bam_stat{
    Array[File] fq_stat
    Array[File] bam_stat
    command <<<
        cd /var/data
        set -x
        echo "====1. merge fq stat===="
        for i in `echo ${sep=" " fq_stat}`;do
        dir=$(basename $i .tar.gz)
        mkdir $dir
        tar zxvf $i -C $dir
        done
        
        python /bioapp/Fastq_StatScripts/merge_stat.py /var/data/*Fastq_Stat/*/*.stat >fastq_statistics.txt
        python /bioapp/Fastq_StatScripts/merge_GC.py  /var/data/*Fastq_Stat/*/raw_*.GC >raw_sample.GC
        python /bioapp/Fastq_StatScripts/merge_GC.py  /var/data/*Fastq_Stat/*/clean_*.GC >clean_sample.GC
        python /bioapp/Fastq_StatScripts/merge_QD.py  /var/data/*Fastq_Stat/*/raw_*.QD >raw_sample.QD
        python /bioapp/Fastq_StatScripts/merge_QD.py  /var/data/*Fastq_Stat/*/clean_*.QD >clean_sample.QD
        python /bioapp/Fastq_StatScripts/merge_QM.py -raw fastq_statistics.txt /var/data/*Fastq_Stat/*/raw_*.QM >raw_sample.QM
        python /bioapp/Fastq_StatScripts/merge_QM.py  -clean fastq_statistics.txt /var/data/*Fastq_Stat/*/clean_*.QM >clean_sample.QM
        Rscript /bioapp/Fastq_StatScripts/plot_raw-clean_GC-QD-QM.R raw_sample.GC raw_sample.QD raw_sample.QM clean_sample.GC clean_sample.QD clean_sample.QM
        Rscript /bioapp/Fastq_StatScripts/plot_raw-clean_GC-QD-QM_pdf.R raw_sample.GC raw_sample.QD raw_sample.QM clean_sample.GC clean_sample.QD clean_sample.QM
        mkdir Fastq_Statistics_Out
        cp fastq_statistics.txt *.GC *.QD *.QM *.png *.pdf Fastq_Statistics_Out/
        tar zcvf Fastq_Statistics_Out.tar.gz Fastq_Statistics_Out/
        echo "====2. merge bam stat===="
        for i in `echo ${sep=" " bam_stat}`;do
        dir=$(basename $i .tar.gz)
        mkdir $dir
        tar zxvf $i -C $dir
        done
        mkdir /var/data/Bam_Statistics_Out
        cd /var/data/Bam_Statistics_Out
        python /bioapp/Bam_StatScripts/merge_bam_stat.py /var/data/*Bam_Stat/*/*.information.xls >mardup_bam_stat.xls
        python /bioapp/Bam_StatScripts/merge_cumu.py /var/data/*Bam_Stat/*/*.cumu.xls >cumu.xls
        python /bioapp/Bam_StatScripts/merge_depth_frequency.py /var/data/*Bam_Stat/*/*.depth_frequency.xls >depth_frequency.xls
        Rscript /bioapp/Bam_StatScripts/cumuPlot.R cumu.xls
        Rscript /bioapp/Bam_StatScripts/histPlot.R  depth_frequency.xls
        cd /var/data
        tar zcvf Bam_Statistics_Out.tar.gz Bam_Statistics_Out/
    >>>
    runtime {
        docker: "genedockdx/genedock_wgs:1.1"
        memory: "8G"
        disk: "400G"
        cpu: 4
    }
    output{
        File fastq_stat="Fastq_Statistics_Out.tar.gz"
        File bam_stats="Bam_Statistics_Out.tar.gz"
    }
}
#HG WDL
task Sentieon_vqsr_ContainRef{
  input{
  File raw_vcf
  String sample_id
  String refname="hg19"
  Int thread=8
  }
  command <<<
    set -ex
    export SENTIEON_LICENSE=node3.genedock.com:16527
    date "+%G-%m-%d %H:%M:%S"
    cd /var/data
    if [ ${refname} == "hg19" ]; then
      fasta="/rdata/genedock/hg19_broad/ucsc.hg19.fasta"
      dir="/rdata/genedock/hg19_broad"
      dbsnp_Mill=$dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
      dbsnp_1000G_omni=$dir/1000G_omni2.5.hg19.sites.vcf
      dbsnp_hapmap=$dir/hapmap_3.3.hg19.sites.vcf
      dbsnp_1000G_phase1=$dir/1000G_phase1.snps.high_confidence.hg19.sites.vcf
      dbsnp_1000G_phase1_indel=$dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
      dbsnp=$dir/dbsnp_138.hg19.vcf
    elif [ ${refname} == "b37" ]; then
      fasta="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
      dir="/rdata/genedock/b37_broad"
      dbsnp_Mill=$dir/Mills_and_1000G_gold_standard.indels.b37.vcf
      dbsnp_1000G_omni=$dir/1000G_omni2.5.b37.vcf
      dbsnp_hapmap=$dir/hapmap_3.3.b37.vcf
      dbsnp_1000G_phase1=$dir/1000G_phase1.snps.high_confidence.b37.vcf
      dbsnp_1000G_phase1_indel=$dir/1000G_phase1.indels.b37.vcf
      dbsnp=$dir/dbsnp_138.b37.vcf
    fi
    raw_vcf=${raw_vcf}

    sentieon util vcfindex $raw_vcf

    echo "===1.snp==="
    echo "create the resource argument"

    resource_text="--resource $dbsnp_1000G_phase1 --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
    resource_text="$resource_text --resource $dbsnp_1000G_omni --resource_param omni,known=false,training=true,truth=true,prior=12.0 "
    resource_text="$resource_text --resource $dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
    resource_text="$resource_text --resource $dbsnp_hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"

    echo "create the annotation argument"

    annotation_array="DP QD MQ MQRankSum ReadPosRankSum FS"
    for annotation in $annotation_array; do
       annotate_text="$annotate_text --annotation $annotation"
    done

    echo "Run the VQSR"
    sentieon driver -r $fasta --algo VarCal -v $raw_vcf $resource_text $annotate_text --var_type SNP  --tranche 90 --tranche 99 --tranche 99.9 --tranche 100 --plot_file vqsr_SNP.hc.plot_file.txt --max_gaussian 4 --srand 47382911 --tranches_file vqsr_SNP.hc.tranches vqsr_SNP.hc.recal

    echo "apply the VQSR"
    sentieon driver -r $fasta --algo ApplyVarCal -v $raw_vcf --var_type SNP --recal vqsr_SNP.hc.recal --tranches_file vqsr_SNP.hc.tranches --sensitivity 99 ${sample_id}_vqsrSNP.hc.recal.vcf

    cat vqsr_SNP.hc.tranches 
    date "+%G-%m-%d %H:%M:%S"

    echo "===2.indel==="
    echo "create the resource argument"
    resource_text="--resource $dbsnp_Mill --resource_param Mills,known=false,training=true,truth=true,prior=12.0 "
    resource_text="$resource_text --resource $dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "

    echo "create the annotation argument"
    annotation_array="DP QD MQ ReadPosRankSum FS"
    annotate_text=""
    for annotation in $annotation_array; do
      annotate_text="$annotate_text --annotation $annotation"
    done

    echo "Run the VQSR"
    sentieon driver -r $fasta --algo VarCal -v ${sample_id}_vqsrSNP.hc.recal.vcf $resource_text $annotate_text --var_type INDEL --plot_file vqsr_SNP_INDEL.hc.plot_file.txt --max_gaussian 4 --srand 47382911  --tranche 90 --tranche 99 --tranche 99.9 --tranche 100 --tranches_file vqsr_SNP_INDEL.hc.tranches vqsr_SNP_INDEL.hc.recal

    echo "apply the VQSR"
    sentieon driver -r $fasta --algo ApplyVarCal -v ${sample_id}_vqsrSNP.hc.recal.vcf --var_type INDEL --recal vqsr_SNP_INDEL.hc.recal --tranches_file vqsr_SNP_INDEL.hc.tranches --sensitivity 99 ${sample_id}_vqsr.hc.recal.vcf

    cat  vqsr_SNP_INDEL.hc.tranches
    date "+%G-%m-%d %H:%M:%S"

>>>
  runtime {
    docker: "public/sentieon:201911"
    memory: "16G"
    disk: "400G" 
    cpu: 8
    }
  output{
    File snp_filter_vcf="${sample_id}_vqsrSNP.hc.recal.vcf"
    File snp_indel_filter_vcf="${sample_id}_vqsr.hc.recal.vcf"
}
}
task cnvnator{
  String chrom
  String refname
  File bam
  String bin_size
  String sample
  command <<<
    if [ ${refname} == "hg19" ] ;then
      ref=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
      ref_chr="chrM chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    fi
    if [ ${refname}  == "b37" ]; then
      ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
      ref_chr="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
    fi
    if [ ${refname}  == "hs37d5" ]; then
      ref="/rdata/genedock/hs37d5_broad/hs37d5.fa"
      ref_chr="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
    fi

    echo Job of ${chrom} started!
    if [[ ${chrom} == "ALL" ]];
    then
        rootname=${chrom}
            chrom=$ref_chr
        for chr in $chrom;do
                /opt/CNVnator_v0.3.2/src/samtools/samtools faidx $ref $chr >$chr.fa
        done
    else
        for chr in ${chrom};do
                /opt/CNVnator_v0.3.2/src/samtools/samtools faidx $ref $chr >$chr.fa
        done
        rootname=$(echo ${chrom}|sed -e 's/ /-/g')

    fi
    cnvnator -root $rootname.root -chrom $chrom -genome $ref -tree ${bam} -unique 2>err
    cnvnator -root $rootname.root -chrom $chrom -his ${bin_size} -eval ${bin_size} -d /var/data/ 2>>err
    cnvnator -root $rootname.root -chrom $chrom -stat ${bin_size} -eval ${bin_size}  2>>err
    cnvnator -root $rootname.root -chrom $chrom -partition ${bin_size} -eval ${bin_size}  2>>err
    cnvnator -root $rootname.root -chrom $chrom -call ${bin_size} -eval ${bin_size}  >tmp 2>>err
    grep -E "duplication|deletion" tmp >/var/data/${sample}_CNVnator.txt
    mv *.root  /var/data/${sample}_CNVnator.root
    echo Job of ${chrom} finished!
    ls
  >>>
  runtime {
    docker: "public/cnvnator:1.0"
    memory: "32768m"
    disk: "400G"
    cpu: 8
  }
  output {
    File out_root="/var/data/${sample}_CNVnator.root"
    File out_txt="/var/data/${sample}_CNVnator.txt"
  }
}
task manta{
  String refname
  String sample_name
  File bam
  File bai
  String enableRemoteReadRetrievalForInsertionsInGermlineCallingModes
  command <<<
  set -ex
    if [ ${refname} == "hg19" ] ;then
      ref=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
    fi
    if [ ${refname}  == "b37" ]; then
      ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
    fi
    if [ ${refname}  == "hs37d5" ]; then
      ref="/rdata/genedock/hs37d5_broad/hs37d5.fa"
    fi
    if [ ${enableRemoteReadRetrievalForInsertionsInGermlineCallingModes} == "0" ];then
      config=/bioapp/build/bin/configManta-0.py.ini
    fi
    if [ ${enableRemoteReadRetrievalForInsertionsInGermlineCallingModes} == "1" ];then
      config=/bioapp/build/bin/configManta.py.ini
    fi
    cd /var/data
    /bioapp/build/bin/configManta.py --bam ${bam} --referenceFasta $ref --runDir analysis --config=$config
    cd analysis
    ./runWorkflow.py -m local -j 15
    python /bioapp/build/libexec/convertInversion.py /usr/local/bin/samtools $ref  results/variants/diploidSV.vcf.gz > /var/data/${sample_name}_manta.vcf
  >>>
  runtime{
    docker: "genedockdx/sv_tools:1.0"
    memory: "65536m"
    disk: "400G"
    cpu: 16
  }
  output {
    File out_vcf="/var/data/${sample_name}_manta.vcf"
  }
}
task erds{
  File bam
  File bai
  File gvcf
  String refname
  String sample_name
  command <<<
    set -ex
    if [ ${refname} == "hg19" ] ;then
      ref=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
    fi
    if [ ${refname}  == "b37" ]; then
      ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
    fi
    if [ ${refname}  == "hs37d5" ]; then
      ref="/rdata/genedock/hs37d5_broad/hs37d5.fa"
    fi
    cd /var/data
    date
    perl /bioapp/ERDS/erds_tcag/src/erds_pipeline.pl -b ${bam} -v ${gvcf} -o /var/data/out -r $ref
    mv /var/data/out/*erds.vcf /var/data/${sample_name}_ERDS.vcf
  >>>
  runtime{
    docker: "xhmdp/erds:2.0"
    memory: "32768m"
    disk: "400G"
    cpu: 8
  }
  output{
    File ERDS_vcf="/var/data/${sample_name}_ERDS.vcf"
  }
}
task mkdup_bam{
  Array[File] bams
  String sample_id
  command {
    set -ex
    picard=/bioapp/picard-2.23.1/picard.jar
    inbam=`echo ${sep=' ' bams}|sed 's/ / I=/g'`
    echo $inbam
    mkdir /var/data/tmp
    java -Xmx8g -jar $picard MarkDuplicates I=$inbam O=${sample_id}_markdup.bam M=out.mkdup.metrics AS=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=/var/data/tmp
  }
  output{
    File MarkdupBam="${sample_id}_markdup.bam"
    File MarkdupBai="${sample_id}_markdup.bai"
  }
  runtime {
      docker: "public/genedock_wgs:1.0"
      memory: "16G"
      disk: "1000G"
      cpu: 4
  }
}
task merge_cnv {
  File cnvnator_txt
  File ERDS_vcf
  String sample_name
  command <<<
    set -ex
    date "+%G-%m-%d %H:%M:%S"
    cd /var/data/
    sed -i 's#>=5#>=10#g'  /bioapp/TCGA/filt_annotSV.pl
    export PATH=/root/miniconda3/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    PATH=/root/miniconda3/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    #ERDS
    mkdir erds
    cd erds
    /root/miniconda3/bin/python /bioapp/TCGA/convert_CNV_calls_to_common_format.py ${ERDS_vcf} ERDS > ERDS.bed
    #annotation
    ln -s /rdata/genedock/AnnotSV/Annotations_Human /bioapp/AnnotSV_2.2/share/doc/AnnotSV
    export ANNOTSV=/bioapp/AnnotSV_2.2
    export PATH=$ANNOTSV/bin:$PATH
    ANNOTSV="/bioapp/AnnotSV_2.2"
    $ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile  ERDS.bed -SVinputInfo 1 -outputFile ./sample.annotated.tsv -svtBEDcol 4
    ##RLCR annotation
    perl /bioapp/TCGA/filt_annotSV.pl  sample.annotated.tsv
    #RLCR #out:filt.uniq.bed.RLCR
    cut -f 1-8 filt.bed |uniq > filt.uniq.bed
    /root/miniconda3/bin/python /bioapp/TCGA/compare_with_RLCR_definition.py /bioapp/TCGA/RLCRs.txt filt.uniq.bed
    perl /bioapp/script/map_population.pl
    perl /bioapp/script/filt_final.pl all_final.xls |cut -f 1-78> /var/data/erds_filt_final.tsv
    cd ../
    
    mkdir cnvnator
    cd cnvnator
    #CNVnator
    /root/miniconda3/bin/python /bioapp/TCGA/convert_CNV_calls_to_common_format.py ${cnvnator_txt} CNVnator > CNVnator.bed
    #annotation
    $ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile CNVnator.bed -SVinputInfo 1 -outputFile ./sample.annotated.tsv -svtBEDcol 4
    ##RLCR annotation
    perl /bioapp/TCGA/filt_annotSV.pl  sample.annotated.tsv
    #RLCR #out:filt.uniq.bed.RLCR
    cut -f 1-8 filt.bed |uniq > filt.uniq.bed
    /root/miniconda3/bin/python /bioapp/TCGA/compare_with_RLCR_definition.py /bioapp/TCGA/RLCRs.txt filt.uniq.bed
    perl /bioapp/script/map_population.pl
    perl /bioapp/script/filt_final.pl all_final.xls |cut -f 1-78> /var/data/cnvnator_filt_final.tsv
    cd ../

    date "+%G-%m-%d %H:%M:%S"
  >>>
  runtime {
    docker: "genedockdx/sv_filt:1.0"
    memory: "16G"
    disk: "100G"
    cpu: 4
  }
  output {
    File erds_filter="/var/data/erds_filt_final.tsv"
    File cnvnator_filter="/var/data/cnvnator_filt_final.tsv"
      }
}
task sv_filt{
  File manta_vcf
  String sample_name
  command <<<
    set -ex
    cd /var/data/
    sed -i 's#>=5#>=10#g'  /bioapp/TCGA/filt_annotSV.pl
    ##filter_population and filter_PR10
    zless ${manta_vcf} > sample_manta.vcf
    grep '#' sample_manta.vcf > sample.vcf
    #cp sample_manta.vcf sample.vcf
    grep -v '#' sample_manta.vcf|awk '{split($9,a,":");split($10,b,":");x=split($9,a,":");for(i=1;i<=x;i++){if(a[i]=="PR"){split(b[i],pr,",")}else if(a[i]=="SR"){split(b[i],sr,",")}else{pass}};if(sr[2]>=10||pr[2]>=10){print $0};sr[2]=0;pr[2]=0}' >> sample.vcf
    export PATH=/root/miniconda3/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    PATH=/root/miniconda3/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    mv /bioapp/script/filt_annotSV.pl   /bioapp/TCGA/filt_annotSV.pl
    mv /bioapp/script/map_population.pl  /bioapp/TCGA/map_population.pl
    #annotation
    ln -s /rdata/genedock/AnnotSV/Annotations_Human /bioapp/AnnotSV_2.2/share/doc/AnnotSV
    export ANNOTSV=/bioapp/AnnotSV_2.2
    export PATH=$ANNOTSV/bin:$PATH
    ANNOTSV="/bioapp/AnnotSV_2.2"
    $ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile  sample.vcf  -SVinputInfo 1 -outputFile ./sample.annotated.tsv -svtBEDcol 4
    ##RLCR annotation
    perl /bioapp/TCGA/filt_annotSV.pl  sample.annotated.tsv
    #RLCR #out:filt.uniq.bed.RLCR
    cut -f 1-8 filt.bed |uniq > filt.uniq.bed
    /root/miniconda3/bin/python /bioapp/TCGA/compare_with_RLCR_definition.py /bioapp/TCGA/RLCRs.txt filt.uniq.bed
    perl /bioapp/TCGA/map_population.pl
    perl /bioapp/script/filt_final.pl all_final.xls |cut -f 1-78> filt_final.tsv
    # add vcfanno end info
    ln -s  /bioapp/script/ref.bed.gz /var/data/refGene.exons.hg19.bed.gz
    ln -s  /bioapp/script/ref.bed.gz.tbi /var/data/refGene.exons.hg19.bed.gz.tbi
    ls /var/data/
    sed -i 's/5/4/g' /bioapp/TCGA/vcfanno_example/conf.toml
    sed -i 's/concat/uniq/g' /bioapp/TCGA/vcfanno_example/conf.toml
    cat /bioapp/TCGA/vcfanno_example/conf.toml
    /bioapp/TCGA/vcfanno_linux64 -lua /bioapp/TCGA/vcfanno_example/custom.lua -ends /bioapp/TCGA/vcfanno_example/conf.toml sample.vcf   > sample.anno_end.vcf
    python /bioapp/script/add_end.py  filt_final.tsv sample.anno_end.vcf sample.annotSV.end.tsv
    ##anno BND
    cut -f 1-81 sample.annotSV.end.tsv|awk  -F '\t' '{if(NR==1){print $0"\tgene_anno\tBND_link"}else{if($79=="."){gene_anno=$5"_intergenic"}else if($79~"exon"){gene_anno=$5"_exon"}else if($79~"gene"&&$79!~"exon"){gene_anno=$5"_intron"};if($5=="BND"){if($6~"]"){split($6,bnd,"]");a="-->"}else{split($6,bnd,"[");a="<--"}if(length(bnd[1])>=1){b=bnd[1]""a}else{b=a""bnd[3]}}else{b="."}{print $0"\t"gene_anno"\t"b}}}' > sample_anno.tsv
    grep -v split  sample_anno.tsv > full.tsv
    grep -w BND full.tsv|cut -f  1-3,6,82|awk -F '\t' '{if($4~"]"){split($4,a,"]")}else{split($4,a,"[")}{print $0"\tchr"$1":"$2"\t"a[2]}}' > bnd.txt
    Rscript /bioapp/script/filter_intergenic.r
    sed 's/^chr//g' intergenic.txt |awk -F ':' '{print $1"\t"$2}'|grep -v -w -F -f - full.tsv|grep -v -E 'INS_intergenic|DEL_intergenic|DUP_intergenic|INV_intergenic' > final_full.tsv
    awk -F '\t' 'BEGIN{OFS="\t"}{$79=substr($79,1,32000);$80=substr($80,1,32000);$81=substr($81,1,32000);$78=substr($78,1,32000);print $0}'  final_full.tsv > ${sample_name}_final_full.xls
    awk -F '\t' 'BEGIN{OFS="\t"}{$79=substr($79,1,32000);$80=substr($80,1,32000);$81=substr($81,1,32000);$78=substr($78,1,32000);print $0}'  sample_anno.tsv > ${sample_name}_anno.xls 
    date
  >>>
  runtime{
    docker: "genedockdx/sv_filt:1.0"
    memory: "16384m"
    disk: "400G"
    cpu: 2
  }
  output {
    File sv_anno_filt_all="/var/data/${sample_name}_anno.xls"
    File sv_anno_filt_full="/var/data/${sample_name}_final_full.xls"
  }
}
workflow Sentieon_WGS_DNAseq{
	Int nt=8
	String sample_id
	String refname="hg19"
	Array[File] lane_fq1
	Array[File] lane_fq2
	File adaptor
	String enableRemoteReadRetrievalForInsertionsInGermlineCallingModes="1"
	String chrom="ALL"
	Int bin_size=100
	scatter(fq in transpose([lane_fq1,lane_fq2])){
		call split_fq{
			input:
				fq=fq,
        fq_name=basename(fq[0],".gz"),
		}
	}

	scatter(Fq in transpose([flatten(split_fq.read1_fq),flatten(split_fq.read2_fq)])){
		call align{
			input:
				sample_id=sample_id,
				thread=8,
				reads=Fq,
				read_name=basename(Fq[0],".gz"),
				refname=refname,
        adaptor=adaptor,
		}
	}
    call GD{
      input:
        ref=refname,
    }
    call DNAseq_HC as hc0{
      input:
         bams=align.bam0,
         sample_name=sample_id,
         refname=refname,
         nt=nt,
         bed=GD.intervals0,
         chrlist="chrlist0"
    }
    call DNAseq_HC as hc1{
      input:
         bams=align.bam1,
         sample_name=sample_id,
         refname=refname,
         nt=nt,
         bed=GD.intervals1,
         chrlist="chrlist1"
    }
     call DNAseq_HC as hc2{
      input:
         bams=align.bam2,
         sample_name=sample_id,
         refname=refname,
         nt=nt,
         bed=GD.intervals2,
         chrlist="chrlist2"
    }
    call DNAseq_HC as hc3{
      input:
         bams=align.bam3,
         sample_name=sample_id,
         refname=refname,
         nt=nt,
         bed=GD.intervals3,
         chrlist="chrlist3"
    }
    call DNAseq_HC as hc4{
      input:
         bams=align.bam4,
         sample_name=sample_id,
         refname=refname,
         nt=nt,
         bed=GD.intervals4,
         chrlist="chrlist4"
    }
    call DNAseq_HC as hc5{
      input:
         bams=align.bam5,
         sample_name=sample_id,
         refname=refname,
         nt=nt,
         bed=GD.intervals5,
         chrlist="chrlist5"
    }
    call merge_fq_bam_stat{
      input:
        fq_stat=align.fqstat,
        bam_stat=[hc0.bam_stat,hc1.bam_stat,hc2.bam_stat,hc3.bam_stat,hc4.bam_stat,hc5.bam_stat,]
    }
    call combinevcf{
      input:
      vcf=[hc0.vcf,hc1.vcf,hc2.vcf,hc3.vcf,hc4.vcf,hc5.vcf],
      sample_id=sample_id,
    }
    call Sentieon_vqsr_ContainRef{
      input:
      	raw_vcf=combinevcf.combinevcf,
      	sample_id=sample_id,
      	refname=refname,
      	thread=nt,
    }
	### merge bam
	call mkdup_bam{
	  input:
		bams=align.bam,
		sample_id=sample_id,
	}
	##### cnvnator
	call cnvnator{
      input:
		chrom=chrom,
		refname=refname,
		bam=mkdup_bam.MarkdupBam,
		bin_size=bin_size,
		sample=sample_id,
	}
	##### manta
	call manta{
	  input:
		refname=refname,
		sample_name=sample_id,
		bam=mkdup_bam.MarkdupBam,
		bai=mkdup_bam.MarkdupBai,
		enableRemoteReadRetrievalForInsertionsInGermlineCallingModes=enableRemoteReadRetrievalForInsertionsInGermlineCallingModes,
	}
	##### erds
	call erds{
	  input:
		bam=mkdup_bam.MarkdupBam,
		bai=mkdup_bam.MarkdupBai,
		gvcf=combinevcf.combinevcf,
		refname=refname,
		sample_name=sample_id,
	}
	##### merge_cnv
	call merge_cnv {
	  input:
		cnvnator_txt=cnvnator.out_txt,
		ERDS_vcf=erds.ERDS_vcf,
		sample_name=sample_id,
	}
	##### sv_anno_filt
	call sv_filt{
	  input:
		manta_vcf=manta.out_vcf,
		sample_name=sample_id,
  }
    output{
        File vcf=combinevcf.combinevcf
        File fqstat=merge_fq_bam_stat.fastq_stat
        File bamstat=merge_fq_bam_stat.bam_stats
        File bam=mkdup_bam.MarkdupBam
        File bai=mkdup_bam.MarkdupBai
        File vqsr_vcf=Sentieon_vqsr_ContainRef.snp_indel_filter_vcf
        File manta_out=manta.out_vcf
        File cnvnator_out=cnvnator.out_txt
        File erds_out=erds.ERDS_vcf
        File cnv_filter_erds=merge_cnv.erds_filter
        File cnv_filter_cnvnator=merge_cnv.cnvnator_filter
        File sv_filter_all=sv_filt.sv_anno_filt_all
        File sv_filter_full=sv_filt.sv_anno_filt_full
    }
}
