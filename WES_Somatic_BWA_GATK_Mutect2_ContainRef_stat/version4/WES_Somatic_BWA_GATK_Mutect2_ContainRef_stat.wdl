#WES_Somatic_BWA_GATK_Mutect2_ContainRef_stat
task soapnuke{
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
        File filter_stat
        File soapnuke_stat
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
                perl ${soapnuke_stat} Basic_Statistics_of_Sequencing_Quality.txt Statistics_of_Filtered_Reads.txt > result/${sample_id}.filter.stat.xls
                perl ${filter_stat} -indir  `pwd`   -output ${sample_id}_output.xls
                mkdir soapnuke_stat
                mv *.txt soapnuke_stat/
                tar -zcvf ${sample_id}_soapnuke.tar.gz soapnuke_stat/
        }
        output{
                File clean_read1="/var/data/${sample_id}_SOAPnuke/${sample_id}_clean.1.fq.gz"
                File clean_read2="/var/data/${sample_id}_SOAPnuke/${sample_id}_clean.2.fq.gz"
                File stat="/var/data/${sample_id}_SOAPnuke/result/${sample_id}.filter.stat.xls"
                File out="/var/data/${sample_id}_SOAPnuke/${sample_id}_output.xls"
                File soapnuke_stats="/var/data/${sample_id}_SOAPnuke/${sample_id}_soapnuke.tar.gz"
        }
        runtime {
                docker: "genedockdx/qc:1.4"
                memory: "8G"
                disk: "400G"
                cpu: 4
        }
}
task BWA{
    input{
    File clean_read1
    File clean_read2
    String sample_id
    String read_id
    String read_sm
    String read_lib="lib"
    String refname
    Int thread
    }
    command <<<
        set -ex
        if [ ${refname} == "hg19" ] ;then
            reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
        elif [ ${refname} == "b37" ] ;then
            reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
        elif [ ${refname} == "hg38" ] ;then
            reference=/rdata/genedock/hg38_broad/hg38.fasta
        elif [ ${refname} == "hs37d5" ] ;then
            reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
        fi        
        echo "===== step1: Bwa alignment and samtools sort ===="
        bwa mem -t ${default=8 thread} -M -R "@RG\tID:${read_id}\tSM:${read_sm}\tPL:illumina\tLB:${read_lib}" $reference ${clean_read1} ${clean_read2} > pe.sam
        samtools sort -@ 8 -o ${sample_id}_sorted.bam -O bam pe.sam 
        samtools index  ${sample_id}_sorted.bam 
    >>>
    output{
        File bam="${sample_id}_sorted.bam"
        File bai="${sample_id}_sorted.bam.bai"
    }
    runtime {
        docker: "public/i1mapping:1.0"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
}
task rmdup_bqsr_normal{
  input{
    File inbam
    File inbai
    File intervals
    String refname
    String sample_id
  }
  command <<<
    set -ex
    if [ ${refname} == "hg19" ] ;then
      reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
      dir="/rdata/genedock/hg19_broad"
      BR="1000G_phase1.indels.hg19.sites.vcf&Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&dbsnp_138.hg19.vcf"
    elif [ ${refname} == "hg38" ] ;then
      reference=/rdata/genedock/hg38_broad/hg38.fasta
      dir="/rdata/genedock/hg38_broad"
      BR="Mills_and_1000G_gold_standard.indels.hg38.vcf&dbsnp_138.hg38.vcf"
    elif [ ${refname} == "b37" ];then
      reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
      dir="/rdata/genedock/b37_broad"
      BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
    elif [ ${refname} == "hs37d5" ]; then
      reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
      dir="/rdata/genedock/b37_broad"
      BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
    fi

    PICARD=/bioapp/picard-2.23.1/picard.jar
    GATK=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
    intervals=${intervals}

    db_BR=`echo "$BR"|sed "s#\&# -knownSites \$dir/#g" |xargs -i echo  "-knownSites $dir/{}"`
    echo "===== step0: picard MarkDuplicates ====="

    java -Xmx8G -jar $PICARD MarkDuplicates I=${inbam} O=dups_marked.bam METRICS_FILE=out_dups_metrics.txt REMOVE_DUPLICATES=false 
    samtools index dups_marked.bam

    echo "===== step1: base recalibration ====="
    java  -Xmx8G -jar $GATK -T BaseRecalibrator -R $reference -I dups_marked.bam $db_BR -L  $intervals -o recal_data.table

    java -Xmx8G  -jar $GATK -T PrintReads -R $reference -I dups_marked.bam -L $intervals -BQSR recal_data.table  -o  ${sample_id}.normal.recal.bam
        >>>
        output{
          File recalBam="${sample_id}.normal.recal.bam"
          File recalBai="${sample_id}.normal.recal.bai"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}
task rmdup_bqsr_tumor{
  input{
    File inbam
    File inbai
    File intervals
    String refname
    String sample_id
  }
  command <<<
    set -ex
    if [ ${refname} == "hg19" ] ;then
      reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
      dir="/rdata/genedock/hg19_broad"
      BR="1000G_phase1.indels.hg19.sites.vcf&Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&dbsnp_138.hg19.vcf"
    elif [ ${refname} == "hg38" ] ;then
      reference=/rdata/genedock/hg38_broad/hg38.fasta
      dir="/rdata/genedock/hg38_broad"
      BR="Mills_and_1000G_gold_standard.indels.hg38.vcf&dbsnp_138.hg38.vcf"
    elif [ ${refname} == "b37" ];then
      reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
      dir="/rdata/genedock/b37_broad"
      BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
    elif [ ${refname} == "hs37d5" ]; then
      reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
      dir="/rdata/genedock/b37_broad"
      BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
    fi

    PICARD=/bioapp/picard-2.23.1/picard.jar
    GATK=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
    intervals=${intervals}
    db_BR=`echo "$BR"|sed "s#\&# -knownSites \$dir/#g" |xargs -i echo  "-knownSites $dir/{}"`
    echo "===== step0: picard MarkDuplicates ====="

    java -Xmx8G -jar $PICARD MarkDuplicates I=${inbam} O=dups_marked.bam METRICS_FILE=out_dups_metrics.txt REMOVE_DUPLICATES=true 
    samtools index dups_marked.bam

    echo "===== step1: base recalibration ====="
    java  -Xmx8G -jar $GATK -T BaseRecalibrator -R $reference -I dups_marked.bam $db_BR -L  $intervals -o recal_data.table

    java -Xmx8G  -jar $GATK -T PrintReads -R $reference -I dups_marked.bam -L $intervals -BQSR recal_data.table  -o  ${sample_id}.tumor.recal.bam
        >>>
        output{
          File recalBam="${sample_id}.tumor.recal.bam"
          File recalBai="${sample_id}.tumor.recal.bai"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}
task realign{
  input{
    File normal
    File normal_bai
    File tumor
    File tumor_bai
    File intervals
    String refname
    String sample_id
  }
  command <<<
    set -ex
    if [ ${refname} == "hg19" ] ;then
      reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
      dir="/rdata/genedock/hg19_broad"
      knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&1000G_phase1.indels.hg19.sites.vcf"
    elif [ ${refname} == "hg38" ] ;then
      reference=/rdata/genedock/hg38_broad/hg38.fasta
      dir="/rdata/genedock/hg38_broad"
      knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.hg38.vcf&1000G_phase1.indels.hg38.vcf"
    elif [ ${refname} == "b37" ];then
      reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
      dir="/rdata/genedock/b37_broad"
      knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.b37.vcf&1000G_phase1.indels.b37.vcf"
    elif [ ${refname} == "hs37d5" ]; then
      reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
      dir="/rdata/genedock/b37_broad"
      knowsites_IndelRealigner="Mills_and_1000G_gold_standard.indels.b37.vcf&1000G_phase1.indels.b37.vcf"
    fi

    db_IR=`echo "$knowsites_IndelRealigner"|sed "s#\&# -known \$dir/#g" |xargs -i echo  "-known $dir/{}"` 
    GATK=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
    intervals=${intervals}
    java -Xmx8G -jar $GATK -T RealignerTargetCreator -R $reference -I ${normal} -I ${tumor} $db_IR -L $intervals -o target_intervals.list
    java -Xmx8G -jar $GATK -T IndelRealigner -R $reference -I ${normal} -I ${tumor}  -targetIntervals target_intervals.list $db_IR --nWayOut .realign.bam
    mv *.tumor.recal.realign.bam ${sample_id}_tumor.recal.realign.bam
    mv *.normal.recal.realign.bam ${sample_id}_normal.recal.realign.bam
    ls 
        >>>
        output{
          File normalBam="${sample_id}_normal.recal.realign.bam"
          File tumorBam="${sample_id}_tumor.recal.realign.bam"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}
task bamstat{
  File bam
  File region
  String sample_id
  command <<< 
   set -ex
   samtools index -@ 8 ${bam}
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
task MuTect2{
  input{
    File normal="genedockdx:/home/admin/yuce_somatic_wes/WES/call-realign/synthetic_normal.recal.realign.bam"
    File tumor="genedockdx:/home/admin/yuce_somatic_wes/WES/call-realign/synthetic_tumor.recal.realign.bam"
    File gnomad="genedockdx:/home/admin/script/somatic-b37_af-only-gnomad.raw.sites.vcf.gz"
    String normal_sample="normal"
    File intervals="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_NGv3/NGv3.bed"
    String refname="b37"
    String sample_id
    File pon
  }
  command <<<
    set -x
    if [ ${refname} == "hg19" ] ;then
      reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
    elif [ ${refname} == "b37" ] ;then
      reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
    elif [ ${refname} == "hg38" ] ;then
      reference=/rdata/genedock/hg38_broad/hg38.fasta
    elif [ ${refname} == "hs37d5" ] ;then
      reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
    fi   
    ls *.bam|while read f;do samtools index $f ;done 
    gatk=/bioapp/gatk-4.1.4.1/gatk
    tabix -p vcf ${gnomad}
    tabix -p vcf ${pon}
    normal=`samtools view -H ${normal}|grep '^@RG'|tr "\t" "\n"|grep SM|awk -F ':' '{print $2}'`
    $gatk Mutect2 -R $reference -I ${normal} -I ${tumor} -normal $normal -L ${intervals}  --germline-resource ${gnomad} --panel-of-normals ${pon}  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter  --f1r2-tar-gz tumor.f1r2.tar.gz  --bam-output  tumor.bamout.bam  -O raw.vcf.gz
    $gatk GetPileupSummaries -R $reference  -I ${tumor} -V  ${gnomad} -L ${intervals}  --interval-set-rule INTERSECTION  -O tumor.getpileupsum.table
    $gatk GetPileupSummaries -R $reference  -I ${normal} -V ${gnomad} -L ${intervals}   --interval-set-rule INTERSECTION  -O normal.getpileupsum.table

    $gatk CalculateContamination -I tumor.getpileupsum.table   -matched normal.getpileupsum.table -O tumor.contamination.table --tumor-segmentation tumor.segmentation.table
    
    $gatk LearnReadOrientationModel -I tumor.f1r2.tar.gz -O tumor.artifact.tar.gz 
    ##filter
    $gatk FilterMutectCalls -R $reference  -V raw.vcf.gz  --contamination-table tumor.contamination.table --tumor-segmentation tumor.segmentation.table  --ob-priors tumor.artifact.tar.gz -stats raw.vcf.gz.stats -O tumor.filtered.vcf.gz
    mv raw.vcf.gz ${sample_id}_raw.vcf.gz
    mv tumor.filtered.vcf.gz ${sample_id}_filtered.vcf.gz
    
    >>>
        output{
          File rawvcf="${sample_id}_raw.vcf.gz"
          File filtvcf="${sample_id}_filtered.vcf.gz"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "400G"
                cpu: 8
        }
}
task vcf2maf{
  File vepvcf
  String normal_id
  String tumor_id
  String sample_id
  String refname
  File cnvkit_bed
  command <<< 
    set -ex
    if [ ${refname} == "hg19" ] ;then
      reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
    elif [ ${refname} == "b37" ] ;then
      reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
    elif [ ${refname} == "hg38" ] ;then
      reference=/rdata/genedock/hg38_broad/hg38.fasta
    elif [ ${refname} == "hs37d5" ] ;then
      reference=/rdata/genedock/hs37d5_broad/hs37d5.fa
    fi
    ###vcf2maf
    perl /bioapp/mskcc-vcf2maf-2235eed/vcf2maf.pl --input-vcf ${vepvcf} --output-maf ${sample_id}_filter.maf --inhibit-vep --ref-fasta $reference --normal-id ${normal_id} --tumor-id ${tumor_id}
    ### vcf_list
    echo """
Variant_Classification
Variant_Type
dbSNP_RS
BIOTYPE
AF
Consequence""" > vcf_list
cat vcf_list
    maf=${sample_id}_filter.maf
    name_id=0
    eval `grep -v '#' $maf|head -1|awk -F '\t' '{for(i=1;i<=NF;i++){if($i=="Variant_Type"){print "variant_type="i}else if($i=="BIOTYPE"){print "biotype="i}}}'`
    grep -v '#' $maf|awk 'BEGIN{col="'$variant_type'"}{if($col=="Variant_Type"||$col=="SNP"){print $0}}' > snp.maf
    grep -v '#' $maf|awk 'BEGIN{col="'$variant_type'"}{if($col=="Variant_Type"||$col!="SNP"){print $0}}' > indel.maf
    ls snp.maf indel.maf|while read file;do
    cat vcf_list|while read name;do 
      name_id=`grep -v '#' $file|head -1|awk -F '\t' '{for(i=1;i<=NF;i++){if($i=="'$name'"){print i}}}'`
      if [ $name_id == 0 ]; then echo "$name is null"
      else
        if [[ $name == "Variant_Classification" || $name == "Consequence" ]]; then
          grep -v '#' $file|grep -v $name|cut -f $biotype,$name_id|tr "," "\n"|awk '{if($2!="protein_coding"){$2="non_protein_coding";print $0}else{print $0}}'|sort|uniq -c|awk '{print $2"\t"$3"\t"$1}' >> $file.$name.txt
        name_id=0
        elif [ $name == "dbSNP_RS" ];then
        grep -v '#' $file|grep -v dbSNP_RS|cut -f $name_id|awk '{if($0!~"^rs"){sum++}}END{print "Fraction of SNPs in dbsnp (%)\t"sum/NR*100}' >>  $file.list.txt
        elif [ $name == "AF" ];then
        grep -v '#' $file|grep -v "dbSNP_RS"|cut -f $name_id|awk '{if($0==""){pass}else{sum++}}END{print "Fraction of SNPs in 1000genomes (%)\t"sum/NR*100}'  >>  $file.list.txt
        fi
        awk 'END{print "total variants\t"NR-1}' $file >> $file.list.txt
      fi
      done
    done
    gatk=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
    java -Xmx8g -jar $gatk -T SelectVariants -R $reference -V ${vepvcf}  -selectType SNP -o snp.vcf
    bcftools stats -s- snp.vcf|grep '^TSTV'|awk '{print "ti/tv\t"$5}' >> snp.maf.list.txt

    grep '#'  ${vepvcf} > cds.vcf
    grep -w -E 'coding_sequence_variant|protein_coding' ${vepvcf} >> cds.vcf
    bcftools stats -s- cds.vcf > bcftools.txt
    grep '^IDD' bcftools.txt|awk 'BEGIN{print "InDel distribution\tlength\tcount"}{print $1"\t"$3"\t"$4}' > cds_indel_LenCount.txt
    
    rm bcftools.txt
    mkdir ${sample_id}_stat
    mv *.txt  ${sample_id}_stat
    tar -zcvf ${sample_id}_stat.tar.gz  ${sample_id}_stat
    awk 'BEGIN{print "Sample\t'${sample_id}'"}{if($4=="Loss"){loss_count++;loss_size+=$3-$2}else{gain_count++;gain_size+=$3-$2}}END{print "gain_count\t"gain_count"\ngain_size\t"gain_size"\nloss_count\t"loss_count"\nloss_size\t"loss_size"\ntotal_count\t"gain_count+loss_count"\ntotal_size\t"loss_size+gain_size}' ${cnvkit_bed} > /var/data/${sample_id}_cnv_info.txt

  >>>
 runtime {
        docker: "genedockdx/qc:1.4"
        memory: "4G"
        disk: "50G" 
        cpu: 2
    }
 output{
   File snp_maf="snp.maf"
   File indel_maf="indel.maf"
   File sample_maf="${sample_id}_filter.maf"
   File sample_stat="${sample_id}_stat.tar.gz"
   File sample_cnv_stat="/var/data/${sample_id}_cnv_info.txt"
  }
}
task cnvkit{
  File tumor_bam
  File normal_bam
  String refname
  File target_bed
  File refFlat
  String sample_id
  String batch_method="hybrid"
  String call_method="clonal"
  Float purity=0.8 
  File cnvkit_gene
  String loss_threshold="-0.9"
  String gain_threshold="1"
  command <<<
  set -ex
  export PATH=/root/miniconda3/bin:/root/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
   if [ ${refname} == "hg19" ] ;then
     reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
     zless ${refFlat} > /var/data/refFlat.txt
   elif [ ${refname} == "hg38" ] ;then
     reference=/rdata/genedock/hg38_broad/hg38.fasta
     zless ${refFlat} > /var/data/refFlat.txt
   elif [ ${refname} == "b37" ];then
     reference=/rdata/genedock/b37_broad/human_g1k_v37.fasta
     zless ${refFlat}|sed 's/chr//g' > /var/data/refFlat.txt
   fi
   /root/miniconda3/bin/cnvkit.py batch  ${tumor_bam} --normal ${normal_bam} --targets ${target_bed} --annotate /var/data/refFlat.txt --fasta $reference --access ${target_bed} --output-reference my_reference.cnn  --output-dir results  -m ${batch_method}
   /root/miniconda3/bin/cnvkit.py call results/*.cns -m  ${call_method} -o ${sample_id}.call.cns
   /root/miniconda3/bin/cnvkit.py scatter results/*.cnr -s ${sample_id}.call.cns -o ${sample_id}.call-scatter.pdf
   /root/miniconda3/bin/cnvkit.py diagram results/*.cnr -s ${sample_id}.call.cns -o ${sample_id}.call-diagram.pdf
   /root/miniconda3/bin/cnvkit.py call results/*.cnr --purity ${purity} -o ${sample_id}.adjust.cnr
   python ${cnvkit_gene} ${sample_id}.adjust.cnr ${loss_threshold} ${gain_threshold}|sort -k 1,1 -k 2,2g > /var/data/${sample_id}.cnv.bed
    >>>
      runtime {
        docker: "genedockgva/cnvkit:2.0"
        memory: "16G"
        disk: "400G"
        cpu: 8
    }
    output{
      File cnvkit_bed="/var/data/${sample_id}.cnv.bed"
      File  cnvkit_cns="/var/data/${sample_id}.call.cns"
      File  scatter_out="/var/data/${sample_id}.call-scatter.pdf"
      File  diagram_out="/var/data/${sample_id}.call-diagram.pdf"
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
        memory: "8G"
        disk: "200G"
        cpu: 4
    }
}
workflow somatic{
  input{
  Array[File] tumor_fq1=["genedockdx:/home/admin/WGS_data/test_1.fq.gz","genedockdx:/home/admin/yucebio/tumor_1.fq.gz"]
  Array[File] tumor_fq2=["genedockdx:/home/admin/WGS_data/test_2.fq.gz","genedockdx:/home/admin/yucebio/tumor_2.fq.gz"]
  Array[File] normal_fq1=["genedockdx:/home/admin/yucebio/normal_1.fq.gz"]
  Array[File] normal_fq2=["genedockdx:/home/admin/yucebio/normal_2.fq.gz"]
  String sample_id
  File intervals="public:/demo-data/WES_Somatic_Synthetic_Challenge_Set3_NGv3/NGv3.bed"
  String adapter_read1="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
  String adapter_read2="CAACTCCTTGGCTCACAGAACGACATGGCTACGATCCGACTT"
  String refname="b37"
  Float misMatch=1
  Float nRate=0.1
  Float lowQual=5
  Float qualRate=0.5
  Int qualSys=2
  Int bwa_thread=8
  String loss_threshold="-0.9"
  String gain_threshold="1"
  Float purity=0.8
  File pon="genedockdx:/home/admin/script/somatic-b37_Mutect2-exome-panel.vcf.gz"
  File cnvkit_gene="genedockdx:/home/admin/script/cnvkit_gene.py"
  File filter_stat="genedockdx:/ref/script/filter_stat.pl"
  File soapnuke_stat="genedockdx:/ref/script/soapnuke_stat.pl"
  String call_method="clonal"
  String batch_method="amplicon"
  File refFlat="genedockdx:/home/admin/script/refFlat.gencode.v11.gz"
  File germline_vcf="genedockdx:/home/admin/script/somatic-b37_af-only-gnomad.raw.sites.vcf.gz"
  }
  String tumor_sample="${sample_id}_tumor"
  String normal_sample="${sample_id}_normal"
  Int array_length_tumor = length(tumor_fq1)
  Int array_length_normal = length(normal_fq1)
    call soapnuke as soapnuke_tumor{
      input:
        read1=tumor_fq1,
        read2=tumor_fq2,
        sample_id=tumor_sample,
        adapter_read1=adapter_read1,
        adapter_read2=adapter_read2,
        misMatch=misMatch,
        lowQual=lowQual,
        nRate=nRate,
        qualRate=qualRate,
        qualSys=qualSys,
        filter_stat=filter_stat,
        soapnuke_stat=soapnuke_stat,
    }
    call BWA as tumor_bwa{
      input:
        clean_read1=soapnuke_tumor.clean_read1,
        clean_read2=soapnuke_tumor.clean_read2,
        sample_id=tumor_sample,
        read_id=tumor_sample,
        read_sm=tumor_sample,
        read_lib="lib",
        refname=refname,
        thread=bwa_thread,
    }
    call rmdup_bqsr_tumor{
      input:
        inbam=tumor_bwa.bam,
        inbai=tumor_bwa.bai,
        intervals=intervals,
        refname=refname,
        sample_id=tumor_sample,
    }
    call soapnuke as soapnuke_normal{
      input:
        read1=normal_fq1,
        read2=normal_fq2,
        sample_id=normal_sample,
        adapter_read1=adapter_read1,
        adapter_read2=adapter_read2,
        misMatch=misMatch,
        lowQual=lowQual,
        nRate=nRate,
        qualRate=qualRate,
        qualSys=qualSys,
        filter_stat=filter_stat,
        soapnuke_stat=soapnuke_stat,
    }
    call BWA as normal_bwa{
      input:
        clean_read1=soapnuke_normal.clean_read1,
        clean_read2=soapnuke_normal.clean_read2,
        sample_id=normal_sample,
        read_id=normal_sample,
        read_sm=normal_sample,
        read_lib="lib",
        refname=refname,
        thread=bwa_thread,
    }
    call rmdup_bqsr_normal{
      input:
        inbam=normal_bwa.bam,
        inbai=normal_bwa.bai,
        intervals=intervals,
        refname=refname,
        sample_id=normal_sample,
    }
  call realign{
    input:
      normal=rmdup_bqsr_normal.recalBam,
      normal_bai=rmdup_bqsr_normal.recalBai,
      tumor=rmdup_bqsr_tumor.recalBam,
      tumor_bai=rmdup_bqsr_tumor.recalBai,
      intervals=intervals,
      refname=refname,
      sample_id=sample_id,
  }
  call bamstat as normal_bamstat{
    input:
      bam=realign.normalBam,
      region=intervals,
      sample_id=normal_sample,

  }
  call bamstat as tumor_bamstat{
    input:
      bam=realign.tumorBam,
      region=intervals,
      sample_id=tumor_sample,
  }
  call MuTect2{
    input:
      normal=realign.normalBam,
      tumor=realign.tumorBam,
      gnomad=germline_vcf,
      intervals=intervals,
      refname=refname,
      normal_sample="normal",
      sample_id=sample_id,
      pon=pon,
  }
  call cnvkit{
    input:
      tumor_bam=realign.tumorBam,
      normal_bam=realign.normalBam,
      refname=refname,
      target_bed=intervals,
      refFlat=refFlat,
      sample_id=sample_id,
      batch_method=batch_method,
      call_method=call_method,
      purity=purity,
      cnvkit_gene=cnvkit_gene,
      loss_threshold=loss_threshold,
      gain_threshold=gain_threshold,
  }
  call VEP{
    input:
      vcf=MuTect2.filtvcf,
      refseq="Yes",
      refname=refname,
      sample_id=sample_id,
  }
  call vcf2maf{
    input:
      vepvcf=VEP.out,
      normal_id=normal_sample,
      tumor_id=tumor_sample,
      sample_id=sample_id,
      refname=refname,
      cnvkit_bed=cnvkit.cnvkit_bed
  }
  output{
    soapnuke_tumor.clean_read1
    soapnuke_tumor.clean_read2
    soapnuke_tumor.stat
    soapnuke_tumor.out
    soapnuke_tumor.soapnuke_stats
    soapnuke_normal.clean_read1
    soapnuke_normal.clean_read2
    soapnuke_normal.stat
    soapnuke_normal.out
    soapnuke_normal.soapnuke_stats
    realign.normalBam
    realign.tumorBam
    MuTect2.rawvcf
    MuTect2.filtvcf
    cnvkit.cnvkit_cns
    cnvkit.cnvkit_bed
    cnvkit.scatter_out
    cnvkit.diagram_out
    VEP.out
    normal_bamstat.outstat
    normal_bamstat.out_txt
    tumor_bamstat.outstat
    tumor_bamstat.out_txt
    vcf2maf.sample_maf
    vcf2maf.snp_maf
    vcf2maf.indel_maf
    vcf2maf.sample_stat
    vcf2maf.sample_cnv_stat
  }
}
