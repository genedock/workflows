task split_fq{
  Array[File] fq
  Int? compression=4
  Int? thread=4
  String fq_name
  command{
    set -ex
    date "+%G-%m-%d %H:%M:%S"
    /bioapp/fastp -i ${sep=' -I ' fq} --split_by_lines 50000000 -A -G -Q -L --thread ${thread} --compression ${compression} -o ${fq_name}.1.new.gz -O ${fq_name}.2.new.gz
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
	File otherBed
	String sample_id
	String refname
	Int? thread
	String read_name
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
		$bwa mem -t ${default=8 thread} -M -R "@RG\tID:$sample_id\tPL:illumina\tPU:bar\tLB:$sample_id\tSM:$sample_id" $reference ${sep=' ' reads} |$samtools view -bST $reference - > $sample_id.bam
		java -Xmx8g -jar $picard  SortSam I=$sample_id.bam O=$sample_id.$read_name.sort.bam SO=coordinate  TMP_DIR=/var/data/tmp VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		if [ ${refname}  == "b37" ]; then
			for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
			do
			{
				samtools view -@ ${default=8 thread} -b $sample_id.$read_name.sort.bam -o $sample_id.$read_name.chr$i.bam $i
			}
			done
			mv $sample_id.$read_name.chrMT.bam $sample_id.$read_name.chrM.bam
		else
			for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
			do 
			{ 
				samtools view -@ ${default=8 thread} -b $sample_id.$read_name.sort.bam -o $sample_id.$read_name.$i.bam $i
			}
			done
		fi
		samtools view -@ ${default=8 thread} -b $sample_id.$read_name.sort.bam -L ${otherBed} -o $sample_id.$read_name.other.bam
	>>>
	output{
		File alignChr1="${sample_id}.${read_name}.chr1.bam" File alignChr2="${sample_id}.${read_name}.chr2.bam" File alignChr3="${sample_id}.${read_name}.chr3.bam" File alignChr4="${sample_id}.${read_name}.chr4.bam" File alignChr5="${sample_id}.${read_name}.chr5.bam" File alignChr6="${sample_id}.${read_name}.chr6.bam"
		File alignChr7="${sample_id}.${read_name}.chr7.bam" File alignChr8="${sample_id}.${read_name}.chr8.bam" File alignChr9="${sample_id}.${read_name}.chr9.bam" File alignChr10="${sample_id}.${read_name}.chr10.bam" File alignChr11="${sample_id}.${read_name}.chr11.bam" File alignChr12="${sample_id}.${read_name}.chr12.bam"
		File alignChr13="${sample_id}.${read_name}.chr13.bam" File alignChr14="${sample_id}.${read_name}.chr14.bam" File alignChr15="${sample_id}.${read_name}.chr15.bam" File alignChr16="${sample_id}.${read_name}.chr16.bam" File alignChr17="${sample_id}.${read_name}.chr17.bam" File alignChr18="${sample_id}.${read_name}.chr18.bam"
		File alignChr19="${sample_id}.${read_name}.chr19.bam" File alignChr20="${sample_id}.${read_name}.chr20.bam" File alignChr21="${sample_id}.${read_name}.chr21.bam" File alignChr22="${sample_id}.${read_name}.chr22.bam" File alignChrX="${sample_id}.${read_name}.chrX.bam" File alignChrY="${sample_id}.${read_name}.chrY.bam"
		File alignOther="${sample_id}.${read_name}.other.bam"
		File alignChrM="${sample_id}.${read_name}.chrM.bam"
		File alignbam="${sample_id}.${read_name}.sort.bam"
	}
	runtime {
		docker: "seqflow/genedock_wgs:1.0"
		memory: "16G"
	 	disk: "100G"
		cpu: 8
  	}
}

task mkdup_bam{
	Array[File] bams
	String sample_id
	command {
		set -x
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
  		docker: "seqflow/genedock_wgs:1.0"
  		memory: "16G"
  		disk: "1000G"
  		cpu: 4
  }
}

task GATK3_HaplotypeCaller{
	Array[File] bams1
	String refname
	String bam_name1
	
command <<<
		set -xe
		cd /var/data
		mkdir /var/data/tmp
		if [ ${refname} == "hg19" ] ;then
			ref=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
			dir=/rdata/genedock/hg19_broad
			BR="1000G_phase1.indels.hg19.sites.vcf&Mills_and_1000G_gold_standard.indels.hg19.sites.vcf&dbsnp_138.hg19.vcf"
		fi
		if [ ${refname} == "hg38" ] ;then
			ref=/rdata/genedock/hg38_broad/hg38.fasta
			dir=/rdata/genedock/hg38_broad
			BR="Mills_and_1000G_gold_standard.indels.hg38.vcf&dbsnp_138.hg38.vcf"
		fi
		if [ ${refname}  == "b37" ]; then
			ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
			dir=/rdata/genedock/b37_broad
			BR="1000G_phase1.indels.b37.vcf&Mills_and_1000G_gold_standard.indels.b37.vcf&dbsnp_138.b37.vcf"
		fi

		gatk=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
    picard=/bioapp/picard-2.23.1/picard.jar
    
		db_BR=`echo $BR|sed "s#\&# -knownSites \$dir/#g" |xargs -i echo  "-knownSites  $dir/{}"`
		inbam=`echo ${sep=' ' bams1}|sed 's/ / I=/g'`
    java -Xmx8g -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $picard MarkDuplicates I=$inbam O=inbam.bam METRICS_FILE=samplename.bam.mat CREATE_INDEX=true TMP_DIR=/var/data/tmp
    L1=`samtools view inbam.bam -h |head -3000|grep -v ^@|head -1|cut -f3` 
    java -Xmx8g -jar $gatk -T BaseRecalibrator -R $ref -I inbam.bam $db_BR -o out_recal_data.table
    java -Xmx8g -jar $gatk -T HaplotypeCaller -R $ref -I inbam.bam -BQSR out_recal_data.table -L $L1 --variant_index_type LINEAR --variant_index_parameter 128000 --emitRefConfidence GVCF -o chr.g.vcf.gz
    java -Xmx8g -jar $gatk -T GenotypeGVCFs -R $ref --variant chr.g.vcf.gz -o ${bam_name1}.vcf.gz -stand_call_conf 10
	>>>
	output{
		File vcf1="${bam_name1}.vcf.gz"
	}
	runtime {
		docker: "seqflow/genedock_wgs:1.0"
    	memory: "16G"
    	disk: "100G"
    	cpu: 8
    }
}

task combinevcf{
	Array[File] vcf
	String refname
	String sample_id
	command <<<
		set -xe
		cd /var/data
		mkdir /var/data/tmp
		if [ ${refname} == "hg19" ] ;then
			ref=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
		fi
		if [ ${refname} == "hg38" ] ;then
			ref=/rdata/genedock/hg38_broad/hg38.fasta
		fi
		if [ ${refname}  == "b37" ]; then
			ref="/rdata/genedock/b37_broad/human_g1k_v37.fasta"
		fi
		gatk=/bioapp/gatk-3.8/GenomeAnalysisTK_3.8.jar
    java -cp $gatk org.broadinstitute.gatk.tools.CatVariants -R $ref -V ${sep=' -V ' vcf} -out ${sample_id}.vcf.gz --assumeSorted
	>>>
	output{
		File combinevcf="${sample_id}.vcf.gz"
	}
	runtime {
		docker: "seqflow/genedock_wgs:1.0"
    	memory: "16G"
    	disk: "100G"
    	cpu: 8
    }
}

workflow WGS_Germline_BWA_GATK3{
	String sample_id="iw160"
	String refname="hg19"
	Array[File] lane_fq1=["genedockdx:/home/admin/wdl_test/rawdata/read1.fq.gz"]
	Array[File] lane_fq2=["genedockdx:/home/admin/wdl_test/rawdata/read2.fq.gz"]
	File otherBed="genedockdx:/home/admin/Database/hg19.other.bed"
  Int? compression=4
  Int? thread=4

	scatter(fq in transpose([lane_fq1,lane_fq2])){
		call split_fq{
			input:
				fq=fq,
        compression=compression,
        thread=thread,
        fq_name=basename(fq[0],".gz"),
		}
	}

	scatter(Fq in transpose([flatten(split_fq.read1_fq),flatten(split_fq.read2_fq)])){
		call align{
			input:
				sample_id=sample_id,
				thread=8,
				reads=Fq,
				otherBed=otherBed,
				read_name=basename(Fq[0],".gz"),
				refname=refname,
		}
	}
	call mkdup_bam{
		input:
			bams=align.alignbam,
			sample_id=sample_id,
	}
	scatter(bam in [align.alignChr1,align.alignChr2,align.alignChr3,align.alignChr4,align.alignChr5,align.alignChr6,align.alignChr7,align.alignChr8,align.alignChr9,align.alignChr10,align.alignChr11,align.alignChr12,align.alignChr13,align.alignChr14,align.alignChr15,align.alignChr16,align.alignChr17,align.alignChr18,align.alignChr19,align.alignChr20,align.alignChr21,align.alignChr22,align.alignChrX,align.alignChrY,align.alignChrM]){
		call GATK3_HaplotypeCaller{
			input:
				bams1=bam,
				bam_name1=basename(bam[0],".bam"),
				refname=refname
		}
	}
	call combinevcf{
		input:
			vcf=GATK3_HaplotypeCaller.vcf1,
			sample_id=sample_id,
			refname=refname
	}
	output{	
		combinevcf.combinevcf
		mkdup_bam.MarkdupBam
		mkdup_bam.MarkdupBai
	}
}
