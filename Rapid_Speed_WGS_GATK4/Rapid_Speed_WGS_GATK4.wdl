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
			mv chrMT.bam ${sample_id}.${read_name}.chrM.bam
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

task GATK4_HaplotypeCaller_speed{
	Array[File] bams1
	Array[File] bams2
	Array[File] bams3
	String refname
	String bam_name1
	String bam_name2
	String bam_name3
	
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

		gatk=/bioapp/gatk-4.1.4.1/gatk
		javaOption="-Xmx4g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit"
		db_BR=`echo $BR|sed "s#\&# --known-sites \$dir/#g" |xargs -i echo  "--known-sites  $dir/{}"`
		inbam1=`echo ${sep=' ' bams1}|sed 's/ / -I /g'`
		inbam2=`echo ${sep=' ' bams2}|sed 's/ / -I /g'`
		inbam3=`echo ${sep=' ' bams3}|sed 's/ / -I /g'`
		$gatk --java-options "$javaOption" MarkDuplicates -I $inbam1 -O inbam1.bam -M inbam1.mat --CREATE_INDEX true && \
		L1=`samtools view inbam1.bam -h |head -3000|grep -v ^@|head -1|cut -f3` && \
		$gatk --java-options "$javaOption" BaseRecalibrator -R $ref -I inbam1.bam $db_BR -O inbam1_recal_data.table && \
		$gatk --java-options "$javaOption" ApplyBQSR -R $ref -bqsr inbam1_recal_data.table -I inbam1.bam -L $L1 -O inbam1_bqsr.bam && \
		rm inbam1.bam && \
		$gatk --java-options "$javaOption" HaplotypeCaller -R $ref -I inbam1_bqsr.bam -L $L1 -ERC GVCF -O chr1.g.vcf.gz  && \
		$gatk --java-options "$javaOption" GenotypeGVCFs -R $ref  -V chr1.g.vcf.gz -O ${bam_name1}.vcf.gz -stand-call-conf 30 &

		$gatk --java-options "$javaOption" MarkDuplicates -I $inbam2 -O inbam2.bam -M inbam2.mat --CREATE_INDEX true && \
		L2=`samtools view inbam2.bam -h |head -3000|grep -v ^@|head -1|cut -f3` && \
		$gatk --java-options "$javaOption" BaseRecalibrator -R $ref -I inbam2.bam $db_BR -O inbam2_recal_data.table && \
		$gatk --java-options "$javaOption" ApplyBQSR -R $ref -bqsr inbam2_recal_data.table -I inbam2.bam -L $L2 -O inbam2_bqsr.bam && \
		rm inbam2.bam &&\
		$gatk --java-options "$javaOption" HaplotypeCaller -R $ref -I inbam2_bqsr.bam -L $L2  -ERC GVCF -O chr2.g.vcf.gz  && \
		$gatk --java-options "$javaOption" GenotypeGVCFs -R $ref  -V chr2.g.vcf.gz -O ${bam_name2}.vcf.gz -stand-call-conf 30 &

		$gatk --java-options "$javaOption" MarkDuplicates -I $inbam3 -O inbam3.bam -M inbam3.mat --CREATE_INDEX true && \
		L3=`samtools view inbam3.bam -h |head -3000|grep -v ^@|head -1|cut -f3` && \
		$gatk --java-options "$javaOption" BaseRecalibrator -R $ref -I inbam3.bam $db_BR  -O inbam3_recal_data.table && \
		$gatk --java-options "$javaOption" ApplyBQSR -R $ref -bqsr inbam3_recal_data.table -I inbam3.bam -L $L3 -O inbam3_bqsr.bam && \
		rm inbam3.bam &&\
		$gatk --java-options "$javaOption" HaplotypeCaller -R $ref -I inbam3_bqsr.bam -L $L3  -ERC GVCF -O chr3.g.vcf.gz  && \
		$gatk --java-options "$javaOption" GenotypeGVCFs -R $ref  -V chr3.g.vcf.gz -O ${bam_name3}.vcf.gz -stand-call-conf 30 &
    
		wait
	>>>
	output{
		File vcf1="${bam_name1}.vcf.gz"
		File vcf2="${bam_name2}.vcf.gz"
		File vcf3="${bam_name3}.vcf.gz"
	}
	runtime {
		docker: "seqflow/genedock_wgs:1.0"
    	memory: "16G"
    	disk: "100G"
    	cpu: 8
    }
}

task GATK4_HaplotypeCaller_speed_chrM{
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

		gatk=/bioapp/gatk-4.1.4.1/gatk
		javaOption="-Xmx4g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit"
		db_BR=`echo $BR|sed "s#\&# --known-sites \$dir/#g" |xargs -i echo  "--known-sites  $dir/{}"`
		inbam=`echo ${sep=' ' bams1}|sed 's/ / -I /g'`

		$gatk --java-options "$javaOption" MarkDuplicates -I $inbam -O inbam1.bam -M inbam1.mat --CREATE_INDEX true && \
    L1=`samtools view inbam1.bam -h |head -3000|grep -v ^@|head -1|cut -f3` && \
		$gatk --java-options "$javaOption" BaseRecalibrator -R $ref -I inbam1.bam $db_BR -O inbam1_recal_data.table && \
		$gatk --java-options "$javaOption" ApplyBQSR -R $ref -bqsr inbam1_recal_data.table -I inbam1.bam -L $L1 -O inbam1_bqsr.bam && \
		rm inbam1.bam && \
		$gatk --java-options "$javaOption" HaplotypeCaller -R $ref -I inbam1_bqsr.bam -L $L1 -ERC GVCF -O chr1.g.vcf.gz  && \
		$gatk --java-options "$javaOption" GenotypeGVCFs -R $ref  -V chr1.g.vcf.gz -O ${bam_name1}.vcf.gz -stand-call-conf 30
  >>>
	output{
		File vcf="${bam_name1}.vcf.gz"
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
		gatk=/bioapp/gatk-4.1.4.1/gatk
		javaOption="-Xmx4g -Djava.io.tmpdir=/var/data/tmp -XX:ParallelGCThreads=4 -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit"
		invcf=`echo ${sep=' ' vcf}|sed 's/ / -I /g'`
		$gatk --java-options "$javaOption" GatherVcfs -R $ref -I $invcf -O ${sample_id}.vcf.gz
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
    
workflow Rapid_Speed_WGS_GATK4{
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
	call GATK4_HaplotypeCaller_speed_chrM{
		input:
			bams1=align.alignChrM,
			bam_name1="chrM",
			refname=refname,
	}
  
call GATK4_HaplotypeCaller_speed as speed1{input:bams1=align.alignChr1,bams2=align.alignChr2,bams3=align.alignChr3,bam_name1="chr1",bam_name2="chr2",bam_name3="chr3",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed2{input:bams1=align.alignChr4,bams2=align.alignChr5,bams3=align.alignChr6,bam_name1="chr4",bam_name2="chr5",bam_name3="chr6",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed3{input:bams1=align.alignChr7,bams2=align.alignChr8,bams3=align.alignChr9,bam_name1="chr7",bam_name2="chr8",bam_name3="chr9",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed4{input:bams1=align.alignChr10,bams2=align.alignChr11,bams3=align.alignChr12,bam_name1="chr10",bam_name2="chr11",bam_name3="chr12",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed5{input:bams1=align.alignChr13,bams2=align.alignChr14,bams3=align.alignChr15,bam_name1="chr13",bam_name2="chr14",bam_name3="chr15",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed6{input:bams1=align.alignChr16,bams2=align.alignChr17,bams3=align.alignChr18,bam_name1="chr16",bam_name2="chr17",bam_name3="chr18",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed7{input:bams1=align.alignChr19,bams2=align.alignChr20,bams3=align.alignChr21,bam_name1="chr19",bam_name2="chr20",bam_name3="chr21",refname=refname,}
	call GATK4_HaplotypeCaller_speed as speed8{input:bams1=align.alignChr22,bams2=align.alignChrX,bams3=align.alignChrY,bam_name1="chr22",bam_name2="chrX",bam_name3="chrY",refname=refname,}
	call combinevcf{
		input:
			vcf=[GATK4_HaplotypeCaller_speed_chrM.vcf,speed1.vcf1,speed1.vcf2,speed1.vcf3,speed2.vcf1,speed2.vcf2,speed2.vcf3,speed3.vcf1,speed3.vcf2,speed3.vcf3,speed4.vcf1,speed4.vcf2,speed4.vcf3,speed5.vcf1,speed5.vcf2,speed5.vcf3,speed6.vcf1,speed6.vcf2,speed6.vcf3,speed7.vcf1,speed7.vcf2,speed7.vcf3,speed8.vcf1,speed8.vcf2,speed8.vcf3,],
			refname=refname,
			sample_id=sample_id,
	}
	output{	
		combinevcf.combinevcf
	}
}
