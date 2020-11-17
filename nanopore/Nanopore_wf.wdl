task split_fq{
	File fq
	Int split_line
	command{
		set -x
		fq_name=$(basename ${fq} .gz)
		echo "fq name is $fq_name"
		pigz -dc ${fq}|/bioapp/split -d -l ${default=50000000 split_line} --filter='pigz > $FILE.gz' - $fq_name.
		rm ${fq}
	}
 	runtime {
		docker: "public/genedock_wgs:1.0"
     	memory: "4G"
		disk: "800G" 
		cpu: 2
  	}
	output{
		Array[File] output_fq = glob("*.gz")
	}
}

task minimap2{
  File fq
  Int thread
  String sm
  String id
  String read_name
  command{
  		reference=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
    	minimap2 -R "@RG\tID:${id}\tSM:${sm}" -t ${thread} -Y --MD -ax map-ont $reference ${fq}|samtools sort -O BAM  -o ${sm}.${read_name}.sort.bam
      #samtools sort -@ ${thread} -O BAM  -o ${sm}.sort.bam out.sam
  }
  runtime {
		docker: "multitest/genedock_nanopore:1.1"
		memory: "32G"
		disk: "1000G" 
		cpu: 8
  	}
  output{
  		File outbam="${sm}.${read_name}.sort.bam"
  }
}
task mergeBam{
	Array [File] bams
	command{
		samtools merge merge.bam ${sep=' ' bams} -c -p
		samtools index -@ 4 merge.bam
		mosdepth stat merge.bam
	}
	runtime {
		docker: "multitest/genedock_nanopore:1.1"
		memory: "16G"
		disk: "800G" 
		cpu: 4
  	}
  	output{
 		File outbam="merge.bam"
 		File outbamstat="stat.mosdepth.summary.txt"
 	}
}
task SV{
	File bam
	Int thread
	Int threshold
	command{
		sniffles -m ${bam} -v output.sv.vcf --genotype -t ${thread} -s ${threshold}
	}
	runtime {
		docker: "multitest/genedock_nanopore:1.1"
		memory: "16G"
		disk: "400G" 
		cpu: 4
  	}
  	output{
  		File outSV="output.sv.vcf"
  	}
}
task SV_Anno{
	File vcf
	command {
		date "+%G-%m-%d %H:%M:%S"
		cd /var/data/
		set -ex
		export PATH=/root/miniconda3/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
		PATH=/root/miniconda3/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
		#annotation
		ln -s /rdata/genedock/AnnotSV/Annotations_Human /bioapp/AnnotSV_2.2/share/doc/AnnotSV
		export ANNOTSV=/bioapp/AnnotSV_2.2
		export PATH=$ANNOTSV/bin:$PATH
		ANNOTSV="/bioapp/AnnotSV_2.2"
		$ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile  ${vcf} -SVinputInfo 1 -outputFile ./SV.annotated.tsv -svtBEDcol 4  
		date "+%G-%m-%d %H:%M:%S" 
	}
	runtime {
		docker: "multitest/cnv_filt:2.0"
		memory: "8G"
		disk: "100G" 
		cpu: 4
  	}
  	output{
  		File outSVanno="/var/data/SV.annotated.tsv"
  	}
}
task Variant{
	File bam
	command{
		ref=/rdata/genedock/hg19_broad/ucsc.hg19.fasta
		export PATH=/usr/local/miniconda/envs/medaka/bin/:$PATH
    samtools index ${bam}
		medaka_variant -i ${bam} -f $ref -t 8 -b 12
		cp medaka_variant/round_1_unfiltered.vcf out.vcf
	}
	runtime {
		docker: "multitest/genedock_nanopore:1.1"
		memory: "16G"
		disk: "800G" 
		cpu: 8
  	}
  	output{
  		File outvcf="out.vcf"
  	}
}
workflow wf_nanopore{
	File fq
	Int split_line=50000000
	Int thread=8
  Int threshold=5
	String sm="bar"
	String id="foo"
	call split_fq{
		input:
			fq=fq,
			split_line=split_line,
	}
	scatter(sfq in split_fq.output_fq){
		call minimap2{
			input:
				fq=sfq,
				thread=thread,
				sm=sm,
				id=id,
				read_name=basename(sfq,".gz")
		}
	}
	call mergeBam{
		input:
			bams=minimap2.outbam
	}
	call SV{
		input:
			bam=mergeBam.outbam,
      thread=thread,
      threshold=threshold
	}
	call SV_Anno{
		input:
			vcf=SV.outSV
	}
	call Variant{
		input:
			bam=mergeBam.outbam
	}
	output{
		mergeBam.outbam
		mergeBam.outbamstat
		SV.outSV
		SV_Anno.outSVanno
		Variant.outvcf
	}
}
