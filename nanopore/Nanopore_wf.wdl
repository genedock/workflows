task split_fq{
	File fq
	command{
    set -ex
    date "+%G-%m-%d %H:%M:%S"
    fq_name=$(echo ${fq}|cut -d . -f1)
    echo "fq name is $fq_name"
    /bioapp/fastp -i ${fq} --split_by_lines 50000000 -A -G -Q -L --thread 4 --compression 4 -o $fq_name.new.gz
    ls
    rm ${fq}
    date "+%G-%m-%d %H:%M:%S"
	}
 	runtime {
		docker: "seqflow/genedock_wgs:1.1"
     	memory: "8G"
		disk: "800G" 
		cpu: 4
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
		docker: "seqflow/genedock_nanopore:1.1"
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
		docker: "seqflow/genedock_nanopore:1.1"
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
		docker: "seqflow/genedock_nanopore:1.1"
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
	command <<<
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
    awk 'BEGIN{OFS="\t"}{x=split($8,sv,";");for(i=1;i<=x;i++){if(sv[i]~"SVTYPE"){split(sv[i],type,"=")}}{if(length($4)>10000){$4="<"type[2]">"}else if(length($5)>10000){$5="<"type[2]">"}else{pass};print $0"\t"length($4)"\t"length($5)}}' ${vcf} > sv.vcf
		$ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile  sv.vcf -SVinputInfo 1 -outputFile ./SV.annotated.tsv -svtBEDcol 4
		date "+%G-%m-%d %H:%M:%S" 
  >>>
	runtime {
		docker: "seqflow/cnv_filt:1.0"
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
		medaka_variant -i ${bam} -f $ref -t 16 -b 30
		cp medaka_variant/round_1_unfiltered.vcf out.vcf
	}
	runtime {
		docker: "seqflow/genedock_nanopore:1.1"
		memory: "32G"
		disk: "800G" 
		cpu: 16
  	}
  	output{
  		File outvcf="out.vcf"
  	}
}
workflow wf_nanopore{
	File fq
	Int thread=8
  Int threshold=5
	String sm="bar"
	String id="foo"
	call split_fq{
		input:
			fq=fq
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
#	call Variant{
#		input:
#			bam=mergeBam.outbam
#	}
	output{
		mergeBam.outbam
		mergeBam.outbamstat
		SV.outSV
		SV_Anno.outSVanno
#		Variant.outvcf
	}
}
