task split_fq{
	File fq
	Int? split_line
	command{
		set -x
		fq_name=$(basename ${fq} .gz)
		echo "fq name is $fq_name"
		gunzip -dc ${fq}|/bioapp/split -d -l ${default=50000000 split_line} --filter='gzip > $FILE.gz' - $fq_name.
		rm ${fq}
	}
 	runtime {
		docker: "public/genedock_wgs:1.0"
     	memory: "4096m"
		disk: "15G" 
		cpu: 2
  	}
	output{
		Array[File] output_fq = glob("*.gz")
	}
}
workflow saw{
	Array[File] proband_fq1
	Array[File] proband_fq2
	Int split_line=50000000
	scatter(fq1 in proband_fq1){
		call split_fq as split_proband_fq1{
			input:
				fq=fq1,
				split_line=split_line,
		}
	}
	scatter(fq2 in proband_fq2){
		call split_fq as split_proband_fq2{
			input:
				fq=fq2,
				split_line=split_line,
		}
	}
}
