task PrinseqLite{
  File read
  File params
  command <<<
    set -xe
    cd /var/data
    gunzip ${read} -c > read.fq
    perl /app/prinseq-lite.pl -out_format 3 -params ${params} -fastq read.fq -out_good prinseq
  >>>
  output{
    File prinseq="prinseq.fastq"
  }
  runtime {
    docker: "seqflow/prinseq:latest"
    memory: "64G"
    disk: "1000G"
    cpu: 16
  }
}

task BMtagger{
  File BitmaskFile
  File InputFastq
  File Reference_fa
  File srprismFiles
  command <<<
    set -xe
    cd /var/data
    mkdir /var/data/srprism
    mkdir /var/data/tmp
    tar -xzvf ${srprismFiles} -C /var/data/srprism
    /app/bmtagger/bmtagger.sh -d ${Reference_fa} -b ${BitmaskFile} -q1 -1 ${InputFastq} -x /var/data/srprism/Desktop/hs37.srprism/hs37.srprism  -o bm -T tmp/ -X
  >>>
  output{
    File bm="bm.fastq"
  }
  runtime {
    docker: "seqflow/bmtagger:latest"
    memory: "64G"
    disk: "1000G"
    cpu: 16
  }
}

task Humann2{
  File InputFastq1
  File InputFastq2
  File NucleotideDatabase
  File ProteinDatabase
  command <<<
    mkdir /var/data/database/
    mkdir /var/data/humann_DIR/
    mkdir /var/data/database/Ndatabase/
    tar -xzvf ${NucleotideDatabase} -C /var/data/database/Ndatabase/
    tar -xzvf ${ProteinDatabase} -C /var/data/database/
    cat ${InputFastq1} {InputFastq2} > /var/data/combined.fastq
    humann2 --nucleotide-database /var/data/database/Ndatabase/ --protein-database /var/data/database/ --input  /var/data/combined.fastq --output /var/data/humann_DIR/ --threads 8  --metaphlan-options='--min_cu_len 100'
    cp /var/data/humann_DIR/combined_genefamilies.tsv combined_genefamilies.tsv
    cp /var/data/humann_DIR/combined_humann2_temp/combined_metaphlan_bugs_list.tsv combined_metaphlan_bugs_list.tsv
    cp  /var/data/humann_DIR/combined_pathabundance.tsv combined_pathabundance.tsv
    cp  /var/data/humann_DIR/combined_pathcoverage.tsv combined_pathcoverage.tsv
>>>
  output{
    File gene="combined_genefamilies.tsv"
    File metaphlan="combined_metaphlan_bugs_list.tsv"
    File pathabundance="combined_pathabundance.tsv"
    File pathcoverage="combined_pathcoverage.tsv"
  }
  runtime {
    docker: "seqflow/humann2:latest"
    memory: "64G"
    disk: "1000G"
    cpu: 16
  }
}

workflow metagenome{
  File read1="public:/demo-data/Metagenome_Prinseq-BMtagger-HUMAnN2/SRR061695_1M_1.fastq.gz"
  File read2="public:/demo-data/Metagenome_Prinseq-BMtagger-HUMAnN2/SRR061695_1M_2.fastq.gz"
  File params="public:/demo-data/prinseq/prinseq.params"
  File BitmaskFile="public:/database/BMtagger/hs37.bitmask"
  File Reference_fa="public:/database/BMtagger/hs37.fa"
  File srprismFiles="public:/database/BMtagger/hs37.srprism.tar.gz"
  File NucleotideDatabase="public:/database/humann2/full_chocophlan.v0.1.1.tar.gz"
  File ProteinDatabase="public:/database/humann2/uniref50_annotated.tar.gz"
  call PrinseqLite as prinseq1{
    input:
    read=read1,
    params=params
  }
  call PrinseqLite as prinseq2{
    input:
      read=read2,
      params=params
  }
  call BMtagger as bm1{
    input:
      InputFastq=prinseq1.prinseq,
      BitmaskFile=BitmaskFile,
      Reference_fa=Reference_fa,
      srprismFiles=srprismFiles,
  }
  call BMtagger as bm2{
    input:
      InputFastq=prinseq2.prinseq,
      BitmaskFile=BitmaskFile,
      Reference_fa=Reference_fa,
      srprismFiles=srprismFiles,
  }
  call Humann2{
    input:
      InputFastq1=bm1.bm,
      InputFastq2=bm2.bm,
      NucleotideDatabase=NucleotideDatabase,
      ProteinDatabase=ProteinDatabase
  }
  output{
    Humann2.gene
    Humann2.metaphlan
    Humann2.pathabundance
    Humann2.pathcoverage
  }
}
