import "depend/WES_Germline_BWA_GATK4_single.wdl" as sample_wes

workflow multi_wes{
  input{
    Array[Array[File]]  lane_fq1=[["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_1.fastq.gz"],["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_1.fastq.gz"]]
    Array[Array[File]]  lane_fq2=[["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_2.fastq.gz"],["public:/demo-data/WES-Germline_NA12878_smallsize/NA12878-NGv3-LAB1360-A_1000000_2.fastq.gz"]]
    File bed="public:/demo-data/WES-Germline_NA12878_NGv3-LAB1360/NA12878-NGv3-LAB1360-A-regions.bed"
    Array[String] sample_lb=["lb1","lb2"]
    String sample_pu="pu"
    Array[String] sample_id=["test1","test2"]
    String sample_pl="illumina"
    Int stand_call_conf=30
    String adapter_read1="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
    String adapter_read2="CAACTCCTTGGCTCACAGAACGACATGGCTACGATCCGACTT"
    String refname="hs37d5"
    Float misMatch=1
    Float nRate=0.1
    Float lowQual=5
    Float qualRate=0.5
    Int qualSys=2
    Int bwa_thread=8
    String refseq="Yes"
}
  Int num = length(sample_id)
  scatter(index in range(num)){
  call sample_wes.wf_WES as WES_sample{
    input:
      lane_fq1=lane_fq1[index],
      lane_fq2=lane_fq2[index],
      bed=bed,
      sample_lb=sample_lb[index],
      sample_pu=sample_pu,
      sample_id=sample_id[index],
      sample_pl=sample_pl,
      refname=refname,
      stand_call_conf=stand_call_conf,
      adapter_read1=adapter_read1,
      adapter_read2=adapter_read2,
      misMatch=misMatch,
      nRate=nRate,
      lowQual=lowQual,
      qualRate=qualRate,
      qualSys=qualSys,
      bwa_thread=bwa_thread,
      refseq=refseq,
  }
  }
  #call 
  output{
  Array[File] sample_bam=WES_sample.bam
  Array[File] sample_bai=WES_sample.bai
  Array[File] sample_gvcf=WES_sample.gvcf
  Array[File] sample_fltindel=WES_sample.fltindel
  Array[File] sample_fltsnp=WES_sample.fltsnp
  Array[File] sample_fltvcf=WES_sample.fltvcf
  Array[File] sample_vep_snp=WES_sample.vep_snp
  Array[File] sample_vep_indel=WES_sample.vep_indel
  Array[File] sample_soapnuke_clean_read1=WES_sample.soapnuke_clean_read1
  Array[File] sample_soapnuke_clean_read2=WES_sample.soapnuke_clean_read2
  Array[File] sample_soapnuke_soapnuke_out=WES_sample.soapnuke_out
  Array[File] sample_soapnuke_soapnuke_stat=WES_sample.soapnuke_stat
  Array[File] sample_bam_stat_txt=WES_sample.bam_stat_txt
  Array[File] sample_bam_stat=WES_sample.bam_stat
  Array[File] sample_vcf_json=WES_sample.vcf_json
  Array[File] sample_anno_indel=WES_sample.anno_indel
  Array[File] sample_anno_snp=WES_sample.anno_snp
  }
} 
