import "single/ExpansionhunterMToolBoxEhmitomerge.wdl" as wf_MT
import "single/cWGS_single_main.wdl" as wf_WGS
import "single/cWGS_single_main_nosplit.wdl" as wf_WGS_nosplit
import "single/Fastq_QC.wdl" as Fastq_QC
import "single/CnvSvLoh.wdl" as wf_CnvSvLoh
import "single/Annotation.wdl" as wf_anno
import "single/report.wdl" as report

workflow single{
  String sample_id
  #Array[File] lane_fq1=["genedockdx:/home/admin/wdl_test/rawdata/iw160R1.fq.gz"]
  #Array[File] lane_fq2=["genedockdx:/home/admin/wdl_test/rawdata/iw160R2.fq.gz"]
  Array[File] lane_fq1
  Array[File] lane_fq2
  String refname="hg19"
  File perl_script="genedockdx:/home/admin/script/markdup_bam_stat_nodepth.pl"
  Int split_line=50000000
  Int? compression=4
  Int? thread=4
  Int bin_size=100
  File config="genedockdx:/home/admin/script/WGS_data_config.tar.gz"
  String chrom="ALL"
  File exon_region="genedockdx:/home/admin/ref/exon.list"
  File region_depth="genedockdx:/home/admin/script/hg19_1Mb_region.txt"
  File intervals="genedockdx:/home/admin/ref/hg19_noN.bed"
  File otherBed="genedockdx:/home/admin/Database/hg19.other.bed"
  String enableRemoteReadRetrievalForInsertionsInGermlineCallingModes="1"
  File repeat_specs="genedockdx:/home/admin/Database/repeat_specs.tgz"
  File fam_mito_merge_py="genedockdx:/home/admin/Database/fam_tab_merge.py"
  String role="child"
  String sample_gender="Male"
  File anno_script="genedockdx:/home/admin/script/anno.sh"
  
   call wf_WGS.wf_WGS as sample_WGS{
    input:
    sample_id=sample_id,
    perl_script=perl_script,
    lane_fq1=lane_fq1,
    lane_fq2=lane_fq2,
    split_line=split_line,
    refname=refname,
    bin_size=bin_size,
    chrom=chrom,
    exon_region=exon_region,
    region_depth=region_depth,
    intervals=intervals,
    compression=compression,
    thread=thread,
    otherBed=otherBed,
    enableRemoteReadRetrievalForInsertionsInGermlineCallingModes=enableRemoteReadRetrievalForInsertionsInGermlineCallingModes,
  }
    call wf_MT.Expansionhunter_MToolBox_eh_mito_merge as sample_MT{
      input:
      Bam=sample_WGS.bam,
      Bai=sample_WGS.bai,
      repeat_specs=repeat_specs,
      sample_id=sample_id,
      refname=refname,
      fam_mito_merge_py=fam_mito_merge_py,
  }
 
  call Fastq_QC.FqQc as fq_qc{
    input:
      fq1=lane_fq1,
      fq2=lane_fq2,
  }
   call wf_CnvSvLoh.CnvSvLoh as CnvSvLoh_sample{
    input:
      Bam=sample_WGS.bam,
      Bai=sample_WGS.bai,
      Gvcf=sample_WGS.gvcf,
      Gvcftbi=sample_WGS.gvcf_tbi,
      fltvcf=sample_WGS.fltvcf,
      rawvcf=sample_WGS.rawvcf,
      sample_id=sample_id,
      refname=refname,
      sample_gender=sample_WGS.gender,
      cnvnator_txt=sample_WGS.cnvnator_txt,
      manta_vcf=sample_WGS.manta_vcf,
  }
  call wf_anno.Annotation as Annotation{
    input:
      vcf=sample_WGS.fltvcf,
      anno_script=anno_script,
      sample_id=sample_id,
  }
  call report.QC_report as report{
    input:
      fq1_size=fq_qc.qc1_fq_size,
      fq2_size=fq_qc.qc2_fq_size,
      bam_size=sample_WGS.bam_size,
      gvcf_size=sample_WGS.gvcf_size,
      gender=sample_WGS.gender,
      Auto=sample_WGS.auto,
      cnvnator_txt=sample_WGS.cnvnator_txt,
      manta_vcf=sample_WGS.manta_vcf,
      fltvcf=sample_WGS.fltvcf,
      FastQC=fq_qc.FQ_json,
      AlignmentQC=sample_WGS.AlignmentQC,
      workdir_anno=Annotation.annotation.workdir,
      role=role,
      sample_id=sample_id,
      erds_vcf=CnvSvLoh_sample.erds,
      config=config,
      refname=refname,
  }
  

  output{
  fq_qc.qc1_qc
  fq_qc.qc2_qc
  fq_qc.FQ_json
  sample_MT.Expansionhunter_MToolBox_sample.ehvcf
  sample_MT.Expansionhunter_MToolBox_sample.mt_anno
  sample_MT.Expansionhunter_MToolBox_sample.mt_tgz
  sample_MT.eh_mito_merge.all_eh_anno_tab
  sample_MT.eh_mito_merge.all_mito_anno_tab
  CnvSvLoh_sample.erds
  CnvSvLoh_sample.merge_cnv_erds
  CnvSvLoh_sample.merge_cnv_cnvnator
  CnvSvLoh_sample.sv_anno_all
  CnvSvLoh_sample.sv_anno_full
  CnvSvLoh_sample.cnv_json
  CnvSvLoh_sample.sv_json
  CnvSvLoh_sample.h3m2
  sample_WGS.bam
  sample_WGS.bai
  sample_WGS.gvcf
  sample_WGS.gvcf_tbi
  sample_WGS.bamstat
  sample_WGS.depth
  sample_WGS.AlignmentQC
  sample_WGS.bam_picard
  sample_WGS.cnvnator_txt
  sample_WGS.manta_vcf
  sample_WGS.VariantQC
  sample_WGS.rawvcf
  sample_WGS.fltvcf
  sample_WGS.snp_vcf
  sample_WGS.indel_vcf
  Annotation.annotation.workdir
  report.report_pdf
  }

}
