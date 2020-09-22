# workflows
FastqQC工作流：使用fastqc对fq文件质控.  
WES工作流：包括数据比对、排序、去重、BQSR、callgvcf、vcfFilter，支持3种参考基因组，使用bwa 7.17+GATK4；
流程支持分块，可在运行时候选择，默认按5千万行拆分.
Sentieon_WES工作流：使用Trimmomatic数据过滤，包括比对、排序、sentieon变异检测、fq和bam质控。
