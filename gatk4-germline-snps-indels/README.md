### 【概述】
根据GATK最佳实践，haplotypecaller-gvcf-gatk4工作流在gvcf模式下在单个样本上运行GATK4HaploypeCaller工具。

执行时，工作流使用间隔列表文件对输入bam样本调用HaplotypeCaller工具。工作流产生的输出将是单个GVCF文件，然后可以与其他几个GVCF文件一起提供给GenomicsBimport，产生多样本VCF。haplotypecaller-gvcf-gatk4工作流默认gvcf模式在有效调用多个样本的变体时非常有用。

在为一个或几个示例调用变量时，可以将make_gvcf输入参数设置为false，使工作流直接调用变量并输出VCF文件。

### 【输入说明】

    "input_bam": "输入的bam文件，支持cram",
    "input_bam_index": "输入的BAM文件的索引",
    "ref_dict": "输入的参考基因组的字典文件",
    "ref_fasta": "输入的参考基因组文件",
    "ref_fasta_index": "输入的参考基因组的索引文件",
    "scattered_calling_intervals_list": "输入的区间列表文件",
