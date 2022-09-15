### 【概述】
体细胞拷贝数变体分析工作流
### 【输入说明】
    "CNVSomaticPairWorkflow.tumor_bam_idx": "输入的肿瘤样本BAM文件的索引文件",
    "CNVSomaticPairWorkflow.normal_bam": "输入的正常样本BAM文件",
    "CNVSomaticPairWorkflow.normal_bam_idx": "输入的正常样本BAM文件的索引文件",
    "CNVSomaticPairWorkflow.tumor_bam": "输入的肿瘤样本BAM文件",

    "CNVSomaticPairWorkflow.ref_fasta": "输入的参考基因组文件",
    "CNVSomaticPairWorkflow.ref_fasta_fai": "输入的参考基因组的索引文件",
    "CNVSomaticPairWorkflow.ref_fasta_dict": "输入的参考基因组的字典文件",
    "CNVSomaticPairWorkflow.common_sites": "输入的公共变异位点文件",
    "CNVSomaticPairWorkflow.read_count_pon": "输入的正常样本的ReadCount文件",
    "CNVSomaticPairWorkflow.intervals": "输入的目标区域文件",

    "CNVSomaticPairWorkflow.gatk_docker": "GATK4的docker镜像",

    "CNVSomaticPairWorkflow.is_run_funcotator": "是否运行Funcotator",
    "CNVSomaticPairWorkflow.funcotator_data_sources_tar_gz": "Funcotator的数据源文件",
    "CNVSomaticPairWorkflow.funcotator_transcript_selection_list": "Funcotator的转录本选择列表",
