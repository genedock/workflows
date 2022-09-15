### 【概述】
使用GATK4进行体细胞短变异分析的工作流
### 【输入说明】
    "Mutect2.gatk_docker": "GATK4的docker镜像",
    "Mutect2.intervals": "输入的目标区域文件",
    "Mutect2.scatter_count": "并行数量",
    "Mutect2.m2_extra_args": "Mutect2的附加参数字符串",
    "Mutect2.filter_funcotations": "是否过滤Funcotator的注释结果",
    "Mutect2.funco_reference_version": "Funcotator的参考基因组版本",
    "Mutect2.funco_data_sources_tar_gz": "Funcotator的数据源",
    "Mutect2.funco_transcript_selection_list": "Funcotator的转录本选择列表",

    "Mutect2.ref_fasta": "输入的参考基因组文件",
    "Mutect2.ref_dict": "输入的参考基因组的字典文件",
    "Mutect2.ref_fai": "输入的参考基因组的索引文件",
    "Mutect2.normal_reads": "输入的正常样本的BAM文件",
    "Mutect2.normal_reads_index": "输入的正常样本的BAM文件的索引文件",
    "Mutect2.tumor_reads": "输入的肿瘤样本的BAM文件",
    "Mutect2.tumor_reads_index": "输入的肿瘤样本的BAM文件的索引文件",

    "Mutect2.pon": "预过滤变异位点集VCF文件",
    "Mutect2.pon_idx": "预过滤变异位点集VCF文件的索引文件",
    "Mutect2.gnomad": "指定人群的germline信息VCF文件",
    "Mutect2.gnomad_idx": "指定人群的germline信息VCF文件的索引文件",
    "Mutect2.variants_for_contamination": "用于计算污染的具有等位基因频率的常见变异的VCF文件",
    "Mutect2.variants_for_contamination_idx": "用于计算污染的具有等位基因频率的常见变异的VCF文件的索引文件",
    "Mutect2.realignment_index_bundle": "GATK BWA 镜像文件，使用BwaMemIndexImageCreator生成"
