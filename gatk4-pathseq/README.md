### 【概述】
使用PathSeq进行计算病原体发现的工作流程，PathSeq是GATK中的一个pipeline，用于检测从宿主生物体获取的测序样本中的微生物
### 【输入说明】
    "PathSeqPipelineWorkflow.input_bam": "输入的BAM文件",
    "PathSeqPipelineWorkflow.sample_name": "输入的BAM文件的样本名称",

    "PathSeqPipelineWorkflow.kmer_file": "输入的Kmer文件",
    "PathSeqPipelineWorkflow.microbe_bwa_image": "输入的Microbe BWA Image文件",
    "PathSeqPipelineWorkflow.microbe_fasta": "输入的Microbe Fasta文件",
    "PathSeqPipelineWorkflow.microbe_fasta_dict": "输入的Microbe Fasta字典文件",
    "PathSeqPipelineWorkflow.filter_bwa_image": "输入的Filter BWA Image文件",
    "PathSeqPipelineWorkflow.taxonomy_file": "输入的Taxonomy文件",
    "PathSeqPipelineWorkflow.is_host_aligned": "是否对宿主进行比对",

    "PathSeqPipelineWorkflow.gatk_docker_override": "GATK Docker镜像",
    "#PathSeqPipelineWorkflow.gatk4_jar_override": "GATK4 Jar包",

    "PathSeqPipelineWorkflow.filter_duplicates": "是否过滤重复序列",
    "#PathSeqPipelineWorkflow.min_clipped_read_length": "最小的剪接比对Read长度",
    "#PathSeqPipelineWorkflow.filter_bwa_seed_length": "过滤 BWA 种子长度",
    "#PathSeqPipelineWorkflow.min_score_identity": "最小的比对得分",
    "#PathSeqPipelineWorkflow.host_min_identity": "宿主最小比对得分",
    "#PathSeqPipelineWorkflow.PathseqPipeline.skip_quality_filters": "是否跳过质量过滤",
    "#PathSeqPipelineWorkflow.PathseqPipeline.default_min_score_identity": "默认的最小比对得分",
    "#PathSeqPipelineWorkflow.identity_margin": "比对得分的边界",
    "#PathSeqPipelineWorkflow.skip_pre_bwa_repartition": "是否跳过BWA重分区",
    "PathSeqPipelineWorkflow.divide_by_genome_length": "是否根据基因组长度进行分区",
    "#PathSeqPipelineWorkflow.bam_partition_size": "BAM分区大小",
    "#PathSeqPipelineWorkflow.PathseqPipeline.default_identity_margin": "默认的比对得分边界",
    "#PathSeqPipelineWorkflow.PathseqPipeline.default_ram_mb": "默认的内存大小",

    "#PathSeqPipelineWorkflow.mem_gb": "内存大小",
    "#PathSeqPipelineWorkflow.cpu": "CPU核数",
    "#PathSeqPipelineWorkflow.increase_disk_size": "是否增加磁盘大小",
    "#PathSeqPipelineWorkflow.PathseqPipeline.default_disk_space_gb": "默认的磁盘大小",
    "#PathSeqPipelineWorkflow.PathseqPipeline.use_ssd": "是否使用SSD",