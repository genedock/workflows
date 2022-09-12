### 【概述】
该WDL工作流根据人类外显子组测序数据中种系SNP和Indel发现的GATK最佳实践实现数据预处理和初始variant calling。

### 【输入说明】
    {
      "WholeGenomeGermlineSingleSample.sample_and_unmapped_bams": {
        "sample_name": "输入的样本名称",
        "base_file_name": "输入的样本名称",
        "flowcell_unmapped_bams": "输入的未映射BAM文件列表",
        "final_gvcf_base_name": "输出的GVCF文件的基本名称",
        "unmapped_bam_suffix": "未映射BAM文件的后缀",
      },
    
      "WholeGenomeGermlineSingleSample.references": {
        "contamination_sites_ud": "污染位点UD文件",
        "contamination_sites_bed": "污染位点BED文件",
        "contamination_sites_mu": "污染位点MU文件",
        "calling_interval_list": "调用区间列表",
        "reference_fasta": {
          "ref_dict": "输入的参考基因组的字典文件",
          "ref_fasta": "输入的参考基因组文件",
          "ref_fasta_index": "输入的参考基因组的索引文件",
          "ref_alt": "输入的参考基因组的ALT文件",
          "ref_sa": "输入的参考基因组的SA文件",
          "ref_amb": "输入的参考基因组的AMB文件",
          "ref_bwt": "输入的参考基因组的BWT文件",
          "ref_ann": "输入的参考基因组的ANN文件",
          "ref_pac": "输入的参考基因组的PAC文件",
        },
        "known_indels_sites_vcfs": "已知的Indels位点VCF文件列表",
        "known_indels_sites_indices": "已知的Indels位点VCF文件的索引列表",
        "dbsnp_vcf": "dbSNP VCF文件",
        "dbsnp_vcf_index": "dbSNP VCF索引文件",
        "evaluation_interval_list": "评估区间列表",
        "haplotype_database_file": "Haplotype数据库文件",
      },
    
      "WholeGenomeGermlineSingleSample.scatter_settings": {
        "haplotype_scatter_count": 50,
        "break_bands_at_multiples_of": 1000000
      },
    
      "WholeGenomeGermlineSingleSample.fingerprint_genotypes_file": "指纹基因型文件",
      "WholeGenomeGermlineSingleSample.wgs_coverage_interval_list": "WGS覆盖区间列表",
    
      "WholeGenomeGermlineSingleSample.papi_settings": {
        "preemptible_tries": 3,
        "agg_preemptible_tries": 3
      }
    }