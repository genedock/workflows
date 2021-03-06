## Rapid\_Speed\_WGS\_GATK4

### 【概述】

本流程使用常用的BWA、GATK软件，对全基因组数据做变异查找分析，即从fastq到vcf的分析（包含SNV和Indel）。

### 【整体步骤】

该流程整体可分为4个部分:

1.序列比对（Mapping）：使用bwa mem进行序列比对，为了加快比对速度，聚道的计算平台将测序数据进行了分布式处理。

2.Bam处理（Bam processing）：比对后的bam文件进行排序和去重处理。最后基于染色体分组进行拆分，准备后续的分布式变异查找分析。 

3.SNV和Indel检测（SNP and Indel calling）：使用最新的GATK软件检测SNV和INDEL。

4.结果合并（Reduce）：将所有的vcf文件合并。

### 【软件版本】

bwa 0.7.17

samtools 1.10

picard 2.23.1

GATK 4.1.4.1

### 【运行时间】

约100G的WGS数据，运行时间约12小时20分。
