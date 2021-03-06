## WES\_Germline\_BWA\_GATK4\_multi

### 【概述】

本流程使用常用的BWA、GATK软件，对外显子组配对测序数据进行变异查找分析，即从fastq到vcf的分析（包含SNV和Indel）。

### 【整体步骤】

该流程整体可分为以下几个部分:

1. 序列比对（Mapping）：适用soapnuke进行比对前处理，对clean reads使用bwa mem进行序列比对。

2. Bam处理（Bam processing）：比对后的bam文件进行去重、BQSR处理。最后基于染色体分组进行拆分，准备后续的分布式变异查找分析。 

3. SNV和Indel检测（SNP and Indel calling）：使用GATK软件（GATK-4.1.4.1）HaplotypeCaller模块检测SNV和INDEL。

4. 过滤(hard filter)：适用GATK 4.1.4.1基于常用过滤阈值对vcf进行hard filter。

5. 注释(annotation)：基于VEP分别注释SNP与INDEL。使用--pick参数，仅保留一条注释转录本。 

6. 指标：输出指标包括soapnuke质控指标统计、bam的质控指标统计、vcf文件的质控指标统计以及VEP注释的SNP/INDEL突变类型统计。
### 【软件版本】

bwa 0.7.17

samtools 1.10

picard 2.23.1

GATK 4.1.4.1

VEP 101

### 【注意事项】
1. 该WES多样本流程允许单样本存在多lane数据，但比对前就已进行数据合并过滤，分lane信息后续不予保留。

2. 该多样本WES流程基于WES单样本流程，质控统计指标未进行统一合并。

3. 由于复杂数据结构难以在前端展示，该流程目前需要自行生成input.json，再基于gwc workflow run提交（可借鉴同级目录下的input.json,仅需对应修改lane_fq1,lane_fq2,bed,sample_lb,sample_id）。

4. 部分脚本内置于聚道系统，使用需联系聚道生信人员。

