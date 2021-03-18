## WES\_Somatic\_BWA-GATK-Mutect2\_ContainRef\_stat

###  【概述】

本流程使用常用的BWA、GATK、Mutect2等软件，对成对的外显子基因组数据，进行从fastq到vcf的分析（包含SNV和Indel) ，同时对原始序列 (FASTQ) 和比对后序列 (BAM) 的基本计量参数进行统计。

### 【整体步骤】

该流程整体可分为多个部分：

1. 序列比对（Mapping）：使用trimmomatic软件进行前处理，使用bwa mem进行比对，使用samtools对bam文件进行排序、格式转换等操作。此步骤tumor和normal样本独立进行。

2. Bam文件前处理之一（Bam processing I）：使用picard Markduplicate进行重复序列处理，使用GATK的BaseRecalibrator和PrintReads进行Base quality score recalibration。此步骤tumor和normal样本独立进行。

3. Bam文件前处理之二（Bam processing II）：使用GATK的RealignerTargetCreator和IndelRealigner模块进行Indel Realignment。此步骤会同时处理tumor和normal样本。

4. SNV、Indel检测（SNV calling & Indel calling）：使用软件Mutect2检测somatic SNV 和Indel突变。此步骤会同时处理tumor和normal样本。

5. CNV检测：使用软件CNVkit检测CNV，同时处理tumor和normal样本。

6. VEP注释：使用软件VEP注释变异检测结果。

7. report报告：基于对fastq，bam，vcf进行统计并自动生成报告。

### 【软件版本】

bwa 0.7.13

samtools 1.3

picard 2.23.1

GATK 3.8

VEP 101
### 【注意事项】

1. 本流程只适用于成对样本。

2. 本流程需要是双端测序。

3. 运行此公共工作流时，有些输入项我们会默认为用户绑定一些数据（主要是参考基因组和测试数据等），用户可以直接使用这些默认数据运行工作流，也可以删除或者替换这些数据后运行工作流。

4. 该公共工作流相关工具中内置了一些常用的reference、knowsites输入文件，您不需要自己上传，只需在对应参数设置中选择您需要的文件的对应名称，使用的reference、knowsites相关文件 均是Broad研究所的Resource bundle，下载自Broad的FTP：ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/。

5. 该流程中参数文件report_config暂时需要聚道为客户授权。
