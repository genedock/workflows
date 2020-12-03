## Metagenome\_Prinseq-BMtagger-HUMAnN2

### 【概述】
本流程使用Prinseq，BMtagger和HUMAnN2等软件，对人相关宏基因组测序的fastq.gz文件进行数据过滤，去重复，去掉疑似人的reads，最后通过HUMAnN2软件分析该样品的物种丰度情况，功能基因家族(Gene families)丰度情况和代谢通路丰度信息(Pathway abundance)和覆盖度（Pathway coverage) 信息。

### 【基本背景】

宏基因组测序（whole metagenome shotgun sequencing）是指使用二代测序技术（NGS）对环境样品的微生物进行测序。相比扩增子类测序（targeted metagenome sequencing），比如16S等,宏基因组测序能发掘环境内的所有微生物的物种丰度信息以及对应基因家族丰度信息和代谢通路信息。

### 【整体步骤】

该流程整体可以分为3部分：

1. 序列过滤和质控（QC）：使用Prinseq软件对测序数据进行处理，过滤低质量的reads和重复reads。

2. 过滤宿主信息（BMtagger）：使用BMtagger与人基因组序列进行比对，去除疑似为人序列的reads。

3. 物种丰度计算，基因丰度计算和代谢通路计算（HUMAnN2) ：使用HUMAnN2计算样品里的微生物物种丰度信息，已经相关基因的丰度信息，代谢通路的丰度信息。最后得到物种丰度的tsv表格，基因丰度的tsv表格，代谢通路丰度的tsv以及该代谢通路的覆盖度tsv文件。

### 【运行时间】

1.45G+1.46G 的fastq.gz文件，北京域：运行时间大约3小时。

### 【准确性评估】

本流程依据HMP项目进行软件选择和参数选择。Prinseq，BMtagger和HUMAnN2这三个软件都是HMP项目指定的软件。

Pringseq homepage : http://prinseq.sourceforge.net/index.html

软件地址和HMP SOP如下:

BMtagger URL: http://hmpdacc.org/doc/HumanSequenceRemoval\_SOP.pdf

ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/

HUMAnN2 URL: http://hmpdacc.org/doc/HUMAnN\_SOP.pdf

https://bitbucket.org/biobakery/humann2/wiki/Home

### 【注意事项】

1. 本流程只适用于illumina双端测序的文件，输入文件是fastq.gz文件。

2. 本流程只适用于人环境相关的，比如人肠道，人口腔等，不适用于其他宿主和其他环境的样品。

3. 运行此公共工作流时，有些输入项我们会默认为用户绑定一些数据（主要是参考基因组,数据库和测试数据等)。 用户可以直接使用这些默认数据运行工作流，也可以删除或者替换测序数据后运行工作流。
