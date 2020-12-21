## Nanopore

### 【概述】
本流程对三代的nanopore数据进行SV分析和点突变分析，使用的软件和步骤与Nanopore官方的EPI2ME平台上流程



### 【整体步骤】
主要包括如下步骤：

1、使用fastp进行fq文件分块；

2、使用minimap2进行比对；

3、使用sniffles进行SV检测；

4、使用AnnotSV对SV结果进行注释；

5、使用medaka_variant进行点突变检测。




### 【运行时间】
测试93G fq文件（public:/demo-data/nanopore/HG002_NOT.fastq.gz）。北京域：运行21h（点突变超时暂停）



### 【软件版本】
fastp   0.21.0
minimap2    2.17
sniffles    1.0.12
annotsv 2.2
medaka-variant  1.0.3



### 【注意事项】
镜像参考https://github.com/nanoporetech/pipeline-structural-variation
