# WDL版本:single_WGS version 1

包括数据比对、排序、去重、BQSR、callgvcf、vcfFilter、VCF注释、线粒体、CNV、SV、LOH以及质控报告。支持参考基因组hg19/,使用bwa 7.17+GATK3.8；
流程支持分块(默认分块)，可在运行时候选择，默认按5千万行拆分. 
使用前请联系聚道工作人员，申请report_config配置文件(报告)授权以及其他配置文件授权。 
