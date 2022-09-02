### 【概述】
该工作流将CRAM文件转换为BAM文件，并输出BAM、BAM索引和验证报告。

之所以选择这种方法而不是直接使用Samtools将CRAM转换为BAM，是因为Samtools 1.3由于包中包含的htslib的旧版本而产生不正确的二进制文件。
### 【输入说明】
    “input_bam”：“输入的CRAM文件”