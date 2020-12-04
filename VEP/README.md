## VEP

###【概述】

使用Variant Effect Predictor软件对变异进行注释，标记变异对基因、转录组、蛋白序列和调控的影响。支持GRCh37/GRCh38。 

###【软件版本】

vep 101。

###【输入文件】

vcf文件。

###【输出文件】

注释后的vcf文件。

### 【参数】

VEP_annotation.refseq：可输入“Yes” 或其他，只有在“Yes”时使用refseq数据库注释，其他情况使用ensembl数据库注释。

VEP_annotation.assembly：支持GRCh37/GRCh38，需与输入vcf一致。

###【数据库及版本】

注释内容，参看http://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html。

GRCh37下额外提供了SPIDEX和dbscSNV数据库，需注意，SPIDEX禁止商业使用。
