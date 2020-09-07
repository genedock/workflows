### fastqc 11.9
task fastqc{
        Array[File] read
        Int thread
        String name
        command {
                fastqc=/bioapp/FastQC_0_11_9/fastqc
                mkdir ${name}
                $fastqc -t ${thread} -o ${name} ${sep=' ' read}
                tar zcvf ${name}.tar.gz ${name}
        } 
        output{
                File Outfastqc="${name}.tar.gz"
        }
        runtime {
                docker: "seqflow/genedock_wgs:1.0"
                memory: "16G"
                disk: "800G"
                cpu: 4
        }
}
# 使用fastqc对fq文件质控，支持多个输入文件，结果以打包形式输出.
workflow FqQc{
  Array[File] fq
  Int thread=4
  String name="QC"
  call fastqc as qc{
    input:
      read=fq,
      thread=thread,
      name=name
  }
  output{
    qc.Outfastqc
  }
}