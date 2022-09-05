version 1.0
## This WDL performs format validation on SAM/BAM files in a list
##
## Requirements/expectations :
## - One or more SAM or BAM files to validate
## - Explicit request of either SUMMARY or VERBOSE mode in inputs.json
##
## Outputs:
## - Set of .txt files containing the validation reports, one per input file
##
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##

# WORKFLOW DEFINITION
workflow ValidateBamsWf {
  input {
    Array[File] bam_array
    String gatk_docker = "genedockdx/gatk:latest"
    String gatk_path = "/gatk/gatk"
  }

  # Process the input files in parallel
  scatter (input_bam in bam_array) {

    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(input_bam, ".bam")

    # Run the validation
    call ValidateBAM {
      input:
        input_bam = input_bam,
        output_basename = bam_basename + ".validation",
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] validation_reports = ValidateBAM.validation_report
  }
}

# TASK DEFINITIONS

# Validate a SAM or BAM using Picard ValidateSamFile
task ValidateBAM {
  input {
    # Command parameters
    File input_bam
    String output_basename
    String? validation_mode
    String gatk_path

    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }

  Int disk_size = ceil(size(input_bam, "GB")) + addtional_disk_space_gb
  String output_name = "${output_basename}_${validation_mode}.txt"

  command {
    ${gatk_path} \
      ValidateSamFile \
      --INPUT ${input_bam} \
      --OUTPUT ${output_name} \
      --MODE ${default="SUMMARY" validation_mode}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File validation_report = "${output_name}"
  }
}