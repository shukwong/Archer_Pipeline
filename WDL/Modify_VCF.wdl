version 1.0

# WDL Pipeline for modifying VCF
# -------------------------------
# This simple pipeline changes the SM tag of a VCF file using bcftools
# More options may be added later

# Created by: Wendy Wong
# Date: 07/27/2022

workflow modify_vcf_workflow {
    input {
        # Sequence Information
        String LIBRARY 
        String SAMPLE_NAME
		File VCF
        File VCF_index
	}

    call change_vcf_sm_tag {
        input:
            vcf = VCF,
            vcf_index = VCF_index,
            library = LIBRARY,
            sample_name = SAMPLE_NAME
    }


    output {
        File output_vcf = change_vcf_sm_tag.output_vcf
        File output_vcf_index = change_vcf_sm_tag.output_vcf_index
    }

}

task change_vcf_sm_tag {
    input {
        File vcf
        File vcf_index
        String library
        String sample_name
    }

    String vcf_base_name = basename(vcf)

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(vcf, "GB")

    runtime {
        docker: "dockerbiotools/bcftools:latest"
        memory: "16GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(2*data_size)
        disks: "local-disk ~{10 + round(2*data_size)} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<

        echo -e '~{library}\t~{sample_name}' >samples.txt 

        bcftools reheader --samples samples.txt -o  ~{vcf_base_name} ~{vcf} 
        tabix -p vcf ~{vcf_base_name}
    >>>

    output {
        File output_vcf = "~{vcf_base_name}"
        File output_vcf_index = "~{vcf_base_name}.tbi"
    }
}
