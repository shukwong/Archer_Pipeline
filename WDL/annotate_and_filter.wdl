version 1.0

# Final WDL Pipeline for TERRAbio
# -------------------------------
# This pipeline is designed to process mutant/wildtype H.sapiens sequencing data from ArcherDX for low VAF variants.
# It features four variant callers (Mutect, Vardict, Lofreq, Pindel) for variant detection and
# performs various false positive filters and detection methods (fp_filter, PoN FishersTest, XGB Model).
# This pipeline also generates VEP style annotations for all called variants as well as additional putative driver
# annotations generated from various database sources (TOPMed, MSK-IMPACT, COSMIC, OncoKB, etc.)

## modify to just do annotate and xgb_model

# Created by: Irenaeus Chan
# Contact: chani@wustl.edu
# Date: 03/17/2022


# Main Workflow
workflow boltonlab_CH_annotate_filter {
    input {
        #sample info
        String tumor_sample_name
        String platform = "Illumina"

        #vcf files for PD annotations
        File mutect_vcf
        File lofreq_vcf
        File vardict_vcf
        String? pon_pvalue = "2.114164905e-6"

        # R Files for Putative Driver Annotations
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File TSG_file
        File oncoKB_curated
        File pd_annotation_file
        File truncating
        File cosmic_dir_zip

        # xgb_model inputs
        File pindel_vcf
        File pon_vcf
    }

   
    call annotatePD {
        input:
            mutect_vcf = mutect_vcf,
            lofreq_vcf = lofreq_vcf,
            vardict_vcf = vardict_vcf,
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            TSG_file = TSG_file,
            oncoKB_curated = oncoKB_curated,
            pd_annotation_file = pd_annotation_file,
            truncating = truncating,
            cosmic_dir_zip = cosmic_dir_zip,
            pon_pvalue = pon_pvalue
    }

    if (platform == 'ArcherDX') {
        call xgb_model {
            input:
            mutect_tsv = annotatePD.mutect_vcf_annotate_pd,
            lofreq_tsv = annotatePD.lofreq_vcf_annotate_pd,
            vardict_tsv = annotatePD.vardict_vcf_annotate_pd,
            pindel_full_vcf = pindel_vcf,
            pon = pon_vcf,
            pon_pvalue = pon_pvalue,
            model = true,
            tumor_sample_name = tumor_sample_name
        }
    }

    if (platform != 'ArcherDX') {
        call xgb_model as no_xgb_model {
            input:
            mutect_tsv = annotatePD.mutect_vcf_annotate_pd,
            lofreq_tsv = annotatePD.lofreq_vcf_annotate_pd,
            vardict_tsv = annotatePD.vardict_vcf_annotate_pd,
            pindel_full_vcf = pindel_vcf,
            pon = pon_vcf,
            pon_pvalue = pon_pvalue,
            model = false,
            tumor_sample_name = tumor_sample_name
        }
    }

    output {

        # R Things
        File mutect_annotate_pd = annotatePD.mutect_vcf_annotate_pd
        File lofreq_annotate_pd = annotatePD.lofreq_vcf_annotate_pd
        File vardict_annotate_pd = annotatePD.vardict_vcf_annotate_pd

        # Model
        File model_output = select_first([xgb_model.model_output, no_xgb_model.model_output])
        File model_raw_output = select_first([xgb_model.model_raw_output, no_xgb_model.model_raw_output])
        File mutect_complex = select_first([xgb_model.mutect_complex, no_xgb_model.mutect_complex])
        File pindel_complex = select_first([xgb_model.pindel_complex, no_xgb_model.pindel_complex])
        File lofreq_complex = select_first([xgb_model.lofreq_complex, no_xgb_model.lofreq_complex])
        File caller_filters = select_first([xgb_model.caller_filters, no_xgb_model.caller_filters])
    }
}






task annotatePD {
    input {
        File mutect_vcf
        File lofreq_vcf
        File vardict_vcf
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File TSG_file
        File oncoKB_curated
        File pd_annotation_file
        File truncating
        File cosmic_dir_zip
        String? pon_pvalue = "2.114164905e-6"
    }

    Float caller_size = size([mutect_vcf, lofreq_vcf, vardict_vcf], "GB")
    Float file_size = size([bolton_bick_vars, mut2_bick, mut2_kelly, matches2, TSG_file, oncoKB_curated, pd_annotation_file, truncating], "GB")
    Float cosmic_size = 3*size(cosmic_dir_zip, "GB")
    Int space_needed_gb = 20 + round(caller_size + file_size + cosmic_size)
    Int cores = 2
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/r_docker_ichan:latest"
      memory: "12GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String cosmic_dir = basename(cosmic_dir_zip, ".zip")

    command <<<
        set -eou pipefail

        unzip -qq ~{cosmic_dir_zip}

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{mutect_vcf} --out $(basename ~{mutect_vcf} .vcf.gz) --caller mutect \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
        echo "Mutect AnnotatePD Finished..."

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{lofreq_vcf} --out $(basename ~{lofreq_vcf} .vcf.gz) --caller lofreq \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
        echo "Lofreq AnnotatePD Finished..."

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{vardict_vcf} --out $(basename ~{vardict_vcf} .vcf.gz) --caller vardict \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
        echo "Vardict AnnotatePD Finished..."
    >>>

    output {
        File mutect_vcf_annotate_pd = basename(mutect_vcf, ".vcf.gz") + ".tsv"
        File lofreq_vcf_annotate_pd = basename(lofreq_vcf, ".vcf.gz") + ".tsv"
        File vardict_vcf_annotate_pd = basename(vardict_vcf, ".vcf.gz") + ".tsv"
    }
}

# Generates an output file that combines all repeated columns (features that are not unique to each caller)
# Adds additional features and runs select features through a trained (on ArcherDX Samples) XGBoost Model.
# False positive filters are run to tag additional artifacts from the final list,

# Additionally performs complex variant matching by looking for complex variants called by Vardict and matching possible
# components found in Mutect or Lofreq. If complex variants are identified they are matched with Pindel as a final stop
# to determine the validity of real complex variants.

# Additional features are added such as the Z-score for the number of reference reads for all variants in a caller.
# The Z-score tells us whether the alignment step correctly worked for all positions as one would except the reference depth to be fairly consistent for all variants.
# An outlier indicates that this region may be over or under aligned. Another feature is the Fisher’s test with the forward and reverse reference / alternate reads.
# A significant value here would indicate that the forward and reverse strands are not consistent and strand bias is occurring. We were not able to calculate a cutoff p-value so this feature is not used in our analysis.

# Finally a predicted probabiliy from the XGBoost Model is calculated for each variant
task xgb_model {
    input {
        File mutect_tsv
        File lofreq_tsv
        File vardict_tsv
        File pindel_full_vcf
        File pon
        String? pon_pvalue = "2.114164905e-6"
        Boolean model = false
        String tumor_sample_name
    }

    Float caller_size = size([mutect_tsv, lofreq_tsv, vardict_tsv, pindel_full_vcf, pon], "GB")
    Int space_needed_gb = 10 + round(caller_size)
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/xgb"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        if [ ~{model} == true ]; then
            /opt/bin/xgbappcompiled/bin/xgbapp ~{lofreq_tsv} ~{mutect_tsv} ~{vardict_tsv} ~{pindel_full_vcf} ~{pon} --pvalue ~{pon_pvalue}
        else
            /opt/bin/xgbappcompiled/bin/xgbapp ~{lofreq_tsv} ~{mutect_tsv} ~{vardict_tsv} ~{pindel_full_vcf} ~{pon} --pvalue ~{pon_pvalue} --nomodel
        fi
        echo "Model Finished..."
    >>>


    output {
        File model_output = "output_~{tumor_sample_name}.tsv.gz"
        File model_raw_output = "output_~{tumor_sample_name}.raw.tsv.gz"
        File mutect_complex = "output_mutect_complex_~{tumor_sample_name}.tsv.gz"
        File pindel_complex = "output_pindel_complex_~{tumor_sample_name}.tsv.gz"
        File lofreq_complex = "output_lofreq_complex_~{tumor_sample_name}.tsv.gz"
        File caller_filters = "Caller_Filters.raw.tsv.gz"
    }
}
