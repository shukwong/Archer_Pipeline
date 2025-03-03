version 1.0

#workflow to scatter samples to run on Archer pipeline 

import "./terra_pipeline.wdl" as archer_pipeline

workflow run_archer_by_sample {
	input {
		File sample_tsv_file

		# Sequence Information
        String platform = "ArcherDX"
        String platform_unit = "Illumina"
        String library = "LIBRARY"
        File fastq_one
        File fastq_two
        Boolean bam_input = true           
        #File? unaligned_bam = ""    
        Boolean is_umi_concensus_unaligned = true       #using consensus bam (unaligned) as input (DEFAULT)
        #File? consensus_unaligned_bam               
        Boolean? aligned = false            # If input is an already aligned BAM file then set this flag
        File? aligned_bam_file
        File? aligned_bam_file_bai
        Array[String] read_structure        # Used for the UMI processing see: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
        #String tumor_sample_name 
        File target_intervals               # Interval List
        Int? mem_limit_override = 6         # Some applications will require more memory depending on BAM size and BED size... (in GB)
                                            # Need to account for these types of errors
                                      
        # Reference
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa

        # FASTQ Preprocessing
        Boolean has_umi = true
        Boolean? umi_paired = true          # If the UMI is paired (R1 and R2) then set this flag
        Int umi_length = 8
        Array[Int] min_reads = [1]          # The minimum number of reads that constitutes a "read family"
        Float? max_read_error_rate = 0.05   # If 5% of the reads within a "read family" does not match, removes the family
        Float? max_base_error_rate = 0.1    # If 10% of the bases are different than the majority, masks the base with N
        Int min_base_quality = 1            # Any base less than QUAL 1 will be masked with N
        Float max_no_call_fraction = 0.5    # Maximum fraction of no-calls (N) in the read after filtering

		# BQSR
        Array[File] bqsr_known_sites
        Array[File] bqsr_known_sites_tbi
        Array[String] bqsr_intervals

        # QC
        # These are essentially coordinates for regions you were able to design probes for in the reagent.
        # Typically the reagent provider has this information available in bed format and it can be converted to an interval_list with Picard BedToIntervalList.
        # Astrazeneca also maintains a repo of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
        File bait_intervals
        Array[LabelledFile] per_base_intervals      # If QC needs to be done at a per-base resolution
        Array[LabelledFile] per_target_intervals    # If QC needs to be done at a per-target resolution
        Array[LabelledFile] summary_intervals       # If QC needs to be done for specific intervals
        File omni_vcf                               # The omni VCF is a list of sites used by verifyBamId for identifying contamination (really mixing of multiple samples)
        File omni_vcf_tbi
        String picard_metric_accumulation_level
        Int? qc_minimum_mapping_quality = 0
        Int? qc_minimum_base_quality = 0
        File chrom_sizes
        File af_only_snp_only_vcf

        # Variant Calling
        Array[Pair[File, File]] pon_bams            # List of BAMs within the Panel of Normals (PoN)
        Float? af_threshold = 0.0001                # Minimum VAF Cut-Off
        String? pon_pvalue = "2.114164905e-6"       # Bonferroni Corrected P-Value for Significance

        # Normal BAMs
        Boolean tumor_only = true                   # Defines if Normal BAMs should be used for Variant Calling
        File normal_bam                             # TODO: Implement Normal
        File normal_bam_bai

        # Pindel
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3

        # See: http://bcb.io/2016/04/04/vardict-filtering/
        # Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

		# PoN2
        # If a variant exists inside two or more of our Panel of Normal Samples (PoN) at 2% VAF or greater, then this variant is assumed to be
        # sequencing noise because our PoN should not have any variants (that are not germline).
        # Instructions:
        # PoN Samples are run through Tumor Only mode and Filtered for each caller.
        # Variants below 2% VAF are removed from PoN Samples
        # Variants that appear within the PoN Samples in two or more samples are kept
        File mutect_pon2_file
        File mutect_pon2_file_tbi
        File lofreq_pon2_file
        File lofreq_pon2_file_tbi
        File vardict_pon2_file
        File vardict_pon2_file_tbi

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

        # VEP Parameters
        File vep_cache_dir_zip                      # WDL does not have a Directory Variable, so the entire cache needs to be ZIP
        String vep_ensembl_assembly
        String vep_ensembl_version
        String vep_ensembl_species
        Array[String] vep_plugins = ["Frameshift", "Wildtype"]
        VepSpliceAIPlugin vep_plugin_spliceAI_files = {}
        File? synonyms_file
        Boolean? annotate_coding_only = true
        Array[VepCustomAnnotation] vep_custom_annotations
        String vep_pick = "pick"
        Boolean everything = true

        # Variants within gnomAD that have a VAF higher than 0.5%
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi

        # Used for Pilot Data. Leave FALSE
        Boolean? pilot = false

	}

    # read table with sample information, we assume that 
    Array[Array[String]] samples_tsv = read_tsv(sample_tsv_file)
    Int n_samples = length(samples_tsv)-1
    #scatter (idx in range(n_samples)) { Array[String] samples_tsv_rows = samples_tsv[(idx+1)] }
    #Map[String, Array[String]] samples_tbl = as_map(zip(samples_tsv[0], transpose(samples_tsv_rows)))

    # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
    scatter (idx in range(n_samples)) {
    	
		#String consensus_unaligned_bam = batch_tbl["bam"][idx]
		#String tumor_sample_name = batch_tbl["tumor_sample_name"][idx]

		String consensus_unaligned_bam = samples_tsv[(idx+1)][3]
		String tumor_sample_name = samples_tsv[(idx+1)][1]

		call archer_pipeline.boltonlab_CH as archer_workflow {
			input:
				# Sequence Information
        		platform = platform,
        	    platform_unit = platform_unit,
         		library = library,
        		fastq_one = fastq_one,
        		fastq_two = fastq_two,
        		bam_input = bam_input,           
        		unaligned_bam = "<unaligned_bam>",    
        	 	is_umi_concensus_unaligned = is_umi_concensus_unaligned,       #using consensus bam (unaligned) as input (DEFAULT)
        		consensus_unaligned_bam = consensus_unaligned_bam,               
        		aligned = aligned,            # If input is an already aligned BAM file then set this flag
        		#File? aligned_bam_file
        		#File? aligned_bam_file_bai
        		read_structure = read_structure,       # Used for the UMI processing see: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
        		tumor_sample_name = tumor_sample_name,
         		target_intervals = target_intervals,               # Interval List
        		mem_limit_override = mem_limit_override,         # Some applications will require more memory depending on BAM size and BED size... (in GB)
                                            # Need to account for these types of errors
                                      
        		# Reference
        	 	reference = reference,
         		reference_fai = reference_fai,
         		reference_dict = reference_dict,
         		reference_amb = reference_amb,
         		reference_ann = reference_ann,
         		reference_bwt = reference_bwt,
         		reference_pac = reference_pac,
         		reference_sa = reference_sa,

        		# FASTQ Preprocessing
         		has_umi = has_umi,
        	    umi_paired = umi_paired,          # If the UMI is paired (R1 and R2) then set this flag
        		umi_length = umi_length,
        		min_reads = min_reads,          # The minimum number of reads that constitutes a "read family"
        		max_read_error_rate = max_read_error_rate,   # If 5% of the reads within a "read family" does not match, removes the family
        		max_base_error_rate = max_base_error_rate,   # If 10% of the bases are different than the majority, masks the base with N
        		min_base_quality = min_base_quality,            # Any base less than QUAL 1 will be masked with N
        		max_no_call_fraction = max_no_call_fraction,    # Maximum fraction of no-calls (N) in the read after filtering

				# BQSR
        		bqsr_known_sites = bqsr_known_sites,
        		bqsr_known_sites_tbi = bqsr_known_sites_tbi,
        		bqsr_intervals = bqsr_intervals,

        		# QC
        		# These are essentially coordinates for regions you were able to design probes for in the reagent.
        		# Typically the reagent provider has this information available in bed format and it can be converted to an interval_list with Picard BedToIntervalList.
        		# Astrazeneca also maintains a repo of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
        		bait_intervals = bait_intervals,
        		per_base_intervals = per_base_intervals,     # If QC needs to be done at a per-base resolution
        		per_target_intervals = per_target_intervals,   # If QC needs to be done at a per-target resolution
        		summary_intervals = summary_intervals,       # If QC needs to be done for specific intervals
         		omni_vcf = omni_vcf,                                # The omni VCF is a list of sites used by verifyBamId for identifying contamination (really mixing of multiple samples)
         		omni_vcf_tbi = omni_vcf_tbi,
         		picard_metric_accumulation_level = picard_metric_accumulation_level,
        		qc_minimum_mapping_quality = qc_minimum_mapping_quality,
        		qc_minimum_base_quality = qc_minimum_base_quality,
         		chrom_sizes = chrom_sizes,
         		af_only_snp_only_vcf = af_only_snp_only_vcf,

        		# Variant Calling
         		pon_bams = pon_bams,            # List of BAMs within the Panel of Normals (PoN)
         		af_threshold = af_threshold,            # Minimum VAF Cut-Off
         		pon_pvalue = pon_pvalue,   # Bonferroni Corrected P-Value for Significance

        		# Normal BAMs
         		tumor_only = tumor_only,                   # Defines if Normal BAMs should be used for Variant Calling
         		normal_bam =  normal_bam,                            # TODO: Implement Normal
         		normal_bam_bai = normal_bam_bai,

        		# Pindel
         		pindel_insert_size = pindel_insert_size,
        		ref_name = ref_name,
         		ref_date = ref_date,
         		pindel_min_supporting_reads = pindel_min_supporting_reads,

        		# See: http://bcb.io/2016/04/04/vardict-filtering/
        		# Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        		bcbio_filter_string = bcbio_filter_string,

				# PoN2
        		# If a variant exists inside two or more of our Panel of Normal Samples (PoN) at 2% VAF or greater, then this variant is assumed to be
        		# sequencing noise because our PoN should not have any variants (that are not germline).
        		# Instructions:
        		# PoN Samples are run through Tumor Only mode and Filtered for each caller.
        		# Variants below 2% VAF are removed from PoN Samples
        		# Variants that appear within the PoN Samples in two or more samples are kept
         		mutect_pon2_file = mutect_pon2_file,
         		mutect_pon2_file_tbi = mutect_pon2_file_tbi,
         		lofreq_pon2_file = lofreq_pon2_file,
         		lofreq_pon2_file_tbi = lofreq_pon2_file_tbi,
         		vardict_pon2_file = vardict_pon2_file,
         		vardict_pon2_file_tbi = vardict_pon2_file_tbi,

        		# R Files for Putative Driver Annotations
         		bolton_bick_vars = bolton_bick_vars,
         		mut2_bick = mut2_bick,
         		mut2_kelly = mut2_kelly,
         		matches2 = matches2,
         		TSG_file = TSG_file,
         		oncoKB_curated = oncoKB_curated,
         		pd_annotation_file = pd_annotation_file,
         		truncating = truncating,
         		cosmic_dir_zip = cosmic_dir_zip,

        		# VEP Parameters
        		vep_cache_dir_zip = vep_cache_dir_zip,                     # WDL does not have a Directory Variable, so the entire cache needs to be ZIP
         		vep_ensembl_assembly = vep_ensembl_assembly,
         		vep_ensembl_version = vep_ensembl_version,
         		vep_ensembl_species = vep_ensembl_species,
        		vep_plugins = vep_plugins,
         		vep_plugin_spliceAI_files = vep_plugin_spliceAI_files,
         		synonyms_file = synonyms_file,
         		annotate_coding_only = annotate_coding_only,
         		vep_custom_annotations = vep_custom_annotations,
         		vep_pick = vep_pick,
         		everything = everything,

        		# Variants within gnomAD that have a VAF higher than 0.5%
         		normalized_gnomad_exclude = normalized_gnomad_exclude,
         		normalized_gnomad_exclude_tbi = normalized_gnomad_exclude_tbi
		}

    }

	output {
        # Alignments
        Array[File?] aligned_bam = archer_workflow.aligned_bam 
        Array[File]? bqsr_bam = archer_workflow.bqsr_bam

        # Tumor QC
        Array[File]? tumor_insert_size_metrics = archer_workflow.tumor_insert_size_metrics
        Array[File]? tumor_hs_metrics =  archer_workflow.tumor_hs_metrics
        Array[Array[File]]? tumor_per_target_coverage_metrics = archer_workflow.tumor_per_target_coverage_metrics
        Array[Array[File]]? tumor_per_target_hs_metrics = archer_workflow.tumor_per_target_hs_metrics
        Array[Array[File]]? tumor_per_base_coverage_metrics = archer_workflow.tumor_per_base_coverage_metrics
        Array[Array[File]]? tumor_per_base_hs_metrics = archer_workflow.tumor_per_base_hs_metrics
        Array[Array[File]]? tumor_summary_hs_metrics = archer_workflow.tumor_summary_hs_metrics
        Array[File]? tumor_flagstats = archer_workflow.tumor_flagstats
        Array[File]? tumor_verify_bam_id_metrics = archer_workflow.tumor_verify_bam_id_metrics
        Array[File]? tumor_verify_bam_id_depth = archer_workflow.tumor_verify_bam_id_depth
        Array[File]? fastqc_html = archer_workflow.fastqc_html
        Array[File]? fastqc = archer_workflow.fastqc
        #File somalier_out = somalier.somalier_out


        # Mutect
        Array[File] mutect_full =  archer_workflow.mutect_full                                # Raw Mutect Ouput
        Array[File] mutect_pon_annotated_vcf = archer_workflow.mutect_pon_annotated_vcf                 # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        Array[File] mutect_vep_annotated_vcf = archer_workflow.mutect_vep_annotated_vcf                  # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Lofreq
        Array[File] lofreq_full = archer_workflow.lofreq_full                                # Raw Lofreq Ouput
        Array[File] lofreq_pon_annotated_vcf = archer_workflow.lofreq_pon_annotated_vcf                     # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        Array[File] lofreq_vep_annotated_vcf = archer_workflow.lofreq_vep_annotated_vcf                  # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Vardict
        Array[File] vardict_full = archer_workflow.vardict_full                              # Raw Vardict Ouput
        Array[File] vardict_pon_annotated_vcf = archer_workflow.vardict_pon_annotated_vcf                  # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        Array[File] vardict_vep_annotated_vcf = archer_workflow.vardict_vep_annotated_vcf               # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Pindel
        Array[File] pindel_full = archer_workflow.pindel_full                                 # Raw Pindel Ouput

        Array[File] pon_total_counts = archer_workflow.pon_total_counts                               # PoN Pileup Results
        Array[File] fpfilter_results = archer_workflow.fpfilter_results
        Array[File] vep_results = archer_workflow.vep_results

        #File gnomAD_exclude = get_gnomad_exclude.normalized_gnomad_exclude

        # R Things
        Array[File] mutect_annotate_pd = archer_workflow.mutect_annotate_pd
        Array[File] lofreq_annotate_pd = archer_workflow.lofreq_annotate_pd
        Array[File] vardict_annotate_pd = archer_workflow.vardict_annotate_pd

        # Model
        Array[File] model_output = archer_workflow.model_output
        Array[File] model_raw_output = archer_workflow.model_raw_output
        Array[File] mutect_complex = archer_workflow.mutect_complex
        Array[File] pindel_complex = archer_workflow.pindel_complex
        Array[File] lofreq_complex = archer_workflow.lofreq_complex
        Array[File] caller_filters = archer_workflow.caller_filters
    }

}

