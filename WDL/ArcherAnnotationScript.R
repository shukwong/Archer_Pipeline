#!/usr/bin/env Rscript

# Irenaeus Chan
# January 7th, 2022

# Loading Libraries
library(optparse)
library(vcfR)
library(dplyr)
library(stringr)
library(sqldf)
library(jsonlite)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input VCF File e.g. mutect.sample_name.final.annotated.vcf.gz", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output Filename e.g. mutect.sample_name.final.annotated.Rscript", metavar="character"),
  make_option(c("-c", "--caller"), type="character", default=NULL,
              help="The VCF from which caller e.g. mutect | lofreq | vardict", metavar="character"),
  make_option("--bolton_bick_vars", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/bick.bolton.vars3.txt",
              help="The directory path where the 'bick.bolton.vars' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--mut2_bick", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/topmed.n2.mutation.1.c.p.txt",
              help="The directory path where the 'bick_topmed' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--mut2_kelly", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/kelly.n2.mutation.1.c.p.txt",
              help="The directory path where the 'kelly_impact' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--matches2", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/matches.2.c.p.txt",
              help="The directory path where the 'bick_kelly' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--TSG_file", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/gene_census_TSG.txt",
              help="The directory path where the 'TSG' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--oncoKB_curated", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/all_curated_genes_v2.0.tsv",
              help="The directory path where the 'oncoKB_genes' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--pd_annotation_file", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/pd_table_kbreview_bick_trunc3.txt",
              help="The directory path where the 'pd_table' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--pan_myeloid", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/panmyeloid_variant_counts.vep.annotated.vcf.tsv",
              help="The directory path where the 'panmyeloid_variant' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--cosmic_dir", type="character", default=NULL,
              help="The directory path where the 'cosmic' files are stored for cosmic function.", metavar="character"),
  make_option("--p_value", type="double", default=2.114164905e-6,
              help="The Bonferroni Corrected P-Value [default = %default]", metavar="double"),
  make_option("--truncating", type="character", default="/Users/irenaeuschan/Documents/Irenaeus/data/BB.truncating.more.than.1.tsv",
              help="The directory path where the 'truncating' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--old", type="logical", default=FALSE,
              help="If this VCF was run before SpliceAI [default = %default]", metavar="logical")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Must supply <input_file_name>", call.=FALSE)
} else if (is.null(opt$out)) {
  print_help(opt_parser)
  stop("Must supply <output_file_name>", call.=FALSE)
} else if (is.null(opt$caller)) {
  print_help(opt_parser)
  stop("Must supply <caller>", call.=FALSE)
}

opt$caller <- tolower(opt$caller)

if (opt$caller == "mutect") {
  vcf <- vcfR2tidy(read.vcfR(opt$input), single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat
} else if (opt$caller == "lofreq") {
  vcf <- vcfR2tidy(read.vcfR(opt$input), info_only = TRUE, info_types = TRUE, format_types = TRUE)$fix
} else if (opt$caller == "vardict") {
  vcf <- vcfR2tidy(read.vcfR(opt$input), single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat
  vcf <- vcf %>% dplyr::rename(QUAL = QUAL...6, ReadQual = QUAL...20)
}

# Append FP_filter Column to FILTER Column
vcf$FP_filter <- str_replace_all(vcf$FP_filter, ",",";")
vcf <- vcf %>% dplyr::mutate(FILTER = ifelse(FP_filter=="PASS", 
                                             FILTER, 
                                             ifelse(FILTER=="PASS",
                                               FP_filter,
                                               paste0(FILTER, ";", FP_filter)
                                             )))

colnames(vcf) <- ifelse(grepl("gt_", colnames(vcf)) |
                        colnames(vcf) == "QUAL" |
                        colnames(vcf) == "FILTER",
                      paste0(colnames(vcf), "_", opt$caller) ,colnames(vcf))
vcf <- vcf %>% mutate(subject = unlist(lapply(strsplit(SAMPLE, "_", fixed = TRUE), "[[", 1)))

# Number of Samples
count_threshold<-ceiling(length(unique((vcf$subject)))*0.1)
count_threshold = 5

minVAF <- 0.005
minAltCount <- 5    # This is for the minimum alt counts filter

# Remove Off-Target & Introns          
vcf <- vcf %>% dplyr::filter(!is.na(CSQ))  

# Create a Universal key to be used to compare across datasets
vcf <- vcf %>% dplyr::mutate(key = paste0(CHROM, ':', POS, ' ', REF, '>', ALT))

# Reformating Ref and Alt Counts
if (opt$caller == "mutect") {
  vcf <- vcf %>% tidyr::separate(gt_AD_mutect, c("gt_AD_ref", "gt_AD_alt"), sep=",", extra = "merge", fill = "right")
  vcf[,"gt_AF_mutect"] <- lapply(vcf[,"gt_AF_mutect"], function(x) as.numeric(as.character(x)))
  vcf[,c("gt_AD_ref", "gt_AD_alt")] <- lapply(vcf[,c("gt_AD_ref", "gt_AD_alt")], function(x) as.numeric(as.character(x)))
  vcf <- vcf %>% dplyr::rename(gt_AD_ref_mutect = gt_AD_ref, gt_AD_alt_mutect = gt_AD_alt)
  # Reformatting StrandBias
  # SB= RefFor, RefRev, AltFor, AltRev
  vcf <- vcf %>% tidyr::separate(gt_SB_mutect, c("RDF_mutect", "RDR_mutect", "ADF_mutect", "ADR_mutect"), sep=",", extra = "merge", fill = "right")
  vcf[,c("RDF_mutect", "RDR_mutect", "ADF_mutect", "ADR_mutect")] <- lapply(vcf[,c("RDF_mutect", "RDR_mutect", "ADF_mutect", "ADR_mutect")], function(x) as.numeric(as.character(x)))
} else if (opt$caller == "lofreq") {
  vcf <- vcf %>% tidyr::separate(DP4, c("RDF_lofreq", "RDR_lofreq", "ADF_lofreq", "ADR_lofreq"), sep=",", extra = "merge", fill = "right")
  vcf[,c("RDF_lofreq", "RDR_lofreq", "ADF_lofreq", "ADR_lofreq")] <- lapply(vcf[,c("RDF_lofreq", "RDR_lofreq", "ADF_lofreq", "ADR_lofreq")], function(x) as.numeric(as.character(x)))
  vcf <- vcf %>% dplyr::mutate(gt_AD_alt_lofreq = ADF_lofreq + ADR_lofreq, gt_AD_ref_lofreq = RDF_lofreq + RDR_lofreq)
  vcf <- vcf %>% dplyr::rename(gt_AF_lofreq = AF)
} else if (opt$caller == "vardict") {
  vcf <- vcf %>% tidyr::separate(gt_AD_vardict, c("gt_AD_ref", "gt_AD_alt"), sep=",", extra = "merge", fill = "right")
  vcf[,"gt_AF_vardict"] <- lapply(vcf[,"gt_AF_vardict"], function(x) as.numeric(as.character(x)))
  vcf[,c("gt_AD_ref", "gt_AD_alt")] <- lapply(vcf[,c("gt_AD_ref", "gt_AD_alt")], function(x) as.numeric(as.character(x)))
  vcf <- vcf %>% dplyr::rename(gt_AD_ref_vardict = gt_AD_ref, gt_AD_alt_vardict = gt_AD_alt)
  # Reformatting StrandBias
  # SB= RefFor, RefRev, AltFor, AltRev
  vcf <- vcf %>% tidyr::separate(gt_RD_vardict, c("RDF_vardict", "RDR_vardict"), sep=",", extra = "merge", fill = "right")
  vcf <- vcf %>% tidyr::separate(gt_ALD_vardict, c("ADF_vardict", "ADR_vardict"), sep=",", extra = "merge", fill = "right")
  vcf[,c("RDF_vardict", "RDR_vardict", "ADF_vardict", "ADR_vardict")] <- lapply(vcf[,c("RDF_vardict", "RDR_vardict", "ADF_vardict", "ADR_vardict")], function(x) as.numeric(as.character(x)))
}

if (opt$old) {
  CSQ_string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|VAR_SYNONYMS|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas" 
} else {
  CSQ_string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|VAR_SYNONYMS|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas" 
}
CSQnames <- str_split(CSQ_string, "\\|")[[1]]
vcf <- vcf %>% tidyr::separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")

# gnomAD 
gnomAD.col <- grep("gnomAD_.*_VEP", colnames(vcf))
vcf[,gnomAD.col] <- apply(vcf[,gnomAD.col], 2, as.numeric)
vcf$max_gnomAD_AF_VEP <- as.numeric(apply(vcf[,gnomAD.col], 1, function(x){max(x, na.rm = T)}))

# gnomADe /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz
gnomADe.col <- grep("gnomADe_.*_VEP", colnames(vcf))
vcf[,gnomADe.col] <- apply(vcf[,gnomADe.col], 2, as.numeric)
vcf$max_gnomADe_AF_VEP <- as.numeric(apply(vcf[,gnomADe.col], 1, function(x){max(x, na.rm = T)}))

# gnomADg http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/r3.0/
gnomADg.col <- grep("gnomADg_.*_VEP", colnames(vcf))
vcf[,gnomADg.col] <- apply(vcf[,gnomADg.col], 2, as.numeric)
vcf$max_gnomADg_AF_VEP <- as.numeric(apply(vcf[,gnomADg.col], 1, function(x){max(x, na.rm = T)}))

vcf$max_gnomAD_AF <- (vcf$max_gnomAD_AF_VEP < 0.005 & vcf$max_gnomADe_AF_VEP < 0.005 & vcf$max_gnomADg_AF_VEP < 0.005)

vcf <- vcf %>%
  mutate(
    # Obtain the 10 Nucleotides Upstream to the Mutation Region
    context_5 = as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, CHROM, start = POS-10, end = POS-1)),
    # Obtain the 10 Nucleotides Downstream to the Mutation Region
    context_3 = as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, CHROM, start = POS + nchar(REF), end = POS + nchar(REF) + 9)),
    # Obtain the 10 Nucleotide Frame surrounding the Mutation Region
    context_10 = as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, CHROM, start = POS-5, end = POS+5))
  ) %>%
  mutate(
    # Calculate the Dust Scores for the Regions Obtained
    dust_score_5 = R453Plus1Toolbox::complexity.dust(Biostrings::DNAStringSet(context_5)),
    dust_score_3 = R453Plus1Toolbox::complexity.dust(Biostrings::DNAStringSet(context_3)),
    dust_score_10 = R453Plus1Toolbox::complexity.dust(Biostrings::DNAStringSet(context_10)),
  ) %>% rowwise() %>%
  # Determine the Dust Score and Shannon Score
  mutate (
    # If either the upstream or downstream regions have low complexity, it can potentially cause sequencing
    # artifacts. Therefore, we want to consider the lowest complexity on either side
    dust_score = max(dust_score_5, dust_score_3, dust_score_10)           # Dust > 7 is "Low Complexity"
  )

# Homopolymer Filter
sameLetters <- function (x) {
  letter <- substr(x, 1, 1)
  for (pos in 1:nchar(x)){
    if (substr(x, pos, pos) != letter) { return(FALSE) }
  }
  return(TRUE)
}

vcf <- vcf %>%
  mutate(
    # Determines if the Downstream sequence has 3 of the same nucleotides in a row
    case_NXXX = str_detect(context_3, regex('^[A]{3,}|^[G]{3,}|^[T]{3,}|^[C]{3,}')) &
      # Checks if the last character of the variant completes the homopolymer
      substr(ALT, nchar(ALT), nchar(ALT)) == substr(context_3, 1, 1),
    # Sandwich variant
    case_XNXX = ifelse(nchar(ALT) == 1, TRUE, FALSE) & # Can only occur with Single Nucleotides
      # Check if the last character of the upstream + the variant + the first two characters of downstream
      # completes the homopolymer
      paste0(substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
             ALT,
             substr(context_3, start = 1, stop = 2)) %>%
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Second Sandwich Case
    case_XXNX = ifelse(nchar(ALT) == 1, TRUE, FALSE) &
      paste0(substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
             ALT,
             substr(context_3, start = 1, stop = 1)) %>%
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Determines if the Upstream sequence ends with 3 of the same nucleotides in a row
    case_XXXN = sameLetters(substr(context_5, start = nchar(context_5)-2, stop = nchar(context_5))) &
      # Checks if the first character of the variant completes the homopolymer
      substr(ALT, 1, 1) == substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
    # Same as Case NXXX but for dinucleotides
    case_NNXX = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & str_detect(context_3, regex('^[A]{2,}|^[G]{2,}|^[T]{2,}|^[C]{2,}')) &
      substr(ALT, nchar(ALT)-1, nchar(ALT)) == substr(context_3, 1, 2),
    # Sandiwch Case
    case_XNNX = ifelse(nchar(ALT) == 2, TRUE, FALSE) & sameLetters(ALT) &
      paste0(
        substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
        ALT,
        substr(context_3, start = 1, stop = 2)) %>%
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Same as Case XXXN but for dinucleotides
    case_XXNN = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & sameLetters(substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5))) &
      substr(ALT, 1, 2) == substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
  ) %>%
  mutate(
    ## first get deletion part
    deletion = ifelse(nchar(REF)>nchar(ALT), substr(REF, 2, nchar(REF)), ""),
    deletion_plus_context_3 = ifelse(
      nchar(REF)>nchar(ALT),
      paste0(deletion, context_3),
      ""
    ),
    case_deletion_in_homopolymer = ifelse(
      nchar(REF)>nchar(ALT),
      paste0(deletion, context_3) %>%
        str_detect(regex('^[A]{4,}|^[G]{4,}|^[T]{4,}|^[C]{4,}')),
      FALSE
    )
  )

# Homopolymer Filter
sameLetters <- function (x) {
  letter <- substr(x, 1, 1)
  for (pos in 1:nchar(x)){
    if (substr(x, pos, pos) != letter) { return(FALSE) }
  }
  return(TRUE)
}

vcf <- vcf %>%
  mutate(
    # Determines if the Downstream sequence has 3 of the same nucleotides in a row
    case_NXXX = str_detect(context_3, regex('^[A]{3,}|^[G]{3,}|^[T]{3,}|^[C]{3,}')) &
      # Checks if the last character of the variant completes the homopolymer
      substr(ALT, nchar(ALT), nchar(ALT)) == substr(context_3, 1, 1),
    # Sandwich variant
    case_XNXX = ifelse(nchar(ALT) == 1, TRUE, FALSE) & # Can only occur with Single Nucleotides
      # Check if the last character of the upstream + the variant + the first two characters of downstream
      # completes the homopolymer
      paste0(substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
             ALT,
             substr(context_3, start = 1, stop = 2)) %>%
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Second Sandwich Case
    case_XXNX = ifelse(nchar(ALT) == 1, TRUE, FALSE) &
      paste0(substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
             ALT,
             substr(context_3, start = 1, stop = 1)) %>%
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Determines if the Upstream sequence ends with 3 of the same nucleotides in a row
    case_XXXN = sameLetters(substr(context_5, start = nchar(context_5)-2, stop = nchar(context_5))) &
      # Checks if the first character of the variant completes the homopolymer
      substr(ALT, 1, 1) == substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
    # Same as Case NXXX but for dinucleotides
    case_NNXX = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & str_detect(context_3, regex('^[A]{2,}|^[G]{2,}|^[T]{2,}|^[C]{2,}')) &
      substr(ALT, nchar(ALT)-1, nchar(ALT)) == substr(context_3, 1, 2),
    # Sandiwch Case
    case_XNNX = ifelse(nchar(ALT) == 2, TRUE, FALSE) & sameLetters(ALT) &
      paste0(
        substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
        ALT,
        substr(context_3, start = 1, stop = 2)) %>%
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Same as Case XXXN but for dinucleotides
    case_XXNN = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & sameLetters(substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5))) &
      substr(ALT, 1, 2) == substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
  )

vcf <- vcf %>%
  mutate(
    pass_homopolymer_filter =
      !(case_NXXX) &
      !(case_XNXX) &
      !(case_XXNX) &
      !(case_XXXN) &
      !(case_NNXX) &
      !(case_XNNX) &
      !(case_XXNN)
  )

# Strand Bias and Alt Alleles Filter
if (opt$caller == "mutect") {
  vcf<-vcf %>% mutate(pass_strand_bias_mutect = ifelse(ADF_mutect/(ADF_mutect+ADR_mutect) < 0.1 | ADF_mutect/(ADF_mutect + ADR_mutect) > 0.9, 0, ifelse((ADF_mutect+ADR_mutect) < 5, 0, 1)))
} else if (opt$caller == "lofreq") {
  vcf<-vcf %>% mutate(pass_strand_bias_lofreq = ifelse(ADF_lofreq/(ADF_lofreq+ADR_lofreq) < 0.1 | ADF_lofreq/(ADF_lofreq + ADR_lofreq) > 0.9, 0, ifelse((ADF_lofreq+ADR_lofreq) < 5, 0, 1)))
} else if (opt$caller == "vardict") {
  vcf<-vcf %>% mutate(pass_strand_bias_vardict = ifelse(ADF_vardict/(ADF_vardict+ADR_vardict) < 0.1 | ADF_vardict/(ADF_vardict + ADR_vardict) > 0.9, 0, ifelse((ADF_vardict+ADR_vardict) < 5, 0, 1)))
}

# Filter PoN and Save Removed to N_samples File
#n_samples <- vcf %>% filter(PON_FISHER > opt$p_value)
if (opt$caller == "mutect") {
  n_samples <- vcf %>% select(key, SAMPLE, gt_AD_alt_mutect, gt_AF_mutect)  
} else if (opt$caller == "lofreq") {
  n_samples <- vcf %>% select(key, SAMPLE, gt_AD_alt_lofreq, gt_AF_lofreq)
} else if (opt$caller == "vardict") {
  n_samples <- vcf %>% select(key, SAMPLE, gt_AD_alt_vardict, gt_AF_vardict)
}
#write.table(n_samples, paste0(opt$out, ".n_samples.tsv"), row.names = FALSE, sep="\t")

AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               '%3D'='=', '='='=')

# Annotation Preparation
vcf$AAchange <- gsub("(.*p\\.)(.*)", "\\2", vcf$HGVSp_VEP)
for (i in 1:length(AminoAcids)) { vcf$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], vcf$AAchange) }
vcf$gene_loci_p <- paste(vcf$SYMBOL_VEP,
                       paste0(sapply(vcf$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1]),
                              as.numeric(str_extract(vcf$AAchange, "\\d+"))),sep = "_")
vcf$gene_loci_c <- paste(vcf$SYMBOL_VEP, gsub(".*:", "", vcf$HGVSc_VEP), sep = "_")
vcf$gene_loci_vep <- ifelse(is.na(vcf$gene_loci_p), vcf$gene_loci_c, vcf$gene_loci_p)
vcf$gene_aachange <- with(vcf, paste(SYMBOL_VEP, AAchange, sep = "_"))
vcf$gene_cDNAchange <- paste(vcf$SYMBOL_VEP, gsub(".*:","",vcf$HGVSc_VEP), sep="_")

vars <- read.table(opt$bolton_bick_vars, sep = "\t", header = T, comment.char = "")
vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")
vars <- vars %>% dplyr::mutate(key = paste0(CHROM, ':', POS, ' ', REF, '>', ALT))

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

dims <- dim(vcf)[[1]]
vcf <- sqldf("SELECT m.*, r.`n.loci.vep`, r.`source.totals.loci`
            FROM `vcf` as m
            LEFT JOIN `vars` as r
            on m.key = r.key OR m.gene_loci_vep = r.gene_loci_vep")
vcf <- vcf[!duplicated(vcf),]
## make sure aachange exists as in doesn't end with an '_'; example: DNMT3A_ for splice 
vcf <- sqldf("SELECT m.*, r.`n.HGVSp`, r.`source.totals.p`
            FROM `vcf` as m
            LEFT JOIN `vars` as r
            on m.key = r.key OR (m.gene_aachange = r.gene_aachange) AND r.gene_aachange NOT LIKE '%_'")
vcf <- vcf[!duplicated(vcf),]
vcf <- sqldf("SELECT m.*, r.`n.HGVSc`, r.`source.totals.c`
            FROM `vcf` as m
            LEFT JOIN `vars` as r
            on m.key = r.key OR (m.gene_cDNAchange = r.gene_cDNAchange) AND r.gene_cDNAchange NOT LIKE '%_'")
vcf <- vcf[!duplicated(vcf),]
paste0("dims match after sqldf: ",dim(vcf)[[1]] == dims)

vcf <- vcf %>% arrange(CHROM, POS, REF, ALT)

# Annotate_PD
# -------------

# mut2-bick, mut2-kelly
# A list of mutations (HGVSp and HGVSc) that appear in Bick or Kelly dataset more than twice
# Matches is essentially bick + bolton >= 2, not just individual dataset (used for ch_pd2)
topmed.mutation.2 <- read.table(opt$mut2_bick, sep = "\t", header = T, comment.char = "")
kelly.mutation.2 <- read.table(opt$mut2_kelly, sep = "\t", header = T, comment.char = "")
matches.2.c.p <- read.table(opt$matches2, sep = "\t", header = T, comment.char = "")

#'@returns dataframe with equal number of rows as input with binded counts columns
#'@param df dataframe: this is the MUTS dataframe; needs to have `HGVSp_VEP`,`var_key`,`sample` columns
#'@param sample.column string: column name for sample name/id
#' MUTS MUST BE ORDERED by CHROM (and maybe POS? just to be safe?? Basically whatever column is paralleled needs to be sorted by...)
cosmic.run <- function(df, sample.column, var_key.column) {
  # COSV52681947
  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  X.1 <- parallel::mclapply(ls(X1), function(x) {
    # X.1 <- lapply(ls(X1), function(x) {
    print(x)
    df.chr <- X1[[x]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP)
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    # df.chr$skey <- with(df.chr, paste(var_key, eid, sep = ":"))
    df.chr$skey <- paste(df.chr[[var_key.column]], df.chr[[sample.column]], sep = ":")
    
    cosmic <- read.table(paste0(opt$cosmic_dir, "/CosmicMutantExport.final.more.minimal.test.",x,".tsv"),
                         # cosmic <- read.table(paste0("~/test/CosmicMutantExport.final.minimal.",x,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
    colnames(cosmic) <- c("COSMIC_ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")
    # colnames(df.chr)[grep("key", colnames(df.chr))] <- "var_key"
    cosmic.test <- sqldf("SELECT l.*, r.COSMIC_ID, r.var_key as var_key_cosmic, HGVSP, r.CosmicCount, 
                    r.heme_cosmic_count, r.myeloid_cosmic_count
                    FROM `df.chr` as l
                    LEFT JOIN `cosmic` as r
                    on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
    cosmic.test <- cosmic.test[,c(colnames(df.chr),"COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count")]
    cosmic.test$CosmicCount <- fillna(cosmic.test$CosmicCount,0)
    cosmic.test$heme_cosmic_count <- fillna(cosmic.test$heme_cosmic_count,0)
    cosmic.test$myeloid_cosmic_count <- fillna(cosmic.test$myeloid_cosmic_count,0)
    cosmic.test <- cosmic.test %>% 
      group_by(var_key,cosmic.test[[sample.column]]) %>%
      slice_max(order_by = heme_cosmic_count, n = 1, with_ties = F)
    cosmic.test <- data.frame(cosmic.test)
    
    #cosmic.test$skey <- with(cosmic.test, paste(var_key, eid, sep = ":"))
    cosmic.test$skey <- paste(cosmic.test[[var_key.column]], cosmic.test[[sample.column]], sep = ":")
    row.names(cosmic.test) <- cosmic.test$skey
    # row.names(cosmic.test) <- cosmic.test$var_key
    ## put back into row order of df.chr
    cosmic.test <- cosmic.test[df.chr$skey,]
    # cosmic.test <- cosmic.test[df.chr$var_key,]
    
    if (all(cosmic.test[,c("CHROM","POS","REF","ALT",sample.column)]==df.chr[,c("CHROM","POS","REF","ALT",sample.column)])) {
      print(paste("good inner", x))
      return(cosmic.test)
    } else {
      print(paste("bad inner", x))
    }
  })
  proc.time() - ptm
  
  if (length(X1) == length(X.1)) { 
    names(X.1) <- ls(X1) 
  } else {
    stop("WRONG dimensions of df split by CHROM")
  }
  ## are all the names in the same order
  all.names <- all(names(X1)==names(X.1))
  all.match <- all(sapply(ls(X1), function(x) {
    dim(X1[[x]])[1] == dim(X.1[[x]])[1]
  }))
  if (all.match & all.names) {
    cosmic.final2 <- bind_rows(X.1)
  }
  if (all(cosmic.final2[,c("CHROM","POS","REF","ALT",sample.column)]==df[,c("CHROM","POS","REF","ALT",sample.column)])) {
    print("COSMIC good #################################################")
    MUTS <- cbind(df, cosmic.final2 %>% dplyr::select(COSMIC_ID,CosmicCount,heme_cosmic_count,myeloid_cosmic_count))
    return(MUTS)
  } else {
    print("COSMIC bad #################################################")
  }
}

annotate.PD <- function(x) {
  MUTS <- x
  
  ## has 9
  ## missing 6
  # additional we can have "upstream_gene_variant", "non_coding_transcript_variant", "non_coding_transcript_exon_variant"
  # total 18 unique(unlist(sapply(MUTS$Consequence_VEP, function(x) { str_split(x,"&",simplify = TRUE)})))
  translate_consequence = c(
    "frameshift_variant" = "frameshift_variant", #1
    
    "inframe_deletion" = "inframe_deletion", #2
    "inframe_insertion" = "inframe_insertion", #3
    "synonymous_variant" = "synonymous_variant", #4
    "missense_variant" = "missense_variant", #5
    "intron_variant" = "intron_variant", #6
    
    "splice_region_variant" = "splice_region_variant", #7
    "splice_donor_variant" = "splice_region_variant", # missing
    "splice_acceptor_variant" = "splice_region_variant", # missing
    
    # "coding_sequence_variant" = "feature_truncation", # manually assigned based on length ref & alt
    # "protein_altering_variant" = "feature_truncation", # manually assigned based on length ref & alt
    "feature_truncation" = "feature_truncation",
    
    "3_prime_UTR_variant" = "3_prime_UTR_variant", # missing
    "5_prime_UTR_variant" = "5_prime_UTR_variant", # missing
    
    
    "stop_gained" = "stop_gained", #8
    "start_lost" = "start_lost", #9
    "stop_lost" = "missense_variant", # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
    "stop_retained_variant" = "synonymous_variant" # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
  )
  
  # translate_type = c(
  #   'SNV' = 'SNV',
  #   'deletion' = 'Del',
  #   'insertion' = 'Ins')
  # MUTS <- MUTS %>%
  #   mutate(VariantClass = translate_consequence[str_split(Consequence_VEP,"&",simplify = TRUE)[,1]]) %>%
  #   mutate(VariantClass = ifelse(VariantClass == 'Frame_Shift_', paste0('Frame_Shift_', translate_type[VARIANT_CLASS_VEP]), VariantClass))
  
  MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP",
                  "HGVSp_VEP","n.HGVSp","HGVSc_VEP","n.HGVSc",
                  "Consequence_VEP","SAMPLE","EXON_VEP","AAchange")]
  
  MUTS <- MUTS %>%
    mutate(VariantClass = str_split(Consequence_VEP,"&",simplify = TRUE)[,1]) %>%
    mutate(VariantClass = case_when((VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) > 0 ~ "inframe_deletion",
                                    (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) < 0 ~ "inframe_insertion",
                                    (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 != 0 ~ "frameshift_variant",
                                    (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF) == nchar(ALT) ~ "missense_variant",
                                    TRUE ~ VariantClass)) %>%
    mutate(VariantClass = translate_consequence[VariantClass])
  
  
  MUTS$aa_ref <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1])
  MUTS$aa_alt <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][2])
  MUTS$aa_pos <- as.numeric(str_extract(MUTS$AAchange, "\\d+"))
  MUTS$var_key = paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  
  
  
  ###create list of tumor suppressor genes####
  ## file 1
  gene_census = read.table(opt$TSG_file, comment.char = "", sep = "\t", quote = "", header = T)
  ## file 2
  oncoKB_curated = read.table(opt$oncoKB_curated, comment.char = "", sep = "\t", quote = "", header = T)
  oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
  TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
  
  
  ###########annotate with COSMIC########
  #source("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/cosmic.parallel.compute1.R")
  ## need to add a parallel
  ## MUTS MUST BE ORDERED by CHROM (and maybe POS? just to be safe?? Basically whatever column is paralleled needs to be sorted by...)
  message("annotating variants with COSMIC...")
  MUTS <- cosmic.run(MUTS, "SAMPLE", "var_key")
  MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  message("COSMIC done\n")
  ##
  
  #annotate oncokb
  #get_oncokb = function(mut) {
  #  
    # apiKey = "a83627cd-47e4-4be0-82dd-8f4cc9e4d6d0"
  #  apiKey = "da03f096-284a-490a-8b3f-95c5f1bf5666"
  #  request_url = paste0("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
  #                       mut['SYMBOL_VEP'], "&alteration=", mut['AAchange'], "&consequence=", mut['VariantClass'])
  #  res=httr::content(httr::GET(request_url, httr::add_headers(Authorization = paste("Bearer", apiKey))))
  #  return(res$oncogenic)
  #}
  
  message("annotating variants with oncoKB...")
  #h = curl::new_handle()
  #curl::handle_setopt(h, http_version = 2)
  #httr::set_config(httr::config(http_version = 0))
  #cl = parallel::makeCluster(4)
  #MUTS$oncoKB = parallel::parApply(cl, MUTS, 1, get_oncokb)
  
  json_query <- apply(MUTS, 1, function(x) {
    my_list <- list(
      alteration=x["AAchange"],
      consequence=x["VariantClass"],
      gene=list(
        hugoSymbol=x["SYMBOL_VEP"]
      )
    )
    return(my_list)
  })
    
  MUTS$oncoKB<-sapply(httr::content(httr::POST(url="https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange", 
                                                     httr::add_headers(c('Content-Type'='application/json', 'Authorization'='Bearer da03f096-284a-490a-8b3f-95c5f1bf5666')), 
                                                     body = jsonlite::toJSON(unname(json_query), auto_unbox = TRUE))), "[[", "oncogenic")
  #parallel::stopCluster(cl)
  message("oncoKB done\n")
  
  ##annotate myeloid-PD
  
  # Variant Annotation
  # Variants were annotated according to evidence for functional relevance in cancer (putative driver or CH-PD)
  # and for relevance to myeloid neoplasms specifically (CH-Myeloid-PD).
  # We annotated variants as oncogenic in myeloid disease (CH-Myeloid-PD) if they fulfilled any of the following criteria:
  #
  #   1. Truncating variants in NF1, DNMT3A, TET2, IKZF1, RAD21, WT1, KMT2D, SH2B3, TP53, CEBPA, ASXL1, RUNX1, BCOR, KDM6A, STAG2, PHF6, KMT2C, PPM1D, ATM, ARID1A, ARID2, ASXL2, CHEK2, CREBBP, ETV6, EZH2, FBXW7, MGA, MPL, RB1,SETD2, SUZ12, ZRSR2 or in CALR exon 9
  #   2. Translation start site mutations in SH2B3
  #   3. TERT promoter mutations
  #   4. FLT3-ITDs
  #   5. In-frame indels in CALR, CEBPA, CHEK2, ETV6, EZH2
  
  # MUTS2 = MUTS
  
  ##read in PD annoation files from papaemmanuil lab##
  ## file 4
  A <- read.table(opt$pd_annotation_file, sep = "\t", header = T, quote = "")[,1:9]
  ## about ~68 genes
  # ch.my.pd
  ch.my.pd.genes <- unique(A$Gene[A$ch_my_pd==1])
  
  A$keys = apply(A, 1, function(x) {
    names(x)[x != '*'] %>% .[!. %in% c('source', 'ch_my_pd', 'ch_pancan_pd')]
  })
  pd_dfs = split(A, paste(A$keys))
  
  # MUTS <- MUTS2
  
  MUTS$aa_pos <- as.character(MUTS$aa_pos)
  MUTS$Exon <- sapply(MUTS$EXON_VEP, function(x) str_split(x,"/")[[1]][1])
  colnames(MUTS)[colnames(MUTS)=="SYMBOL_VEP"] <- c("Gene")
  MUTS.temp <- MUTS
  matched = list()
  ch_my_variants = NULL
  
  message("annotating variants with PD Table...")
  
  ## anything truncating in PD_table
  truncating <- c("frameshift_variant", "splice_region_variant", "stop_gained")
  for (i in 1:length(pd_dfs)){
    print(i)
    temp = as.data.frame(pd_dfs[i])
    colnames(temp) = colnames(A)
    keys = unlist(unique(temp$keys))
    if ('Exon' %in% keys) {
      MUTS.temp <- MUTS.temp %>% filter(!(var_key %in% matched))
      temp <- temp %>% filter(Exon != '*')
    } else {
      MUTS.temp <- MUTS.temp %>% filter(!(var_key %in% matched))
      temp <- temp %>% filter(Exon == '*')
    }
    # are all keys is both tables?
    print(all(keys %in% colnames(temp), keys %in% colnames(MUTS.temp)))
    pds = left_join(MUTS.temp, temp[,c(keys,'ch_my_pd','ch_pancan_pd','source')], by=keys)
    # pds = pds %>% filter(ch_my_pd>0)
    pds_ch_my_pd = pds %>% filter(ch_my_pd > 0)
    pds_ch_pd = pds %>% filter(ch_pancan_pd >=1) # & VariantClass %in% truncating
    # matched = append(matched, pds$var_key)
    matched = append(matched, pds_ch_my_pd$var_key)
    matched = append(matched, pds_ch_pd$var_key)
    # ch_my_variants = rbind(ch_my_variants, pds)
    ch_my_variants = rbind(ch_my_variants, pds_ch_my_pd, pds_ch_pd)
  }
  
  # just so we have the gene for reviewing
  ch_my_variants = ch_my_variants %>% 
    dplyr::rename(gene = Gene) %>%
    dplyr::select(source, var_key, ch_my_pd, ch_pancan_pd, gene) %>% unique()
  
  # MUTS.temp <- MUTS
  # matched = list()
  # ch_my_variants2 = NULL
  # for (i in 1:length(pd_dfs)){
  #   print(i)
  #   temp = as.data.frame(pd_dfs[i])
  #   colnames(temp) = colnames(A)
  #   keys = unlist(unique(temp$keys))
  #   if ('Exon' %in% keys) {
  #     MUTS.temp <- MUTS.temp %>% filter(!(var_key %in% matched))
  #     temp <- temp %>% filter(Exon != '*')
  #   } else {
  #     MUTS.temp <- MUTS.temp %>% filter(!(var_key %in% matched))
  #     temp <- temp %>% filter(Exon == '*')
  #   }
  #   # are all keys is both tables?
  #   print(all(keys %in% colnames(temp), keys %in% colnames(MUTS.temp)))
  #   pds = left_join(MUTS.temp, temp[,c(keys,'ch_my_pd','ch_pancan_pd','source')], by=keys)
  #   pds = pds %>% filter(ch_my_pd>0)
  #   matched = append(matched, pds$var_key)
  #   ch_my_variants2 = rbind(ch_my_variants2, pds)
  # }
  # length(intersect(ch_my_variants$var_key, ch_my_variants2$var_key))
  # length(setdiff(ch_my_variants$var_key, ch_my_variants2$var_key))
  # 
  # length(setdiff(ch_my_variants2$var_key, ch_my_variants2$var_key))
  # ch_my_variants_test <- ch_my_variants[! ch_my_variants$var_key %in% ch_my_variants2$var_key, ]
  # ch_my_variants_testa <- ch_my_variants[ch_my_variants$var_key %in% ch_my_variants2$var_key, ]
  # ch_my_variants_test2 <- ch_my_variants2[! ch_my_variants2$var_key %in% ch_my_variants$var_key, ]
  # ch_my_variants_test2a <- ch_my_variants2[ ch_my_variants2$var_key %in% ch_my_variants$var_key, ]
  
  # ch_my_variants2 = ch_my_variants2 %>% 
  #   dplyr::rename(gene = Gene) %>%
  #   dplyr::select(source, var_key, ch_my_pd, ch_pancan_pd, gene) %>% unique()
  
  
  
  tryCatch(
    expr = {
      ch_my_variants <- aggregate(source ~ var_key + ch_my_pd + ch_pancan_pd, data = ch_my_variants, FUN = paste, collapse = ",")
    },
    error = function(e){ 
      message(e)
    }
  )
  message("PD Table done\n")
  
  ## easier to just not return anything from TC and call this below it
  MUTS = left_join(MUTS, ch_my_variants, by="var_key")
  ## for the ifelse statement below, ch_my_pd CANNOT be NA!!!!!!!!!!!
  MUTS$ch_my_pd <- fillna(MUTS$ch_my_pd, 0)
  ## update ch_my_pd for those ~71 genes listed in Slack if B/B hotspot (12/22/2021 9:41 PM)
  MUTS$n.HGVSp <- fillna(MUTS$n.HGVSp, 0) 
  MUTS$n.HGVSc <- fillna(MUTS$n.HGVSc, 0) 
  MUTS$ch_my_pd <- ifelse((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes, 1, MUTS$ch_my_pd)  
  ## initial ch_pd just get pancan and truncating, myeloid stuff comes below
  MUTS$ch_pd <- ifelse(MUTS$VariantClass %in% truncating & MUTS$ch_pancan_pd>=1, 1, 0) 
  
  MUTS <- MUTS %>% 
    mutate(WHY_CH_ch_my_pd = ifelse(MUTS$ch_my_pd>0, "PD_table;", ""))
  
  
  # 6. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, paste0(WHY_CH_ch_my_pd, " Cosmic_heme;"), WHY_CH_ch_my_pd))
  # MUTS$ch_my_pd = ifelse(MUTS$heme_cosmic_count>=10,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  
  # 7. Any variant noted as potentially oncogenic in an in-house dataset of 7,000 individuals with myeloid neoplasm greater than or equal to 5 times
  # panmyeloid_variant_counts = read.table(opt$pan_myeloid, comment.char = "", sep = "\t", quote = "", header = T)
  # panmyeloid_variant_counts$var_key = with(panmyeloid_variant_counts, paste(CHROM,POS,REF,ALT, sep = ":"))
  # panmyeloid_variant_counts <- panmyeloid_variant_counts[panmyeloid_variant_counts$Annotation == "ONCOGENIC" | 
  #                                                          panmyeloid_variant_counts$Annotation == "SNP",]
  # 
  # MUTS$Gene_HGVSp_VEP <- paste(MUTS$Gene, gsub(".*p.","", MUTS$HGVSp_VEP), sep = "_")
  # panmyeloid_variant_counts$Gene_HGVSp_VEP <- paste(panmyeloid_variant_counts$SYMBOL_VEP, gsub(".*p.","", panmyeloid_variant_counts$HGVSp_VEP), sep = "_")
  
  # message("annotating variants with panmyeloid...")
  # MUTS <- sqldf("SELECT l.*, r.MDS, r.AML, r.MPN, r.Annotation
  #           FROM `MUTS` as l
  #           LEFT JOIN `panmyeloid_variant_counts` as r
  #           on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
  # MUTS <- MUTS %>%
  #   mutate(
  #     MDS = fillna(MDS, 0),
  #     AML = fillna(AML, 0),
  #     MPN = fillna(MPN, 0),
  #     Annotation = fillna(Annotation, ''),
  #     n_panmyeloid = ifelse(Annotation == 'ONCOGENIC', MDS + AML + MPN, 0))
  # message("panmyeloid done\n")
  # 
  # MUTS <- MUTS %>%
  #   mutate(ch_my_pd = ifelse(n_panmyeloid>=5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
  #          WHY_CH_ch_my_pd = ifelse(MUTS$n_panmyeloid>=5, paste0(WHY_CH_ch_my_pd, " Pan_Myeloid;"), WHY_CH_ch_my_pd))
  # MUTS$ch_my_pd = ifelse(MUTS$n_panmyeloid>=5,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", 1, ch_pd),
           WHY_CH_ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", "OncoKB;", ""))
  
  MUTS$isOncogenic <- MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic" | MUTS$oncoKB=="Predicted Oncogenic"
  MUTS$isTSG <- MUTS$Gene %in% TSG
  
  # 2. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  # Genes not listed in the cancer census or OncoKB were reviewed in the literature to determine if they were potentially tumor suppressor genes.# Annotate PD based on prevelance in cancer databases (COSMIC, Oncokb)
  
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(isTSG & VariantClass %in% truncating, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(isTSG & VariantClass %in% truncating, 
                                 paste0(WHY_CH_ch_pd, " TSG & Truncating;"), WHY_CH_ch_pd))
  
  ## is gene a truncating hotpot gene and is the variant in truncating class
  vars.truncating <- read.table(opt$truncating, sep = '\t', header = T)
  
  truncating_genes <- unique(vars.truncating[vars.truncating$n>=10,"SYMBOL_VEP"])
  MUTS$isTruncatingHotSpot <- ifelse(MUTS$Gene %in% truncating_genes & MUTS$VariantClass %in% truncating, 1, 0)
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(isTruncatingHotSpot, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(isTruncatingHotSpot, paste0(WHY_CH_ch_pd, " BB truncating-HS;"), WHY_CH_ch_pd),
           ch_my_pd = ifelse(isTruncatingHotSpot & Gene %in% ch.my.pd.genes, 1, ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(isTruncatingHotSpot & Gene %in% ch.my.pd.genes, paste0(WHY_CH_ch_my_pd, " BB truncating-HS;"), WHY_CH_ch_my_pd))
  
  #3. Any variant reported as somatic at least 20 times in COSMIC
  MUTS$CosmicCount <- as.numeric(MUTS$CosmicCount)
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(CosmicCount>=20, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(CosmicCount>=20, paste0(WHY_CH_ch_pd, " Cosmic;"), WHY_CH_ch_pd))
  
  #4. Any variant meeting criteria for CH-Myeloid-PD as above.
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, paste0(WHY_CH_ch_pd, " ch_my_pd>=1"), WHY_CH_ch_pd))
  
  annotate_PD.R = MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  
  # sample.x=grep("SAMPLE",colnames(x))
  # sample.annorate=grep("SAMPLE",colnames(annotate_PD.R))  # sample.annorate=grep("SAMPLE",colnames(annotate_PD.R))
  
  if(all(annotate_PD.R[,c("CHROM","POS","REF","ALT","SAMPLE")]==x[,c("CHROM","POS","REF","ALT","SAMPLE")])) {
    message("good same dim as orig")
    annotate_PD.R <- annotate_PD.R %>% # annotate_PD.R is the MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
      dplyr::mutate(gene_loci_p = paste(Gene, paste0(aa_ref, aa_pos), sep = "_"),
                    cDNAchange = gsub(".*:","", HGVSc_VEP),
                    gene_loci_c = paste(Gene, cDNAchange, sep = "_"),
                    gene_aachange = paste(Gene, AAchange, sep = "_"),
                    gene_loci_vep = ifelse(is.na(AAchange) | AAchange == "", gene_loci_c, gene_loci_p),
                    ch_pd2 = case_when(ch_pd==1 ~ 1,
                                       (gene_loci_vep %in% vars$gene_loci_vep | 
                                          gene_aachange %in% matches.2.c.p$exact_match |
                                          gene_loci_c %in% matches.2.c.p$exact_match) ~ 1,
                                       TRUE ~ 0),
                    WHY_CH_ch_pd2 = ifelse(ch_pd==1, "CH_pd=1;", ""),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_vep %in% vars$gene_loci_vep, paste0(WHY_CH_ch_pd2, " Gene_loci>=5;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_aachange %in% matches.2.c.p$exact_match, paste0(WHY_CH_ch_pd2, " gene_aachange>=2 BB;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_c %in% matches.2.c.p$exact_match, paste0(WHY_CH_ch_pd2, " gene_loci_c>=2 BB;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_aachange %in% topmed.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_aachange>=2 Bick;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_c %in% topmed.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_loci_c>=2 Bick;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_aachange %in% kelly.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_aachange>=2 Bolton;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_c %in% kelly.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_loci_c>=2 Bolton;"), WHY_CH_ch_pd2))
    
    
    annotate_PD.R$WHY_CH <- apply(annotate_PD.R, 1, function(x) {
      return(toJSON(list('ch_my_pd'=x["WHY_CH_ch_my_pd"],
                         'ch_pd2'=x["WHY_CH_ch_pd2"],
                         'ch_pd'=x["WHY_CH_ch_pd"])))
    })
    return(
      annotate_PD.R %>% # 
        dplyr::select("CHROM","POS","REF","ALT","SAMPLE",
                      "n.HGVSp","n.HGVSc",
                      "COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count",
                      "oncoKB","isOncogenic","isTSG", "isTruncatingHotSpot",
                      "ch_my_pd","ch_pd","ch_pd2",
                      "VariantClass","AAchange","Gene","WHY_CH") 
    )
  } else {
    message("Annotate PD failed: something wrong with dims")
  }
}

if (nrow(vcf) != 0) {
  annotation.vcf <- annotate.PD(vcf) 
} else {
  annotation.vcf <- data.frame(CHROM=character(), POS=integer(), REF=character(), ALT=character(), SAMPLE=character(),
                                n.HGVSp=integer(), n.HGVSc=integer(),
                                COSMIC_ID=character(), CosmicCount=integer(), heme_cosmic_count=integer(), myeloid_cosmic_count=integer(),
                                oncoKB=character(), isOncogenic=logical(), isTSG=logical(), isTruncatingHotSpot=logical(),
                                ch_my_pd=integer(), ch_pd=integer(), ch_pd2=integer(),
                                VariantClass=character(), AAchange=character(), Gene=character(), WHY_CH=character())
}

vcf <-  left_join(vcf, annotation.vcf, by=c("CHROM","POS","REF","ALT","SAMPLE"))
colnames(vcf)[(which(sapply(vcf, class)=="list"))]

vars$aa.pos <- as.numeric(str_extract(vars$loci.vep, "\\d+"))
vars$CHROM.POS <- with(vars, paste0(CHROM,":",POS))
vars$GENE.AA.POS <- with(vars, paste(SYMBOL_VEP,aa.pos,sep=":"))
vars$gene_aachange <- with(vars, paste(SYMBOL_VEP,AAchange2,sep=":"))
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")

## new lines to handle DNMT3A:NA, ASXL1:NA etc.
vars$GENE.AA.POS <- with(vars, ifelse(is.na(vars$aa.pos),NA,paste(SYMBOL_VEP,aa.pos,sep=":"))) # new line
vars$GENE.AA.POS[vars$GENE.AA.POS=="NA:NA"] = NA # new line
vars$gene_cDNAchange <-gsub("_",":", vars$gene_cDNAchange)
vars$gene_aachange <- with(vars, paste(SYMBOL_VEP,AAchange2,sep=":"))

vcf$aa.pos <- as.numeric(str_extract(vcf$AAchange.x, "\\d+"))
vcf$gene_aachange <- with(vcf, paste(SYMBOL_VEP, AAchange.x,sep=":"))

vcf$near.BB.HS <- apply(vcf[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange")], 1, function(x) {
  p = c(-3:0,1:3)
  n = c(-9:0,1:9)
  
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% vars$GENE.AA.POS
  
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% vars$CHROM.POS
  
  if (any(any.in.p)) {
    res <- unique(vars[vars$GENE.AA.POS %in% vector.p, c("gene_aachange","source.totals.p")])
    my_pre_return <- c(x[["gene_aachange"]], paste0(res$gene_aachange, "(", res$source.totals.p,")"))
    if(grepl("Ter", x["gene_aachange"])) {
      return(my_pre_return[grepl("Ter", my_pre_return)])
    } else {
      return(my_pre_return[!grepl("Ter", my_pre_return)])
    }
  } else if (any(any.in.n)) {
    res <- unique(vars[vars$CHROM.POS %in% vector.n, c("CHROM.POS", "source.totals.c", "gene_cDNAchange")])
    my_pre_return <- c(x[["gene_cDNAchange"]], paste0(res$gene_cDNAchange, "(", res$source.totals.c,")"))
    if(grepl("del|ins|dup", x["gene_cDNAchange"])) {
      return(my_pre_return[grepl("del|ins|dup", my_pre_return)])
    } else {
      return(my_pre_return[!grepl("del|ins|dup", my_pre_return)])
    }
  } else {
    return("")
  }
})

# toJSON
vcf$near.BB.HS <- sapply(vcf$near.BB.HS, function(x) {
  if (x == "") {
    return(x)
  } else {
    return(toJSON(x))
  }
})
vcf$near.BB.HS <- ifelse(grepl(",", vcf$near.BB.HS), vcf$near.BB.HS, "")

write.table(vcf, paste0(opt$out, ".tsv"), row.names = FALSE, sep="\t")
message(paste0("Finished: Writing ", opt$caller, " TSV file: ", opt$out, ".tsv"))
