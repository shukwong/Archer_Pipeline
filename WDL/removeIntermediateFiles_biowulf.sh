folder_names=( "call-ScatterAt805_18" "call-ScatterAt618_18" "call-ScatterAt777_18" "call-lofreqTumorTask" "call-vep" "call-interval_to_bed"   "call-mergeCallers"       "call-mutect_isec_complement_gnomAD"  "call-pindelSanitizeVcf"    "call-vardict_pon2" "call-bcbio_filter" "call-lofreq_annotate_vcf"  "call-mutectNormalize" "call-pindelToVcf"    "call-vardictSanitizeVcf" "call-lofreq_call_R_fisher"       "call-mutect_pon2" "call-pindelTumorCat"     "call-selectVariants"  "call-vardictTumorTask" "call-lofreq_isec_complement_gnomAD"  "call-mutectSanitizeVcf" "call-pindelTumorTask"    "call-split_bed_to_chr"  "call-lofreqNormalize"  "call-mutectTumorTask" "call-realign" "call-vardict_annotate_vcf"  "call-lofreq_pon2"  "call-reformat"  "call-vardict_call_R_fisher" "call-filterClipAndCollectMetrics"  "call-lofreqSanitizeVcf"  "call-mutect_annotate_vcf"   "call-pileup_merge"  "call-removeEndTags"  "call-vardict_isec_complement_gnomAD" "call-fpFilter" "call-mutect_call_R_fisher"  "call-pindelNormalize"  "call-vardictNormalize" )

for folder in "${folder_names[@]}"
do
   echo $folder 

   rm -rf call-archer*/*/boltonlab_CH/*/${folder}
done



#rm -rf shard*/boltonlab_CH/*/call-ScatterAt805_18/*
