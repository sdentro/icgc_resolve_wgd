library(readr)
library(GenomicRanges)
source("~/repo/icgc_resolve_wgd/functions.R")

#############################################################
# Init
#############################################################
# setwd("/Users/sd11/Documents/Projects/icgc/consensus_SNVs/assessment_consCNA_through_consSNV_6")
setwd("/Users/sd11/Documents/Projects/icgc/consensus_clonal_copynumber/201609/")

output_file = "icgc_copynumber_wgd_adjustments.txt"
purities_file = "interim_consensus/icgc_consensus_dp_master_file.txt" # file in summary table format
corrections_file = "interim_consensus/correction_summary.20160923.txt"
hum_genome_mb = 3235

# Jeffs interim consensus files
interim_diploid_prefix = "interim_consensus/consensus.best_correction_combo_thresh3/"
interim_tetraploid_prefix = "interim_consensus/consensus.discard_diploid/"
segments_file_postfix = "_segments.txt" # files in calibration format

# DPClust files
dpclust_prefix = "dpclust_clust_locs/"
dpclust_postfix = "_subclonal_structure.txt.gz" # files in calibration format

# purities = load_purites_adjustments(purities_file)
applied_corrections = parseAppliedCorrections(corrections_file)

# Create the output vector with a statement what should be done with the tumour
sample_labels = rep(NA, nrow(applied_corrections))

# Thresholds for various metrics
CONSENSUS_PROPORTION = 0.3 # Fraction of the genome agreed upon between methods
LARGEST_HOM_DEL_ALLOWED_METRIC = 10 # in Mb
CLUSTER_0.5CCF_CLONAL_ALLOWED_SITANCE = 0.02 # in frcation of CCF
LARGEST_50_50_SUBCLONAL_SEGMENT_ALLOWED = 20 # in Mb

#############################################################
# Metric 1 - complete agreement
#############################################################
# Metric 1. All samples for which there has not been any correction we accept
not_corrected = apply(applied_corrections[,2:6], 1, function(x) { all(x[!is.na(x)] == "None") })
sample_labels = setSampleLabelMultiple(not_corrected, "no_correction", sample_labels, applied_corrections)

#############################################################
# Parsing
#############################################################
# Load the breakpoint counts - no longer needed after standardisation of the breakpoints
# breakpoints = getPreCalcNrBreakpoints()
# breakpoints.m = melt(breakpoints, id="samplename")
# Get purities and power for each tumour
purities = readr::read_tsv(purities_file)

# Calculate all the metrics
all_metrics = data.frame()
for (samplename in applied_corrections$samplename[is.na(sample_labels)]) {
  
  diploid_cnprofile = paste0(interim_diploid_prefix, samplename, segments_file_postfix)
  tetraploid_cnprofile = paste0(interim_tetraploid_prefix, samplename, segments_file_postfix)
  
  metrics = calc_metrics(diploid_cnprofile, tetraploid_cnprofile, samplename, paste0(dpclust_prefix, samplename, dpclust_postfix), purities, hum_genome_mb)
  all_metrics = rbind(all_metrics, metrics)
}

#############################################################
# Metric 2 - Hom dels
#############################################################
homdel_samples = all_metrics[all_metrics$homdel_mb>=LARGEST_HOM_DEL_ALLOWED_METRIC,]$samplename
if (length(homdel_samples) > 0) {
  sample_labels = setSampleLabelMultiple(applied_corrections$samplename %in% homdel_samples, "correct_to_wgd", sample_labels, applied_corrections)
}
metric2_samplenames = homdel_samples

#############################################################
# Metric 3 - Subclone close to 0.5*CCF
#############################################################
sel = all_metrics$halve_clone_dist<=CLUSTER_0.5CCF_CLONAL_ALLOWED_SITANCE & !all_metrics$clust_between_clonal_and_halve
sel[is.na(sel)] = F
half_ccf_clone_cluster_samples = all_metrics[sel, ]$samplename
if (length(half_ccf_clone_cluster_samples) > 0) {
  sample_labels = setSampleLabelMultiple(applied_corrections$samplename %in% half_ccf_clone_cluster_samples, "correct_to_wgd", sample_labels, applied_corrections)
}
metric3_samplenames = half_ccf_clone_cluster_samples

#############################################################
# Metric 4 - Odd nr allele specific copy number in tetraploid
#############################################################
sel = all_metrics$largest_segment >= LARGEST_50_50_SUBCLONAL_SEGMENT_ALLOWED
sel[is.na(sel)] = F
odd_nr_copynumber_samples = all_metrics[sel, ]$samplename
if (length(odd_nr_copynumber_samples) > 0) {
  sample_labels = setSampleLabelMultiple(which(applied_corrections$samplename %in% odd_nr_copynumber_samples), "correct_to_wgd", sample_labels, applied_corrections)
}
metric4_samplenames = odd_nr_copynumber_samples

#############################################################
# Metric 5 - After filters above: Single method diploid -> tetraploid
#############################################################
# All samples for which a single method calls it diploid should be tetraploid, set label to correct_wgd
metric5_samplenames = c()
for (samplename in applied_corrections$samplename[is.na(sample_labels)]) {
  r = applied_corrections[applied_corrections$samplename==samplename,2:6]
  r = r[!is.na(r)]
  one_diploid_call = (sum(r == "halve" | r == "subtract_one")==(length(r)-1) & length(r) >= 3) 
  
  if (one_diploid_call) {
    sample_labels = setSampleLabel(samplename, "correct_to_wgd", sample_labels, applied_corrections)
    metric5_samplenames = c(metric5_samplenames, samplename)
  }
}

#############################################################
# Metric 6 - After filters above: accept tetraploid call from Broad, Peifer and Van Loo/Wedge
#############################################################

get_method_uniq_adjustments = function(remain, method_index, other_methods_index) {
  method_only = rep(F, nrow(remain))
  for (i in 1:nrow(remain)) {
    single_method_remaining = remain[i, method_index]!="None"
    other_methods = remain[i,other_methods_index]
    
    method_only[i] = (single_method_remaining & all(other_methods[!is.na(other_methods)]=="None"))
  }
  return(method_only)
}

# Make inventory where a single method disagrees
remain = applied_corrections[is.na(sample_labels),]
broad_only = get_method_uniq_adjustments(remain, 2, c(3:6))
dkfz_only = get_method_uniq_adjustments(remain, 3, c(2,4:6))
clonehd_only = get_method_uniq_adjustments(remain, 4, c(2:3,5:6))
peifer_only = get_method_uniq_adjustments(remain, 5, c(2:4,6))
vanloo_wedge_only = get_method_uniq_adjustments(remain, 6, c(2:5))

metric61_samplenames = c()
metric62_samplenames = c()
for (samplename in remain[broad_only | vanloo_wedge_only | peifer_only,]$samplename) {
  sample_labels[applied_corrections$samplename==samplename] = "correct_to_wgd"
  metric61_samplenames = c(metric61_samplenames, samplename)
}

for (samplename in remain[dkfz_only | clonehd_only,]$samplename) {
  sample_labels[applied_corrections$samplename==samplename] = "correct_to_diploid"
  metric62_samplenames = c(metric62_samplenames, samplename)
}

#############################################################
# Save the output
#############################################################
# View(applied_corrections[is.na(sample_labels),])

samples_to_be_corrected_wgd = !is.na(sample_labels) & sample_labels=="correct_to_wgd" #& sample_labels!="no_correction"
samples_to_be_corrected_diploid = !is.na(sample_labels) & sample_labels=="correct_to_diploid"
colnames(applied_corrections) = c("samplename", paste0("input_", colnames(applied_corrections)[2:ncol(applied_corrections)]))

applied_corrections$corrections_broad = NA
sel = !is.na(applied_corrections$input_broad) & applied_corrections$input_broad=="None"
applied_corrections$corrections_broad[samples_to_be_corrected_wgd & sel] = "correct_to_wgd"
sel = !is.na(applied_corrections$input_broad) & applied_corrections$input_broad!="None"
applied_corrections$corrections_broad[samples_to_be_corrected_diploid & sel] = "correct_to_diploid"

applied_corrections$corrections_dkfz = NA
sel = !is.na(applied_corrections$input_dkfz) & applied_corrections$input_dkfz=="None"
applied_corrections$corrections_dkfz[samples_to_be_corrected_wgd & sel] = "correct_to_wgd"
sel = !is.na(applied_corrections$input_dkfz) & applied_corrections$input_dkfz!="None"
applied_corrections$corrections_dkfz[samples_to_be_corrected_diploid & sel] = "correct_to_diploid"

applied_corrections$corrections_mustonen095 = NA
sel = !is.na(applied_corrections$input_mustonen095) & applied_corrections$input_mustonen095=="None"
applied_corrections$corrections_mustonen095[samples_to_be_corrected_wgd & sel] = "correct_to_wgd"
sel = !is.na(applied_corrections$input_mustonen095) & applied_corrections$input_mustonen095!="None"
applied_corrections$corrections_mustonen095[samples_to_be_corrected_diploid & sel] = "correct_to_diploid"

applied_corrections$corrections_peifer = NA
sel = !is.na(applied_corrections$input_peifer) & applied_corrections$input_peifer=="None"
applied_corrections$corrections_peifer[samples_to_be_corrected_wgd & sel] = "correct_to_wgd"
sel = !is.na(applied_corrections$input_peifer) & applied_corrections$input_peifer!="None"
applied_corrections$corrections_peifer[samples_to_be_corrected_diploid & sel] = "correct_to_diploid"

applied_corrections$corrections_vanloo_wedge = NA
sel = !is.na(applied_corrections$input_vanloo_wedge) & applied_corrections$input_vanloo_wedge=="None"
applied_corrections$corrections_vanloo_wedge[samples_to_be_corrected_wgd & sel] = "correct_to_wgd"
sel = !is.na(applied_corrections$input_vanloo_wedge) & applied_corrections$input_vanloo_wedge!="None"
applied_corrections$corrections_vanloo_wedge[samples_to_be_corrected_diploid & sel] = "correct_to_diploid"

triggered_metrics = as.data.frame(matrix(F, nrow(applied_corrections), 6))
colnames(triggered_metrics) = c(paste0("metric_", c(2:5, 6.1, 6.2)))
triggered_metrics$metric_2 = applied_corrections$samplename %in% metric2_samplenames
triggered_metrics$metric_3 = applied_corrections$samplename %in% metric3_samplenames
triggered_metrics$metric_4 = applied_corrections$samplename %in% metric4_samplenames
triggered_metrics$metric_5 = applied_corrections$samplename %in% metric5_samplenames
triggered_metrics$metric_6.1 = applied_corrections$samplename %in% metric61_samplenames
triggered_metrics$metric_6.2 = applied_corrections$samplename %in% metric62_samplenames



applied_corrections$consensus_purity = NA
applied_corrections$consensus_ploidy = NA
consensus_proportion = read.table(corrections_file, header=T, stringsAsFactors=F)$consensus_proportion
for (samplename in applied_corrections$samplename) {
  sample_index = which(applied_corrections$samplename==samplename)
  vlw_input_correction = applied_corrections[sample_index, "input_vanloo_wedge"]
  vlw_output_correction = applied_corrections[sample_index, "corrections_vanloo_wedge"]
  if (is.na(vlw_input_correction) | consensus_proportion[sample_index] < CONSENSUS_PROPORTION) {
    # Currently no way to calculate the hint as there is no Battenberg output
    next
  }
  else if (is.na(vlw_output_correction)) {
    # take Battenberg solution as hint - if there is one
    if (sum(purities$samplename==samplename) > 0) {
      applied_corrections$consensus_purity[sample_index] = purities[purities$samplename==samplename, "purity"]
      applied_corrections$consensus_ploidy[sample_index] = purities[purities$samplename==samplename, "ploidy"]
    }
  } else {
    # Battenberg corrected, so adjust the purity/ploidy
    # read in the baflogr and the segments of the tetraploid consensus
    baf_file = paste0("battenberg_baf_logr/", samplename, "_baflogr.txt.gz")
    segments_file = paste0(interim_tetraploid_prefix, "/", samplename, "_segments.txt")
    
    if (file.exists(segments_file) & file.exists(baf_file)) {
      baflogr = makeGRangesFromDataFrame(read.table(baf_file, header=T, stringsAsFactors=F), keep.extra.columns=T)
      segments = makeGRangesFromDataFrame(read.table(segments_file, header=T, stringsAsFactors=F), keep.extra.columns=T)
      segments$size = end(segments)-start(segments)
      largest_segment = which.max(segments$size)
      
      # Obtain the BAF/logR of the largest segment
      overlap = findOverlaps(segments[largest_segment,], baflogr)
      max_overlap = which.max(width(intersect(segments[largest_segment,], baflogr[subjectHits(overlap),])))
      segment_baf = as.data.frame(baflogr[subjectHits(overlap)[max_overlap],])$baf
      segment_logr = as.data.frame(baflogr[subjectHits(overlap)[max_overlap],])$logr
  
      # Calculate the purity and ploidy
      purity_ploidy = calcPurityPloidy(as.data.frame(segments[largest_segment,])$major_cn, 
                                       as.data.frame(segments[largest_segment,])$minor_cn, 
                                       segment_baf, 
                                       segment_logr)
      
      applied_corrections$consensus_purity[sample_index] = purity_ploidy$purity
      applied_corrections$consensus_ploidy[sample_index] = purity_ploidy$ploidy
    } else {
      print(paste0("No baflogr or segments found for ", samplename))
    }
  }
}

output = data.frame(applied_corrections, triggered_metrics)
write.table(output, file=output_file, row.names=F, quote=F, sep="\t")



