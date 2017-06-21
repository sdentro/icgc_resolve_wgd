library(readr)
library(GenomicRanges)
source("~/repo/icgc_resolve_wgd/functions.R")

#############################################################
# Init
#############################################################
# setwd("/Users/sd11/Documents/Projects/icgc/consensus_SNVs/assessment_consCNA_through_consSNV_6")
setwd("/Users/sd11/Documents/Projects/icgc/consensus_clonal_copynumber/201609/")

# output_file = "icgc_copynumber_wgd_adjustments.txt"
output_file = "interim_consensus_20161011/icgc_copynumber_wgd_adjustments.txt"
purities_file = "interim_consensus/icgc_consensus_dp_master_file.txt" # file in summary table format
corrections_file = "interim_consensus/correction_summary.20160923.txt"
hum_genome_mb = 3235

# Jeffs interim consensus files
# interim_diploid_prefix = "interim_consensus/consensus.best_correction_combo_thresh3/"
interim_diploid_prefix = "interim_consensus_20161011/consensus.best_correction_combo_thresh3/"
# interim_tetraploid_prefix = "interim_consensus/consensus.discard_diploid/"
interim_tetraploid_prefix = "interim_consensus_20161011/consensus.discard_diploid/"
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
  
  #' Note: 080ecc31-756a-4a1b-a51e-d632ac8219f7 All metrics NA
  #' 0bfd1043-7343-fdd0-e050-11ac0c484cab No clustering result?
  #' 0bfd1068-3fcf-a95b-e050-11ac0c4860c3 No clustering result?
  
  
}

# Get a table that contains whether a method was corrected or not
all_metrics_corrections = do.call(rbind, lapply(as.character(all_metrics$samplename), function(samplename) { applied_corrections[as.character(applied_corrections$samplename)==samplename, !grepl("samplename", colnames(applied_corrections))] }))
for (i in 1:nrow(all_metrics_corrections)) {
  all_metrics_corrections[i, is.na(all_metrics_corrections[i,])] = F
}

for (i in 1:nrow(all_metrics_corrections)) {
  all_metrics_corrections[i, all_metrics_corrections[i,]!="None"] = T
  all_metrics_corrections[i, all_metrics_corrections[i,]=="None"] = F
}

for (i in 1:ncol(all_metrics_corrections)) {
  all_metrics_corrections[,i] = as.logical(all_metrics_corrections[,i])
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

#' Note: This should only be triggered if the sample that has been triggered is also corrected. It therefore needs to be aware of the individual CN profiles
#' 
#' Removed after adjustments
#' 03ad38a6-0902-4aaa-84a3-91ea88fa9883 <- overrule dkfz
#' 0448206f-3ade-4087-b1a9-4fb2d14e1367 <- dkfz, should be overruled
#' 0624eb1d-3aff-4037-a3c5-fc363a9edd02 <- dkfz overrule (large segment)
#' 0e54cea2-d568-4a33-b9db-b698844e6ad9 <- dkfz different solution, overrule
#' 0ed2e2e1-2fe3-43eb-8cad-34f3f21a7169 <- sclust overrule
#' 
#' Still there
#' 0168a2a6-c3af-4d58-a51c-d33f0fc7876d <- BB overrule, do different refit suggestion: chr3 30M 2+0
#' 650fe009-da01-4717-89df-9c95fafe3d7e <- BB overrule, do different refit suggestion: chr3 30M 2+1
#' 
#' 
#' kep = c("02c97e2b-914e-4afc-bf50-78f0cfbfa67b", "0332b017-17d5-4083-8fc4-9d6f8fdbbbde", "03b5268e-881e-49e7-824f-170c3fc8b11b", "09508a0d-ebe0-4fa1-b7b2-1710814181cd", "09cb8bc5-13ac-44ac-9b7d-6de143373570", "0a9c9db0-c623-11e3-bf01-24c6515278c0","0bfd1043-7ed6-9ccc-e050-11ac0c481957","0bfd1043-8170-e3e4-e050-11ac0c4860c5","0bfd1043-8172-e3e4-e050-11ac0c4860c5","148536ce-ee2a-4952-a19d-10d6f44146b9","43dadc68-c623-11e3-bf01-24c6515278c0")
#' View(all_metrics[all_metrics$samplename %in% kep, c("samplename", "largest_segment", "total_genome", "frac_rep_genome", "frac_genome", "total_five_largest", "avg_five_largest")])
#' 
#' overrul = c("0168a2a6-c3af-4d58-a51c-d33f0fc7876d", "03ad38a6-0902-4aaa-84a3-91ea88fa9883", "0448206f-3ade-4087-b1a9-4fb2d14e1367", "0624eb1d-3aff-4037-a3c5-fc363a9edd02", "0e54cea2-d568-4a33-b9db-b698844e6ad9", "0ed2e2e1-2fe3-43eb-8cad-34f3f21a7169", "650fe009-da01-4717-89df-9c95fafe3d7e")
#' View(all_metrics[all_metrics$samplename %in% overrul, c("samplename", "largest_segment", "total_genome", "frac_rep_genome", "frac_genome", "total_five_largest", "avg_five_largest")])
#' 
#' from pass Also removed
#' 09508a0d-ebe0-4fa1-b7b2-1710814181cd <- accept (only aceseq calls as diploid), also triggered by metric 5
#' 09cb8bc5-13ac-44ac-9b7d-6de143373570 <- accept BB solution (mcn peak at 1), also triggered by metric 3
#' 
#' 02c97e2b-914e-4afc-bf50-78f0cfbfa67b <- accept BB solution, correct all others
#' 0332b017-17d5-4083-8fc4-9d6f8fdbbbde <- accept ABSOLUTE, refit BB as chr4 100M 3+1
#' 03b5268e-881e-49e7-824f-170c3fc8b11b <- accept (only clonehd calls as diploid)
#' 0a9c9db0-c623-11e3-bf01-24c6515278c0 <- accept
#' 0bfd1043-7ed6-9ccc-e050-11ac0c481957 <- accept
#' 0bfd1043-8170-e3e4-e050-11ac0c4860c5 <- accept
#' 0bfd1043-8172-e3e4-e050-11ac0c4860c5 <- accept
#' 148536ce-ee2a-4952-a19d-10d6f44146b9 <- accept (multiple methods, overrule diploid aceseq)
#' 43dadc68-c623-11e3-bf01-24c6515278c0 <- accept
#' 
#' 09537dce-c797-4b60-962a-d4c3cd6ab00a <- not sure, asked david
#' 
#' 
#' 
#' exta BB refits
#' 
#' 005794f1-5a87-45b5-9811-83ddf6924568 3 30M 1 0
#' 0168a2a6-c3af-4d58-a51c-d33f0fc7876d 3 30M 2 0
#' 650fe009-da01-4717-89df-9c95fafe3d7e 3 30M 2 1
#' 0332b017-17d5-4083-8fc4-9d6f8fdbbbde 4 100M 3 1
#' 01dc6872-c623-11e3-bf01-24c6515278c0 8 30M 1 0
#' 0bfd1043-7343-fdd0-e050-11ac0c484cab 8 30M 2 0
#' 0bfd1043-7346-fdd0-e050-11ac0c484cab 1 50M 2 0
#' 2ccd028d-e7e0-4f77-a512-f658a31819a4 remove refit suggestion
#' 3e7a30e6-2202-40bd-bfc1-0383604050da dont know, very small segments
#' 

#' Further looked at samples
#' 
#' Adding remove DKFZ only removes these (among others)
#' 01df36af-3617-40fc-9892-f54ce433cf71 overrule dkfz
#' 03ad38a6-0902-4aaa-84a3-91ea88fa9883 overrule dkfz
#' 0448206f-3ade-4087-b1a9-4fb2d14e1367 overrule dkfz
#' 33ea81f2-db2c-4567-bd7b-4cb9aadfef88 overrule clonehd
#' 
#' 
#' Rest
#' 00db4dc2-3ec7-4ff9-9233-d69c8c8a607f overrule clonehd
#' 01dc6872-c623-11e3-bf01-24c6515278c0 overrule bb (or refit)
#' 0567d3e6-6278-4d0a-81ae-c084d73c6dd3 overrule bb (no aberrations, yet fit as tetraploid ?)
#' 0624eb1d-3aff-4037-a3c5-fc363a9edd02 overrule dkfz
#' 080ecc31-756a-4a1b-a51e-d632ac8219f7 only broad calls, may be able to fit more clonal with this profile, but not sure
#' 0bfd1043-7343-fdd0-e050-11ac0c484cab accept broad
#' 0bfd1043-7346-fdd0-e050-11ac0c484cab accept clonehd
#' 0bfd1068-3fd3-a95b-e050-11ac0c4860c3 could see a wgd solution work here, not sure whether to accept broad here
#' 0bfd1068-3fe1-a95b-e050-11ac0c4860c3 overrule clonehd
#' 0ed2e2e1-2fe3-43eb-8cad-34f3f21a7169 overrule martin (PRAD heavy wave probably)
#' 126ee433-d345-4cac-882a-c91831a24690 overrule martin (PRAD heavy wave probably)
#' 1d91f9c7-67ba-4606-9f0a-01ec6fc08262 overrule clonehd (PRAD heavy wave probably)
#' 1f5e70c1-c5de-49e7-941a-46e11a4f4416 dont know, BB would fit chr16 with subcl hom del of some size, but doesn't look obviously tetraploid
#' 2399ab13-abfa-480e-9fda-7947edc420be overrule clonehd
#' 2ccd028d-e7e0-4f77-a512-f658a31819a4 overrule bb
#' 3b41cb48-c623-11e3-bf01-24c6515278c0 accept bb
#' 3db3b7b1-da1d-4b9c-a92a-c60fecf4328c not much advantage, overrule broad
#' 3e7a30e6-2202-40bd-bfc1-0383604050da accept broad
#' 42af8f74-fd4b-486d-bc11-db53cc471d62 overrule clonehd
#' 
#' 
#' 
#' Multiple methods disagreeing
#' 09e1fe3e-bfd8-4175-ac42-0e1bf0ba5523 (2 wgd) dont see the benefit of adding the wgd here, overrule broad and martin
#' 0cd60b96-eb2d-4687-9709-d1455ec45de7 (3 different fits) dont see why not a wgd here, should accept (overrule clonehd and martin)
#' 0ead45d8-d785-4404-8319-2ef951e02e03 (2 methods call wgd) think accept based on bb 1+0 and almost 2+0 segments (would become subcl hom dels (largest ~20Mb)), overrule broad and dkfz
#' 1c00925b-7328-4db0-b930-04aab2d80719 (3 call wgd) accept based on chr8 large 1+0 segment, overrule broad and dkfz
#' 1c10ab52-01a3-11e4-8395-af1f6b7ba88c (2 call wgd) see no reason to add wgd, overrule clonehd and martin
#' 20e02396-e676-412d-9724-44a428919cdb (3 call wgd) accept, overrule clonehd and dkfz
#' 25f07374-313a-4100-9a60-3d21d2988fca (2 call wgd) accept, dont see how diploid would fit, overrule broad and dkfz
#' 3e4d0e50-8cf4-4eb0-a00a-ccf0484ecc2f (3 different fits), wgd makes sense, but not the additional one from martin. Overrule broad and dkfz, and martin
#' 418e916b-7a4e-4fab-8616-15dcec4d79f8 (3 call wgd), accept, overrule broad and dkfz
#' 
#' 
#' 
#' 0e6654c9-cd5e-4f94-a6d4-54f2bb16de1f overrule clonehd
#' 1d4a091d-fe65-49c0-8810-5a95243b108a overrule broad
#' 
#' got to 322


# Disable for now
#' # sel = all_metrics$largest_segment >= LARGEST_50_50_SUBCLONAL_SEGMENT_ALLOWED
#' sel = all_metrics$total_five_largest >= 22
#' sel[is.na(sel)] = F
#' 
#' #' Add a further filter to exclude if it's only ACEseq or Clonehd
#' only_dkfz = (all_metrics_corrections$dkfz & rowSums(all_metrics_corrections[,!grepl("dkfz", colnames(all_metrics_corrections))])==0)
#' # only_clonehd = (all_metrics_corrections$mustonen095 & rowSums(all_metrics_corrections[,!grepl("mustonen", colnames(all_metrics_corrections))])==0)
#' sel = sel & !only_dkfz #& !only_clonehd
#' 
#' odd_nr_copynumber_samples = all_metrics[sel, ]$samplename
#' if (length(odd_nr_copynumber_samples) > 0) {
#'   sample_labels = setSampleLabelMultiple(which(applied_corrections$samplename %in% odd_nr_copynumber_samples), "correct_to_wgd", sample_labels, applied_corrections)
#' }
#' metric4_samplenames = odd_nr_copynumber_samples


#############################################################
# Metric TEMP - accept all where ABSOLUTE and Battenberg agree
#############################################################
# Disabled as this metric is not capturing any characteristics of the profile

# length(get_samplenames_with_no_solution(all_metrics, applied_corrections, sample_labels))
# br = applied_corrections$broad!="None"
# br[is.na(br)] = F
# vlw = applied_corrections$vanloo_wedge!="None"
# vlw[is.na(vlw)] = F
# 
# sel = br & vlw & (!is.na(applied_corrections$broad) & applied_corrections$broad!="None")
# 
# # The reverse does not yield any cases
# # applied_corrections[br & vlw & (!is.na(applied_corrections$broad) & applied_corrections$broad=="None"),]
# 
# metricTEMP_samplenames = applied_corrections[sel,]$samplename
# if (length(metricTEMP_samplenames) > 0) {
#   sample_labels = setSampleLabelMultiple(which(applied_corrections$samplename %in% metricTEMP_samplenames), "correct_to_wgd", sample_labels, applied_corrections)
# }


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
# Disabled as this metric is not capturing any characteristics of the profile

# Make inventory where a single method disagrees
# remain = applied_corrections[is.na(sample_labels),]
# broad_only = get_method_uniq_adjustments(remain, 2, c(3:6))
# dkfz_only = get_method_uniq_adjustments(remain, 3, c(2,4:6))
# clonehd_only = get_method_uniq_adjustments(remain, 4, c(2:3,5:6))
# peifer_only = get_method_uniq_adjustments(remain, 5, c(2:4,6))
# vanloo_wedge_only = get_method_uniq_adjustments(remain, 6, c(2:5))
# 
# metric61_samplenames = c()
# metric62b_samplenames = c()
# for (samplename in remain[broad_only | vanloo_wedge_only | peifer_only,]$samplename) {
#   sample_labels[applied_corrections$samplename==samplename] = "correct_to_wgd"
#   metric61_samplenames = c(metric61_samplenames, samplename)
# }
# 
# for (samplename in remain[dkfz_only | clonehd_only,]$samplename) {
#   sample_labels[applied_corrections$samplename==samplename] = "correct_to_diploid"
#   metric62b_samplenames = c(metric62b_samplenames, samplename)
# }

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
triggered_metrics$metric_4 = NA #applied_corrections$samplename %in% metric4_samplenames
triggered_metrics$metric_5 = applied_corrections$samplename %in% metric5_samplenames
triggered_metrics$metric_6.1 = NA #applied_corrections$samplename %in% metric61_samplenames
triggered_metrics$metric_6.2 = NA #applied_corrections$samplename %in% metric62_samplenames



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



