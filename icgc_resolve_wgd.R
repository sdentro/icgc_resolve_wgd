# load_purites_adjustments = function(purities_file) {
#   dat = read_tsv(purities_file)
# 
#   known_methods = c("broad","mustonen095","peifer", "vanloo_wedge")
# 
#   # Get overview of which method was corrected for which sample
#   methods = data.frame(matrix(NA, ncol=length(known_methods), nrow=nrow(dat)))
#   corrections = data.frame(matrix(NA, ncol=length(known_methods), nrow=nrow(dat)))
#   colnames(corrections) = known_methods
# 
#   for (i in 1:nrow(dat)) {
#     if (i %% 100 == 0) { print(paste0(i, " / ", nrow(dat))) }
#     curr_methods = unlist(strsplit(dat$methods[i], ","))
#     curr_corrections = unlist(strsplit(dat$corrections[i], ","))
#     names(curr_corrections) = curr_methods
#     indices = match(curr_methods, colnames(corrections))
#     for (j in 1:length(indices)) {
#       if (!is.na(indices[j])) {
#         corrections[i,indices[j]] = curr_corrections[j]
#         methods[i,indices[j]] = curr_methods[j]
#       }
#     }
#   }
#   return(list(purity_table=dat, methods=methods, corrections=corrections, samplenames=dat$tumor))
# }

createPng = function(p, filename, height, width) { png(filename=filename, height=height, width=width); print(p); dev.off() }

parseAppliedCorrections = function(corrections_file) {
  dat = read.table(corrections_file, header=F, stringsAsFactors=F)
  
  known_methods = c("broad","dkfz","mustonen095","peifer","vanloo_wedge")
  corrections = data.frame(matrix(NA, ncol=length(known_methods), nrow=nrow(dat)))
  colnames(corrections) = known_methods
  
  for (i in 1:nrow(dat)) {
    if (i %% 100 == 0) { print(paste0(i, " / ", nrow(dat))) }
    entries = dat[i,2]
    entries = unlist(strsplit(entries, ","))
    
    splitEntry = function(entry, known_methods) { s = unlist(strsplit(item, "=")); return(list(column=which(known_methods==s[1]), value=s[2])) }
    
    for (item in entries) {
      r = splitEntry(item, known_methods)
      corrections[i, r$column] = r$value
    }
  }
  
  return(data.frame(samplename=dat[,1], corrections))
}

# No longer needed after standardising the breakpoints
# getPreCalcNrBreakpoints = function() {
#   nr.methods = 5
#   broad = read.table("segcounts/broad_segcounts.tsv", header=F, stringsAsFactors=F)
#   dkfz = read.table("segcounts/dkfz_segcounts.tsv", header=F, stringsAsFactors=F)
#   mustonen = read.table("segcounts/mustonen095_segcounts.tsv", header=F, stringsAsFactors=F)
#   peifer = read.table("segcounts/peifer_segcounts.tsv", header=F, stringsAsFactors=F)
#   vanloo_wedge = read.table("segcounts/vanloo_wedge_segcounts.tsv", header=F, stringsAsFactors=F)
#   
#   known_samplenames = Reduce(union, list(b=broad[,1], d=dkfz[,1], m=mustonen[,1], p=peifer[,1], v=vanloo_wedge[,1]))
#   
#   breakpoint_counts = data.frame(matrix(NA, nrow=length(known_samplenames), ncol=nr.methods+1))
#   colnames(breakpoint_counts) = c("samplename", "broad", "dkfz", "mustonen", "peifer", "vanloo_wedge")
#   breakpoint_counts$samplename = known_samplenames
#   
#   
#   breakpoint_counts[match(broad[,1], breakpoint_counts$samplename), "broad"] = broad[,2]
#   breakpoint_counts[match(dkfz[,1], breakpoint_counts$samplename), "dkfz"] = dkfz[,2]
#   breakpoint_counts[match(mustonen[,1], breakpoint_counts$samplename), "mustonen"] = mustonen[,2]
#   breakpoint_counts[match(peifer[,1], breakpoint_counts$samplename), "peifer"] = peifer[,2]
#   breakpoint_counts[match(vanloo_wedge[,1], breakpoint_counts$samplename), "vanloo_wedge"] = vanloo_wedge[,2]
#   
#   return(breakpoint_counts)
# }

#' Check for no cluster between 1 and 0.5 cluster
getClusterMetricsSample = function(filename, samplename, purities) {
  
  if (file.exists(filename)) {
    
    clusters = readr::read_tsv(filename)
    clusters$ccf = clusters$proportion / purities[purities$tumor==samplename, "purity_consensus"]
    clonal = which.min(abs(clusters$ccf - 1))
    
    dist_to_halve = abs(clusters$ccf - (clusters$ccf[clonal] * 0.5))
    dist_to_quarter = abs(clusters$ccf - (clusters$ccf[clonal] * 0.25))
    
    closest_to_halve = which.min(dist_to_halve)
    closest_to_quarter = which.min(dist_to_quarter)
    
    clust_between_clonal_and_halve = !((clonal-closest_to_halve)==1)
    output = data.frame(samplename=samplename, 
                        clone_ccf=clusters$ccf[clonal],
                        clone_cluster=clusters$cluster[clonal], 
                        halve_clone_ccf=clusters$ccf[closest_to_halve], 
                        halve_clone_cluster=clusters$cluster[closest_to_halve], 
                        halve_clone_dist=dist_to_halve[closest_to_halve], 
                        quart_clone_ccf=clusters$ccf[closest_to_quarter], 
                        quart_clone_cluster=clusters$cluster[closest_to_quarter], 
                        quart_clone_dist=dist_to_quarter[closest_to_quarter], 
                        clust_between_clonal_and_halve=clust_between_clonal_and_halve)
  } else {
    output = data.frame(samplename=samplename, 
                        clone_ccf=NA,
                        clone_cluster=NA, 
                        halve_clone_ccf=NA, 
                        halve_clone_cluster=NA, 
                        halve_clone_dist=NA, 
                        quart_clone_ccf=NA, 
                        quart_clone_cluster=NA, 
                        quart_clone_dist=NA, 
                        clust_between_clonal_and_halve=NA)
  }
  return(output)
}

getClusterMetrics = function(path_prefix, vector_of_samplenames, dat) {
  res = lapply(1:length(vector_of_samplenames), function(i) {
    getClusterMetricsSample(paste0(path_prefix, vector_of_samplenames[i], "_subclonal_structure.txt.gz"),
                            vector_of_samplenames[i],
                            dat)
  })
  res = do.call(rbind, res)
  return(res)
}

calc_metrics = function(diploid_cnprofile, tetraploid_cnprofile, samplename, clust_locs_table_file) {
  if (file.exists(diploid_cnprofile)) {
    diploid_cnprofile = read.table(diploid_cnprofile, header=T, stringsAsFactors=F)
    
    # Get homozygous deletions
    diploid_cnprofile$length = (diploid_cnprofile$end-diploid_cnprofile$start)/1000000
    homdel = diploid_cnprofile[diploid_cnprofile$major_cn==0 & diploid_cnprofile$minor_cn==0,]$length
    if (length(homdel)==0) {
      homdel = 0
    } else {
      homdel = max(homdel)
    }
  } else {
    homdel = NA
  }
   
  
  # Get 0.5*CCF metrics
  clustering_metrics = getClusterMetricsSample(clust_locs_table_file, samplename, purities)
  
  
  if (file.exists(tetraploid_cnprofile)) {
    tetraploid_cnprofile = read.table(tetraploid_cnprofile, header=T, stringsAsFactors=F)
    
    tetraploid_cnprofile$seg_length = (tetraploid_cnprofile$end-tetraploid_cnprofile$start)/1000000
    odd_segs = tetraploid_cnprofile[tetraploid_cnprofile$major_cn %% 2 == 1 | tetraploid_cnprofile$minor_cn %% 2 == 1, ]
    
    five_largest = order(odd_segs$seg_length, decreasing=T)[1:5]
    tetra_cons = data.frame(largest_segment=max(odd_segs$seg_length),
                      total_genome=sum(odd_segs$seg_length),
                      frac_rep_genome=sum(odd_segs$seg_length)/sum(tetraploid_cnprofile$seg_length),
                      frac_genome=sum(odd_segs$seg_length)/hum_genome_mb,
                      total_five_largest=sum(five_largest),
                      avg_five_largest=mean(five_largest))
  } else {
    tetra_cons = data.frame(largest_segment=NA,
                            total_genome=NA,
                            frac_rep_genome=NA,
                            frac_genome=NA,
                            total_five_largest=NA,
                            avg_five_largest=NA)
  }
  
  return(data.frame(clustering_metrics, tetra_cons, homdel_mb=homdel))
}

setSampleLabel = function(samplename, label, sample_labels, applied_corrections) {
  sample_labels[applied_corrections$samplename==samplename] = label
  return(sample_labels)
}

setSampleLabelMultiple = function(selection_vector, label, sample_labels, applied_corrections) {
  sample_labels[selection_vector] = label
  return(sample_labels)
}


parseReplicationTimingData = function() {
  tumourR2 = read.table("broad_replication_timing_2016_06/Tumor_R2_with_project.tsv", header=T, stringsAsFactors=F)
  normalR2 = read.table("broad_replication_timing_2016_06/Normals of CORR_COEFF_2.tsv", header=T, stringsAsFactors=F)
  return(list(tumourR2=tumourR2, normalR2=normalR2))
}

#############################################################
# Init
#############################################################
library(readr)
setwd("/Users/sd11/Documents/Projects/icgc/consensus_SNVs/assessment_consCNA_through_consSNV_6")

output_file = "icgc_copynumber_wgd_adjustments.txt"
purities_file = "interim_consensus/consensus_purities.20160319.txt"
corrections_file = "applied_corrections_gt30_agreement.txt"
hum_genome_mb = 3235

# Jeffs interim consensus files
interim_diploid_prefix = "interim_consensus/clonal_consensus_001/"
interim_tetraploid_prefix = "interim_consensus/consensus.discard_diploid/"
segments_file_postfix = "_segments.txt"

# DPClust files
dpclust_prefix = "dpclust_clust_locs/"
dpclust_postfix = "_subclonal_structure.txt.gz"

# purities = load_purites_adjustments(purities_file)
applied_corrections = parseAppliedCorrections(corrections_file)

# Create the output vector with a statement what should be done with the tumour
sample_labels = rep(NA, nrow(applied_corrections))

# Thresholds for various metrics
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
  
  metrics = calc_metrics(diploid_cnprofile, tetraploid_cnprofile, samplename, paste0(dpclust_prefix, samplename, dpclust_postfix))
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
View(applied_corrections[is.na(sample_labels),])

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

matches = match(applied_corrections$samplename, purities$tumor)
consensus_purities = purities[matches, c("tumor", "purity_median", "uncorrected_purities")]

# Adjust to tetraploid take the median of all corrected cases
consensus_purities$consensus_purity = NA

# TODO: This code does not yet work with the current adjustment table
# adjust_to_wgd = applied_corrections$samplename[!is.na(sample_labels) & sample_labels=="correct_to_wgd"]
# # samplename = "fc95d5ce-6899-62f1-e040-11ac0c486011"
# for (samplename in adjust_to_wgd) {
#   sample_corrections = applied_corrections[applied_corrections$samplename==samplename, grepl("input", colnames(applied_corrections))]
#   uncorrected_purities_index = which(!is.na(sample_corrections) & sample_corrections!="None")
#   uncorrected_purities = as.numeric(unlist(strsplit(purities[purities$tumor==samplename, c("uncorrected_purities")], ",")))
#   consensus_purities$consensus_purity[consensus_purities$tumor==samplename] = median(uncorrected_purities[uncorrected_purities_index])
# }
# 
# adjust_to_diploid = applied_corrections$samplename[!is.na(sample_labels) & sample_labels=="correct_to_diploid"]
# # samplename = adjust_to_diploid[1]
# for (samplename in adjust_to_diploid) {
#   sample_corrections = applied_corrections[applied_corrections$samplename==samplename, grepl("input", colnames(applied_corrections))]
#   uncorrected_purities_index = which(!is.na(sample_corrections) & sample_corrections=="None")
#   uncorrected_purities = as.numeric(unlist(strsplit(purities[purities$tumor==samplename, c("uncorrected_purities")], ",")))
#   consensus_purities$consensus_purity[consensus_purities$tumor==samplename] = median(uncorrected_purities[uncorrected_purities_index])
# }
consensus_purities = consensus_purities[match(consensus_purities$tumor, applied_corrections$samplename),]
consensus_purities = consensus_purities[,colnames(consensus_purities)!="tumor"]
output = data.frame(applied_corrections, triggered_metrics, consensus_purities)

write.table(output, file=output_file, row.names=F, quote=F, sep="\t")



