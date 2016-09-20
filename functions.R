
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

calc_metrics = function(diploid_cnprofile, tetraploid_cnprofile, samplename, clust_locs_table_file, purities) {
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