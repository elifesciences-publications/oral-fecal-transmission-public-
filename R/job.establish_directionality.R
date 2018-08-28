#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Establish putative directionality of transmission (oral -> gut vs gut -> oral)
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load SNV tables per taxon
#=> identify "personal SNVs" per subject
#=> quantify longitudinal stability of allele profiles per subject & body site
#=> for each individual, track (asymmetric) overlap of alleles between sites over time
#=> save and export individual datafiles in R format
#
#2017-10-17
#sebastian.schmidt@embl.de
################################################################################


################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
################################################################################
################################################################################


################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.base <- "################"
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.logs <- paste0(PARAM$folder.base, "logs/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Set parameters
#=> minimum number of shared observations to quantify overlap
PARAM$min.shared_obs <- 20

#Read current taxon of interest from command line
args <- commandArgs(TRUE)
tax.id <- args[2]

#Report
writeLines(paste(date(), "=> Processing", tax.id))
################################################################################
################################################################################


################################################################################
################################################################################
#Load allele data
file.allele_data <- paste0(PARAM$folder.SNV_data, tax.id, ".allele_data.filtered.RData")
if (! file.exists(file.allele_data)) {writeLines(paste(tax.id, "file does not exist:", file.allele_data)); q("no")}
load(file.allele_data)

#Load sample data
load(paste0(PARAM$folder.parameters, "data.sample.RData"))
samples <- colnames(mat.allele)
n.samples <- length(samples)
data.sample <- data.sample[samples, ]

#Load timepoint data
load(paste0(PARAM$folder.parameters, "data.timepoint.RData"))

#Load specI abundance data
load(paste0(PARAM$folder.data, "data.specI.RData"))
curr.specI <- data.specI[tax.id, samples]

#Load coverage data
load(paste0(PARAM$folder.data, "data.coverage.RData"))
curr.hor_cov <- data.hor_cov[tax.id, samples]
curr.ver_cov <- data.ver_cov[tax.id, samples]

#Load taxonomy data
load(paste0(PARAM$folder.data, "data.taxonomy.RData"))
################################################################################
################################################################################


################################################################################
################################################################################
#Define a function to obtain 
# get.allele_overlap <- function(pair) {
#   inc.1 <- inc.allele[,pair[1]]; inc.2 <- inc.allele[, pair[2]]
#   union.inc <- which((inc.1 & obs.allele[,pair[2]]) | (inc.2 & obs.allele[,pair[1]]))
#   if (sum(union.inc) < PARAM$min.shared_obs) {return(NA)}
#   shared.inc <- which(inc.1[union.inc] & inc.2[union.inc])
#   length(shared.inc) / length(union.inc)
# }
# ###################################
# 
# ###################################
# #Define a function to detect "longitudinal asymmetries" in allele profiles between sites
# #=> consider two body sites A & B, across timepoints 1 & 2
# detect.longitudinal_asymmetries <- function(obs.A.1, obs.A.2, obs.B.1, obs.B.2, inc.A.1, inc.A.2, inc.B.1, inc.B.2) {
#   #Get alleles with observations in all four samples
#   #=> this is the baseline of alleles to consider in the first place
#   obs.all_shared <- obs.A.1 & obs.A.2 & obs.B.1 & obs.B.2
#   
#   #Get "trivial" cases that are all 0-0 or 1-1 in *both* sites
#   inc.trivial <- obs.all_shared & inc.A.1 == inc.A.2 & inc.B.1 == inc.B.2
#   
#   #Get A->B trends
#   #=> alleles observed in site A but not B at timepoint t0, and observed in site B at timepoint t1
#   inc.A_B <- obs.all_shared & (!inc.trivial) & inc.A.1 & (!inc.B.1) & inc.B.2
#   
#   #Get B->A trends
#   inc.B_A <- obs.all_shared & (!inc.trivial) & (!inc.A.1) & inc.B.1 & inc.A.2
#   
#   return(list(shared_obs=sum(obs.all_shared), trivial_inc=sum(inc.trivial), A_B=sum(inc.A_B), B_A=sum(inc.B_A)))
# }
###################################
################################################################################
################################################################################


################################################################################
################################################################################
#Get current habitat vector
habitats <- data.sample[colnames(mat.allele), "material"]
subjects <- data.sample[colnames(mat.allele), "subject"]
subjects.timepoints <- data.sample[colnames(mat.allele), "subject_timepoint"]

#Get current taxon name
taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
print.taxon <- gsub(" ", "_", taxon)
print.taxon <- gsub("/", "_", print.taxon)

#Generate allele observation and incidence matrices
obs.allele <- mat.allele > 0
inc.allele <- mat.allele > 10^-12
rm(mat.allele); invisible(gc())
###################################

###################################
#Get background frequencies for "allele flips" over time, per time bin
#=> for saliva & stool at t=30, 60 & 90 only (too little data for the other timepoints)
###################################
#Subset data.delta_t
sub.delta_t <- data.delta_t[data.delta_t$site.1 == data.delta_t$site.2, ]
mat.allele.flip <- matrix(data=NA, nrow=nrow(obs.allele), ncol=nrow(sub.delta_t), dimnames=list(rownames(obs.allele), rownames(sub.delta_t)))
for (j in 1:nrow(sub.delta_t)) {
  #Get current samples
  smpl.1 <- sub.delta_t$sample.1[j]
  smpl.2 <- sub.delta_t$sample.2[j]
  
  #Skip unless both samples are available for current pair
  if (!smpl.1 %in% colnames(obs.allele) | !smpl.2 %in% colnames(obs.allele)) {next()}
  
  #Get coupled sites' sample(s) at t0, if available
  curr.coupled <- c(
    data.timepoint[data.timepoint$subject == sub.delta_t$subject[j] & data.timepoint$timepoint == sub.delta_t$timepoint.1[j], "sample.oral"],
    data.timepoint[data.timepoint$subject == sub.delta_t$subject[j] & data.timepoint$timepoint == sub.delta_t$timepoint.1[j], "sample.plaque"],
    data.timepoint[data.timepoint$subject == sub.delta_t$subject[j] & data.timepoint$timepoint == sub.delta_t$timepoint.1[j], "sample.gut"]
  )
  names(curr.coupled) <- c("saliva", "dental_plaque", "stool")
  
  #Get current shared observations across time within individual
  curr.shared_obs <- obs.allele[,smpl.1] & obs.allele[,smpl.2]
  
  #Get current longitudinal profile of alleles (which alleles are stable, lost or gained)
  curr.longitudinal_profile <- inc.allele[curr.shared_obs, smpl.2] - inc.allele[curr.shared_obs, smpl.1]
  mat.allele.flip[curr.shared_obs, j] <- curr.longitudinal_profile
  
  #Store the global allele stability, loss and gain stats for current timepoint pair
  sub.delta_t[j, "n.shared_obs"] <- sum(curr.shared_obs)
  sub.delta_t[j, "n.stable"] <- sum(curr.longitudinal_profile == 0)
  sub.delta_t[j, "n.loss"] <- sum(curr.longitudinal_profile == -1)
  sub.delta_t[j, "n.gain"] <- sum(curr.longitudinal_profile == 1)
  
  #Iterate through potential coupling partners
  for (k in names(curr.coupled)) {
    #Skip within-site comparison
    if (k == sub.delta_t$site.1[j]) {next()}
    #Skip if no coupled sample is available
    if (is.na(curr.coupled[k]) | !curr.coupled[k] %in% colnames(obs.allele)) {next()}
    
    #Get current alleles with shared observations between all three samples
    curr.shared_obs.coupled <- curr.shared_obs & obs.allele[, unname(curr.coupled[k])]
    sub.delta_t[j, paste(k, "n.shared_obs", sep=".")] <- sum(curr.shared_obs.coupled)
    
    #Skip if there is too little information available (too few shared observations)
    if (sum(curr.shared_obs.coupled) < PARAM$min.shared_obs) {next()}
    
    #Get current incidence in coupled sites at t0
    curr.coupled.obs <- inc.allele[curr.shared_obs.coupled, curr.coupled[k]]
    #Get longitudinal stability info for focal alleles
    curr.long_stab <- mat.allele.flip[curr.shared_obs.coupled, j]
    
    #Get a contingency table for current observations
    ct <- matrix(nrow=2, ncol=2)
    #=> allele gain and incidence in coupled site at t0
    ct[1,1] <- sum(curr.coupled.obs & curr.long_stab > 0)
    #=> allele gain, but no incidence in coupled site at t0
    ct[1,2] <- sum(!curr.coupled.obs & curr.long_stab > 0)
    #=> allele stable or lost and incidence at coupled site at t0
    ct[2,1] <- sum(curr.coupled.obs & curr.long_stab <= 0)
    #=> allele stable or lost and no incidence at coupled sitet at t0
    ct[2,2] <- sum(!curr.coupled.obs & curr.long_stab <= 0)
    
    #Perform Fisher's exact test
    curr.test <- fisher.test(ct)
    
    #Store info
    sub.delta_t[j, paste(k, "n.gain.coupled", sep=".")] <- ct[1,1]
    sub.delta_t[j, paste(k, "n.gain.not_coupled", sep=".")] <- ct[1,2]
    sub.delta_t[j, paste(k, "n.not_gained.coupled", sep=".")] <- ct[2,1]
    sub.delta_t[j, paste(k, "odds_ratio", sep=".")] <- curr.test$estimate
    sub.delta_t[j, paste(k, "p", sep=".")] <- curr.test$p.value
  }
}
###################################

###################################
#Store data
save(sub.delta_t, file=paste0(PARAM$folder.results, "establish.directionality/", tax.id, ".", print.taxon,  ".directionality.fisher_test.RData"))

#Report
writeLines(paste(date(), "=> Done with", tax.id, taxon))
################################################################################
################################################################################

q(save="no")



#Subset data to alleles with sufficient number of observations
# na.freq <- rowSums(is.na(mat.allele.flip))
# na.all.idx <- na.freq >= ncol(mat.allele.flip) - 5
# if (sum(na.all.idx) < 100) {writeLines(paste(date(), tax.id, taxon, "=> Too few alleles remain after filtering.")); q(save="no")}
# mat.allele.flip.use <- mat.allele.flip[!na.all.idx, ]
# obs.allele.use <- obs.allele[!na.all.idx, ]
# inc.allele.use <- inc.allele[!na.all.idx, ]
# #Tidy up
# rm(mat.allele); rm(mat.allele.flip); rm(obs.allele); rm(inc.allele); invisible(gc())
# 
# #Get background frequencies for stable, lost and gained alleles
# bg.freq.stable <- rowSums(mat.allele.flip.use == 0, na.rm=T) / (ncol(mat.allele.flip.use) - na.freq[!na.all.idx])
# bg.freq.loss <- rowSums(mat.allele.flip.use == -1, na.rm=T) / (ncol(mat.allele.flip.use) - na.freq[!na.all.idx])
# bg.freq.gain <- rowSums(mat.allele.flip.use == 1, na.rm=T) / (ncol(mat.allele.flip.use) - na.freq[!na.all.idx])
###################################




