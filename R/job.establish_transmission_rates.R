#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Establish rates of transmission (oral -> gut vs gut -> oral)
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load SNV tables per taxon
#=> identify "personal SNVs" per subject
#=> quantify longitudinal stability of allele profiles per subject & body site
#=> for each individual, track (asymmetric) overlap of alleles between sites over time
#=> from this, estimate total strain turnover per day, per site
#=> estimate fraction of strain turnover attributable to transmission (=> relative transmission rate)
#=> from a between-subject background, estimate the random error on strain turnover and transmission rate
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
PARAM$folder.base <- "###############"
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
sub.delta_t <- data.delta_t[data.delta_t$site.1 == data.delta_t$site.2 & data.delta_t$site.1 != "dental_plaque", ]
sub.delta_t.bg.list <- list()
#mat.allele.flip <- matrix(data=NA, nrow=nrow(obs.allele), ncol=nrow(sub.delta_t), dimnames=list(rownames(obs.allele), rownames(sub.delta_t)))
for (j in 1:nrow(sub.delta_t)) {
  #Get current samples
  smpl.1 <- sub.delta_t$sample.1[j]
  smpl.2 <- sub.delta_t$sample.2[j]
  
  #Skip unless both samples are available for current pair
  if (!smpl.1 %in% colnames(obs.allele) | !smpl.2 %in% colnames(obs.allele)) {next()}
  
  writeLines(paste(date(), "=> Processing taxon", taxon, "subject", sub.delta_t$subject[j]))
  
  #Get coupled sites' sample(s) at t0, if available
  curr.coupled <- c(
    data.timepoint[data.timepoint$subject == sub.delta_t$subject[j] & data.timepoint$timepoint == sub.delta_t$timepoint.1[j], "sample.oral"],
    #data.timepoint[data.timepoint$subject == sub.delta_t$subject[j] & data.timepoint$timepoint == sub.delta_t$timepoint.1[j], "sample.plaque"],
    data.timepoint[data.timepoint$subject == sub.delta_t$subject[j] & data.timepoint$timepoint == sub.delta_t$timepoint.1[j], "sample.gut"]
  )
  #names(curr.coupled) <- c("saliva", "dental_plaque", "stool")
  names(curr.coupled) <- c("saliva", "stool")
  
  #Get current shared observations across time within individual
  curr.obs_allele <- obs.allele[,smpl.1]
  curr.inc_allele <- inc.allele[,smpl.1]
  curr.shared_obs <- curr.obs_allele & obs.allele[,smpl.2]
  
  #Get current longitudinal profile of alleles (which alleles are stable, lost or gained)
  curr.flip <- logical(length = nrow(obs.allele))
  curr.longitudinal_profile <- inc.allele[curr.shared_obs, smpl.2] - curr.inc_allele[curr.shared_obs]
  curr.flip[curr.shared_obs] <- curr.longitudinal_profile
  
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
    
    sub.delta_t[j, "source.site"] <- k
    sub.delta_t[j, "sink.site"] <- sub.delta_t$site.1[j]
    if (k == "saliva") {sub.delta_t[j, "source_sink"] <- "saliva_to_stool"} else {sub.delta_t[j, "source_sink"] <- "stool_to_saliva"}
    
    #Get current alleles with shared observations between all three samples
    curr.shared_obs.coupled <- curr.shared_obs & obs.allele[, unname(curr.coupled[k])]
    sub.delta_t[j, "cross.n.shared_obs.coupled"] <- sum(curr.shared_obs.coupled)
    
    #Skip if there is too little information available (too few shared observations)
    if (sum(curr.shared_obs.coupled) < PARAM$min.shared_obs) {next()}
    
    #Get current incidence in coupled sites at t0
    curr.coupled.obs <- inc.allele[curr.shared_obs.coupled, unname(curr.coupled[k])]
    #Get longitudinal stability info for focal alleles
    curr.long_stab <- curr.flip[curr.shared_obs.coupled]
    
    #Get a contingency table for current observations
    ct <- matrix(nrow=2, ncol=2)
    #=> allele gain and incidence in coupled site at t0
    ct[1,1] <- sum(curr.coupled.obs & curr.long_stab > 0)
    #=> allele gain, but no incidence in coupled site at t0
    ct[1,2] <- sum(!curr.coupled.obs & curr.long_stab > 0)
    #=> allele stable or lost and incidence at coupled site at t0
    ct[2,1] <- sum(curr.coupled.obs & curr.long_stab <= 0)
    #=> allele stable or lost and no incidence at coupled site at t0
    ct[2,2] <- sum(!curr.coupled.obs & curr.long_stab <= 0)
    
    #Perform Fisher's exact test
    curr.test <- fisher.test(ct)
    
    #Get current turnover data
    curr.turnover.raw <- 1 - (sum(curr.longitudinal_profile == 0) / sum(curr.shared_obs))
    curr.turnover <- curr.turnover.raw / sub.delta_t$delta_t[j] 
    curr.turnover.transmission.raw = ct[1,1] / (ct[1,1] + ct[1,2])
    curr.turnover.transmission <- curr.turnover * curr.turnover.transmission.raw
    curr.turnover.non_transmission <- curr.turnover - (curr.turnover * curr.turnover.transmission.raw)
    
    #Store info
    sub.delta_t[j, "cross.n.gained.coupled"] <- ct[1,1]
    sub.delta_t[j, "cross.n.gained.not_coupled"] <- ct[1,2]
    sub.delta_t[j, "cross.n.not_gained.coupled"] <- ct[2,1]
    sub.delta_t[j, "odds_ratio"] <- curr.test$estimate
    sub.delta_t[j, "p"] <- curr.test$p.value
    sub.delta_t[j, "turnover.raw"] <- curr.turnover.raw
    sub.delta_t[j, "turnover"] <- curr.turnover
    sub.delta_t[j, "turnover.transmission.raw"] <- curr.turnover.transmission.raw
    sub.delta_t[j, "turnover.transmission"] <- curr.turnover.transmission
    sub.delta_t[j, "turnover.non_transmission.raw"] <- 1 - curr.turnover.transmission.raw
    sub.delta_t[j, "turnover.non_transmission"] <- curr.turnover.non_transmission
    
    #Calculate background for longitudinal profile
    #=> allele overlap between sample at t0 and random samples from same cohort
    bg.samples.longit <- rownames(data.sample)[
      data.sample$material == sub.delta_t$site.1[j] &
        as.character(data.sample$subject) != as.character(sub.delta_t$subject[j]) &
        as.character(data.sample$family) != as.character(sub.delta_t$family[j]) &
        as.character(data.sample$cohort) == as.character(sub.delta_t$cohort[j]) &
        data.sample$timepoint %in% c(0, 7, 600)
    ]
    if (length(bg.samples.longit) >= 50) {bg.samples.longit <- sample(bg.samples.longit, 50, replace = F)}
    #=> background samples for coupled site
    bg.samples.coupled <- rownames(data.sample)[
      data.sample$material == k &
        as.character(data.sample$subject) != as.character(sub.delta_t$subject[j]) &
        as.character(data.sample$family) != as.character(sub.delta_t$family[j]) &
        as.character(data.sample$cohort) == as.character(sub.delta_t$cohort[j]) &
        data.sample$timepoint %in% c(0, 7, 600)
      ]
    if (length(bg.samples.coupled) >= 50) {bg.samples.coupled <- sample(bg.samples.coupled, 50, replace = F)}
    
    #Preallocate
    curr.bg <- data.frame()
    
    #Loop through longitudinal background samples
    for (b in bg.samples.longit) {
      #Shared observations across "time"
      bg.shared_obs <- curr.obs_allele & obs.allele[, b]
      
      #"Longitudinal" profile
      bg.flip <- logical(length = nrow(obs.allele))
      bg.longitudinal_profile <- inc.allele[bg.shared_obs, b] - curr.inc_allele[bg.shared_obs]
      bg.flip[bg.shared_obs] <- bg.longitudinal_profile
      
      #Get current alleles with shared obs btw all three samples
      bg.shared_obs.coupled <- bg.shared_obs & obs.allele[, unname(curr.coupled[k])]
      if (sum(bg.shared_obs.coupled) < PARAM$min.shared_obs) {next()}
      
      #Get current incidence in coupled sites at t0
      bg.coupled.obs <- inc.allele[bg.shared_obs.coupled, unname(curr.coupled[k])]
      #Get longitudinal stability info for focal alleles
      bg.long_stab <- bg.flip[bg.shared_obs.coupled]
      
      #Get a contingency table for current observations
      ct <- matrix(nrow=2, ncol=2)
      #=> allele gain and incidence in coupled site at t0
      ct[1,1] <- sum(bg.coupled.obs & bg.long_stab > 0)
      #=> allele gain, but no incidence in coupled site at t0
      ct[1,2] <- sum(!bg.coupled.obs & bg.long_stab > 0)
      #=> allele stable or lost and incidence at coupled site at t0
      ct[2,1] <- sum(bg.coupled.obs & bg.long_stab <= 0)
      #=> allele stable or lost and no incidence at coupled sitet at t0
      ct[2,2] <- sum(!bg.coupled.obs & bg.long_stab <= 0)
      
      #Perform Fisher's exact test
      bg.test <- fisher.test(ct)
      
      #Get current turnover data
      bg.turnover.raw <- bg.turnover <- 1 - (sum(bg.longitudinal_profile == 0) / sum(bg.shared_obs))
      bg.turnover.transmission.raw = ct[1,1] / (ct[1,1] + ct[1,2])
      bg.turnover.transmission <- bg.turnover * bg.turnover.transmission.raw
      bg.turnover.non_transmission <- bg.turnover - (bg.turnover * bg.turnover.transmission.raw)
      
      #Append data to current data.frame
      curr.bg <- rbind(curr.bg, data.frame(
        subject = as.character(sub.delta_t$subject[j]),
        cohort = as.character(sub.delta_t$cohort[j]),
        family = as.character(sub.delta_t$family[j]),
        site.1 = sub.delta_t$site.1[j],
        site.2 = sub.delta_t$site.1[j],
        timepoint.1 = sub.delta_t$timepoint.1[j],
        timepoint.2 = -1, #for background
        delta_t = -1,
        sample.1 = smpl.1,
        sample.2 = b,
        delta_t.bin = "bg.longitudinal",
        n.shared_obs = sum(bg.shared_obs),
        n.stable = sum(bg.longitudinal_profile == 0),
        n.loss = sum(bg.longitudinal_profile == -1),
        n.gain = sum(bg.longitudinal_profile == 1),
        source.site = k,
        sink.site = sub.delta_t$site.1[j],
        source_sink = sub.delta_t[j, "source_sink"],
        cross.n.shared_obs.coupled = sum(bg.shared_obs.coupled),
        cross.n.gained.coupled = ct[1,1],
        cross.n.gained.not_coupled = ct[1,2],
        cross.n.not_gained.coupled = ct[2,1],
        odds_ratio = bg.test$estimate,
        p = bg.test$p.value,
        turnover.raw = bg.turnover.raw,
        turnover = bg.turnover,
        turnover.transmission.raw = bg.turnover.transmission.raw,
        turnover.transmission = bg.turnover.transmission,
        turnover.non_transmission.raw = 1 - bg.turnover.transmission.raw,
        turnover.non_transmission = bg.turnover.non_transmission
      ))
    }
    
    #Loop through coupled background samples
    for (b in bg.samples.coupled) {
      #Get current alleles with shared obs btw all three samples
      bg.shared_obs.coupled <- curr.shared_obs & obs.allele[, b]
      if (sum(bg.shared_obs.coupled) < PARAM$min.shared_obs) {next()}
      
      #Get current incidence in randomly picked coupled site
      bg.coupled.obs <- inc.allele[bg.shared_obs.coupled, b]
      #Get longitudinal stability info for focal alleles
      bg.long_stab <- curr.flip[bg.shared_obs.coupled]
      
      #Get a contingency table for current observations
      ct <- matrix(nrow=2, ncol=2)
      #=> allele gain and incidence in coupled site at t0
      ct[1,1] <- sum(bg.coupled.obs & bg.long_stab > 0)
      #=> allele gain, but no incidence in coupled site at t0
      ct[1,2] <- sum(!bg.coupled.obs & bg.long_stab > 0)
      #=> allele stable or lost and incidence at coupled site at t0
      ct[2,1] <- sum(bg.coupled.obs & bg.long_stab <= 0)
      #=> allele stable or lost and no incidence at coupled sitet at t0
      ct[2,2] <- sum(!bg.coupled.obs & bg.long_stab <= 0)
      
      #Perform Fisher's exact test
      bg.test <- fisher.test(ct)
      
      #Get current turnover data
      bg.turnover.transmission.raw = ct[1,1] / (ct[1,1] + ct[1,2])
      bg.turnover.transmission <- curr.turnover * bg.turnover.transmission.raw
      bg.turnover.non_transmission <- curr.turnover - (curr.turnover * bg.turnover.transmission.raw)
      
      #Append data to current data.frame
      curr.bg <- rbind(curr.bg, data.frame(
        subject = as.character(sub.delta_t$subject[j]),
        cohort = as.character(sub.delta_t$cohort[j]),
        family = as.character(sub.delta_t$family[j]),
        site.1 = sub.delta_t$site.1[j],
        site.2 = sub.delta_t$site.1[j],
        timepoint.1 = sub.delta_t$timepoint.1[j],
        timepoint.2 = -2, #for background
        delta_t = -2,
        sample.1 = smpl.1,
        sample.2 = smpl.2,
        delta_t.bin = "bg.coupled",
        n.shared_obs = sum(curr.shared_obs),
        n.stable = sum(curr.longitudinal_profile == 0),
        n.loss = sum(curr.longitudinal_profile == -1),
        n.gain = sum(curr.longitudinal_profile == 1),
        source.site = k,
        sink.site = sub.delta_t$site.1[j],
        source_sink = sub.delta_t[j, "source_sink"],
        cross.n.shared_obs.coupled = sum(bg.shared_obs.coupled),
        cross.n.gained.coupled = ct[1,1],
        cross.n.gained.not_coupled = ct[1,2],
        cross.n.not_gained.coupled = ct[2,1],
        odds_ratio = bg.test$estimate,
        p = bg.test$p.value,
        turnover.raw = curr.turnover.raw,
        turnover = curr.turnover,
        turnover.transmission.raw = bg.turnover.transmission.raw,
        turnover.transmission = bg.turnover.transmission,
        turnover.non_transmission.raw = 1 - bg.turnover.transmission.raw,
        turnover.non_transmission = bg.turnover.non_transmission
      ))
    }
    
    #Store current background data
    sub.delta_t.bg.list[[j]] <- curr.bg
  }
}

#Summarise bg calculation results
sub.delta_t.bg <- do.call(rbind, sub.delta_t.bg.list)
###################################

###################################
#Store data
save(sub.delta_t, sub.delta_t.bg, sub.delta_t.bg.list, file=paste0(PARAM$folder.results, "establish.transmission_rates/", tax.id, ".", print.taxon,  ".transmission_rates.RData"))

#Report
writeLines(paste(date(), "=> Done with", tax.id, taxon))
################################################################################
################################################################################

q(save="no")
###################################




