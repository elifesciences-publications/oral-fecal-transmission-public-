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
PARAM$folder.base <- "#############"
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.logs <- paste0(PARAM$folder.base, "logs/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Set parameters
#=> minimum number of shared observations to quantify overlap
PARAM$min.shared_obs <- 100
PARAM$max.combn <- 500
PARAM$max.alleles <- 100000

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
rm(data.allele); invisible(gc())

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
#Define accessory function(s)
get.allele_overlap <- function(pair) {
  inc.1 <- inc.allele[,pair[1]]; inc.2 <- inc.allele[, pair[2]]
  union.inc <- which((inc.1 & obs.allele[,pair[2]]) | (inc.2 & obs.allele[,pair[1]]))
  if (sum(as.numeric(union.inc)) < PARAM$min.shared_obs) {return(NA)}
  shared.inc <- which(inc.1[union.inc] & inc.2[union.inc])
  length(shared.inc) / length(union.inc)
}
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

#Downsample to maximum number of alleles
if (nrow(mat.allele) > PARAM$max.alleles) {mat.allele <- mat.allele[sample(nrow(mat.allele), size=PARAM$max.alleles, replace=F), ]}

#Generate allele observation and incidence matrices
obs.allele <- mat.allele > 0
inc.allele <- mat.allele > 10^-12
rm(mat.allele); invisible(gc())
###################################
#Pre-calculate all pairwise allele profile similarities within cohorts, within body sites
collect.allele_overlap <- list()
results.allele_overlap <- data.delta_t
data.sample$cohort.bg <- data.sample$cohort; data.sample$cohort.bg[data.sample$cohort.bg == "MM"] <- "FR-CRC"
data.timepoint$cohort.bg <- data.timepoint$cohort; data.timepoint$cohort.bg[data.timepoint$cohort.bg == "MM"] <- "FR-CRC"
#Iterate through cohorts with longitudinal data
for (coh in c("FR-CRC", "Zhang-RA", "LU-T1D")) {
  #Get unique combinations of subjects and timepoints for current cohort
  combn.s_t <- t(combn(rownames(data.timepoint)[data.timepoint$cohort.bg == coh], 2))
  combn.samples <- t(combn(rownames(data.sample)[data.sample$cohort.bg == coh], 2))
  
  #Annotate current sample combinations to habitats
  combn.samples.habitats <- cbind(as.character(data.sample[combn.samples[,1], "material"]), as.character(data.sample[combn.samples[,2], "material"]))
  combn.samples.subject <- cbind(as.character(data.sample[combn.samples[,1], "subject"]), as.character(data.sample[combn.samples[,2], "subject"]))
  combn.samples.s_t <- cbind(as.character(data.sample[combn.samples[,1], "subject_timepoint"]), as.character(data.sample[combn.samples[,2], "subject_timepoint"]))
  combn.samples.family <- cbind(as.character(data.sample[combn.samples[,1], "family"]), as.character(data.sample[combn.samples[,2], "family"]))
  
  #Preallocate
  collect.allele_overlap[[coh]] <- list()
  
  #Loop through body sites and get pairwise allele frequency overlaps
  for (hab in c("saliva", "dental_plaque", "stool")) {
    #Select current list of within-site comparisons
    #=> within site, but not within family
    curr.combn <- combn.samples.habitats[,1] == hab & combn.samples.habitats[,2] == hab & combn.samples.family[,1] != combn.samples.family[,2]
    
    #Calculate pairwise allele profile agreement
    if (sum(curr.combn) > 2) {
      if (sum(curr.combn) > PARAM$max.combn) {
        curr.combn.choose <- sample(which(curr.combn), PARAM$max.combn, replace=T)
        curr.overlap <- apply(combn.samples[curr.combn.choose, ], 1, get.allele_overlap)
      } else {
        curr.combn.choose <- which(curr.combn)
        curr.overlap <- apply(combn.samples[curr.combn, ], 1, get.allele_overlap)
      }
      collect.allele_overlap[[coh]][[hab]] <- curr.overlap
      
      # curr.overlap <- apply(combn.samples[curr.combn, ], 1, get.allele_overlap)
      # curr.overlap.long <- rep.int(NA, nrow(combn.samples))
      # curr.overlap.long[curr.combn] <- curr.overlap
      # curr.overlap.mat <- matrix(NA, nrow=sum(data.timepoint$cohort.bg == coh), ncol=sum(data.timepoint$cohort.bg == coh), dimnames=list(rownames(data.timepoint)[data.timepoint$cohort.bg == coh], rownames(data.timepoint)[data.timepoint$cohort.bg == coh]))
      # curr.overlap.mat[cbind(combn.samples.s_t[curr.combn, 1], combn.samples.s_t[curr.combn, 2])] <- curr.overlap
      # diag(curr.overlap.mat) <- 1
      # collect.allele_overlap[[coh]][[hab]] <- curr.overlap.mat
    } else {
      collect.allele_overlap[[coh]][[hab]] <- NA
      next()
    }
    
    #Get current generic background overlap (which is basically what was calculated above)
    curr.overlap.bg.generic <- collect.allele_overlap[[coh]][[hab]]
    
    #Iterate through subjects and timepoint-combinations and store current overlaps
    curr.subjects <- unique(c(combn.samples.subject))[unique(c(combn.samples.subject)) %in% as.character(data.delta_t$subject)]
    for (subj in curr.subjects) {
      #Get all timepoint combinations for current subject
      if (length(unique(data.sample[data.sample$subject == subj, "timepoint"])) < 2) {next()}
      curr.tp <- t(combn(unique(data.sample[data.sample$subject == subj, "timepoint"]), 2))
      for (tt in 1:nrow(curr.tp)) {
        #Get current overlap
        curr.samples.within <- rownames(data.sample)[data.sample$subject_timepoint %in% paste(subj, curr.tp[tt,], sep=".") & data.sample$material == hab]
        #SKip if data is not available for both current subject_timepoints
        if (length(curr.samples.within) < 2) {next()}
        curr.overlap.within <- get.allele_overlap(curr.samples.within)
        #Get current subject-specific background
        # curr.overlap.bg.subj <- c(
        #   curr.overlap.mat[paste(subj, curr.tp[tt,1], sep="."), data.timepoint[colnames(curr.overlap.mat), "subject"] != subj],
        #   curr.overlap.mat[paste(subj, curr.tp[tt,2], sep="."), data.timepoint[colnames(curr.overlap.mat), "subject"] != subj]
        # )
        #Get index of current pair in delta_t
        curr.idx.delta_t <- which(data.delta_t$site.1 == hab & data.delta_t$site.2 == hab & data.delta_t$subject == subj & ((data.delta_t$timepoint.1 == curr.tp[tt,1] & data.delta_t$timepoint.2 == curr.tp[tt,2]) | (data.delta_t$timepoint.1 == curr.tp[tt,2] & data.delta_t$timepoint.2 == curr.tp[tt,1])))
        if (length(curr.idx.delta_t) != 1) {next()}
        #Store
        results.allele_overlap[curr.idx.delta_t, paste(hab, "allele.overlap", sep=".")] <- curr.overlap.within
        results.allele_overlap[curr.idx.delta_t, paste(hab, "allele.overlap.bg_generic.mean", sep=".")] <- mean(curr.overlap.bg.generic, na.rm=T)
        results.allele_overlap[curr.idx.delta_t, paste(hab, "allele.overlap.bg_generic.sd", sep=".")] <- sd(curr.overlap.bg.generic, na.rm=T)
        #results.allele_overlap[curr.idx.delta_t, paste(hab, "allele.overlap.bg_subject.mean", sep=".")] <- mean(curr.overlap.bg.subj, na.rm=T)
        #results.allele_overlap[curr.idx.delta_t, paste(hab, "allele.overlap.bg_generic.sd", sep=".")] <- sd(curr.overlap.bg.subj, na.rm=T)
      }
    }
    
    #Report
    writeLines(paste(date(), "=> Completed", tax.id, taxon, coh, hab))
  }
}
###################################


###################################
#Store data
save(results.allele_overlap, collect.allele_overlap, file=paste0(PARAM$folder.results, "longitudinal.stability/", tax.id, ".", print.taxon,  ".longitudinal_allele_profiles.RData"))

#Report
writeLines(paste(date(), "=> Done with", tax.id, taxon))
################################################################################
################################################################################

q(save="no")




