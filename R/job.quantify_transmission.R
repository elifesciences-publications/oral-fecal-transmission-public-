#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Quantify putative oral-gut transmission from SNV data.
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load SNV tables per taxon
#=> identify "personal SNVs" per subject
#=> for each individual, track (asymmetric) overlap of personal SNVs btw. oral & gut
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
PARAM$folder.base <- "##################"
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.logs <- paste0(PARAM$folder.base, "logs/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Set parameters
#PARAM$min.shared_obs <- 20

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

###################################
#Preallocate per-individual results collector
data.s_t.transmission <- data.frame()

#Iterate through subject_timepoints
for (s_t in unique(as.character(subjects.timepoints))) {
  #Get focal samples for current subject_timepoint
  curr.sample.saliva <- rownames(data.sample)[subjects.timepoints == s_t & habitats == "saliva"]
  curr.sample.plaque <- rownames(data.sample)[subjects.timepoints == s_t & habitats == "plaque"]
  curr.sample.stool <- rownames(data.sample)[subjects.timepoints == s_t & habitats == "stool"]
  
  #Skip if only one body site has data
  if (sum(c(length(curr.sample.saliva), length(curr.sample.plaque), length(curr.sample.stool))) <= 1) {next()}
  
  #Get current subject
  curr.subject <- unique(as.character(subjects[subjects.timepoints == s_t]))
  
  #Identify "personal SNVs" for current subject_timepoint
  #=> exclude other timepoints for same subject from background (obviously)
  curr.comb_alleles.s_t <- which(rowSums(inc.allele[, na.omit(c(curr.sample.saliva, curr.sample.plaque, curr.sample.stool))]) > 0)
  curr.comb_alleles.background <- which(rowSums(inc.allele[, subjects != curr.subject]) > 0)
  curr.personal_snv <- setdiff(curr.comb_alleles.s_t, curr.comb_alleles.background)
  
  #Report
  writeLines(paste(date(), "=>", print.taxon, "has", length(curr.personal_snv), "personal SNVs in", s_t))
  
  #Collect per-body site data
  #=> saliva
  if (length(curr.sample.saliva) == 0) {
    curr.sample.saliva <- curr.specI.saliva <- curr.hor_cov.saliva <- curr.ver_cov.saliva <- NA
    curr.obs.saliva <- curr.inv.saliva <- F
    curr.inc.saliva <- curr.obs.saliva <- curr.personal_snv.inc.saliva <- curr.personal_snv.obs.saliva <- numeric()
  } else {
    #specI abundance & coverage
    curr.specI.saliva <- data.specI.rel[tax.id, curr.sample.saliva]
    curr.hor_cov.saliva <- data.hor_cov[tax.id, curr.sample.saliva]
    curr.ver_cov.saliva <- data.ver_cov[tax.id, curr.sample.saliva]
    #Global allele incidence and observed position count
    curr.obs.saliva <- obs.allele[, curr.sample.saliva]
    curr.inc.saliva <- inc.allele[, curr.sample.saliva]
    #Personal SNVs in body site
    curr.personal_snv.inc.saliva <- curr.inc.saliva[curr.personal_snv]
    curr.personal_snv.obs.saliva <- curr.obs.saliva[curr.personal_snv]
  }
  #=> plaque
  if (length(curr.sample.plaque) == 0) {
    curr.sample.plaque <- curr.specI.plaque <- curr.hor_cov.plaque <- curr.ver_cov.plaque <- NA
    curr.obs.plaque <- curr.inv.plaque <- F
    curr.inc.plaque <- curr.obs.plaque <- curr.personal_snv.inc.plaque <- curr.personal_snv.obs.plaque <- numeric()
  } else {
    curr.specI.plaque <- data.specI.rel[tax.id, curr.sample.plaque]
    curr.hor_cov.plaque <- data.hor_cov[tax.id, curr.sample.plaque]
    curr.ver_cov.plaque <- data.ver_cov[tax.id, curr.sample.plaque]
    curr.inc.plaque <- inc.allele[, curr.sample.plaque]
    curr.obs.plaque <- obs.allele[, curr.sample.plaque]
    curr.personal_snv.inc.plaque <- curr.inc.plaque[curr.personal_snv]
    curr.personal_snv.obs.plaque <- curr.obs.plaque[curr.personal_snv]
  }
  #=> stool
  if (length(curr.sample.stool) == 0) {
    curr.sample.stool <- curr.specI.stool <- curr.hor_cov.stool <- curr.ver_cov.stool <- NA
    curr.obs.stool <- curr.inv.stool <- F
    curr.inc.stool <- curr.obs.stool <- curr.personal_snv.inc.stool <- curr.personal_snv.obs.stool <- numeric()
  } else {
    curr.specI.stool <- data.specI.rel[tax.id, curr.sample.stool]
    curr.hor_cov.stool <- data.hor_cov[tax.id, curr.sample.stool]
    curr.ver_cov.stool <- data.ver_cov[tax.id, curr.sample.stool]
    curr.inc.stool <- inc.allele[, curr.sample.stool]
    curr.obs.stool <- obs.allele[, curr.sample.stool]
    curr.personal_snv.inc.stool <- curr.inc.stool[curr.personal_snv]
    curr.personal_snv.obs.stool <- curr.obs.stool[curr.personal_snv]
  }
  
  
  #Quantify saliva <-> stool
  if (!is.na(curr.sample.saliva) & !is.na(curr.sample.stool)) {
    #Which alleles have observations in both oral and gut?
    curr.obs.shared.saliva_stool <- intersect(which(curr.obs.saliva), which(curr.obs.stool))
    curr.obs.total.saliva_stool <- union(which(curr.obs.saliva), which(curr.obs.stool))
    
    #Reduce to shared observations
    curr.shared_obs.inc.saliva <- curr.inc.saliva[curr.obs.shared.saliva_stool]
    curr.shared_obs.inc.stool <- curr.inc.stool[curr.obs.shared.saliva_stool]
    #Quantify agreement between oral and gut on shared pos
    curr.shared.double_zeroes.saliva_stool <- !curr.shared_obs.inc.saliva & !curr.shared_obs.inc.stool
    curr.shared_alleles.saliva_stool <- which(curr.shared_obs.inc.saliva == curr.shared_obs.inc.stool & !curr.shared.double_zeroes.saliva_stool)
    
    #How many personal SNVs are shared between the oral and gut samples?
    curr.personal_snv.shared.saliva_stool <- sum(curr.personal_snv.inc.saliva & curr.personal_snv.inc.stool)
    
    #Get some summary stats
    n.aSNV.obs.shared.saliva_stool <- length(curr.obs.shared.saliva_stool)
    n.aSNV.obs.total.saliva_stool <- length(curr.obs.total.saliva_stool)
  } else {
    curr.shared_obs.inc.saliva <- curr.shared.double_zeroes.saliva_stool <- curr.shared_alleles.saliva_stool <- curr.shared.double_zeroes.saliva_stool <- F 
    curr.obs.shared.saliva_stool <- curr.obs.total.saliva_stool <- curr.personal_snv.shared.saliva_stool <- n.aSNV.obs.shared.saliva_stool <- n.aSNV.obs.total.saliva_stool <- NA
  }
  #Quantify plaque <-> stool
  if (!is.na(curr.sample.plaque) & !is.na(curr.sample.stool)) {
    #Which alleles have observations in both oral and gut?
    curr.obs.shared.plaque_stool <- intersect(which(curr.obs.plaque), which(curr.obs.stool))
    curr.obs.total.plaque_stool <- union(which(curr.obs.plaque), which(curr.obs.stool))
    
    #Reduce to shared observations
    curr.shared_obs.inc.plaque <- curr.inc.plaque[curr.obs.shared.plaque_stool]
    curr.shared_obs.inc.stool <- curr.inc.stool[curr.obs.shared.plaque_stool]
    #Quantify agreement between oral and gut on shared pos
    curr.shared.double_zeroes.plaque_stool <- !curr.shared_obs.inc.plaque & !curr.shared_obs.inc.stool
    curr.shared_alleles.plaque_stool <- which(curr.shared_obs.inc.plaque == curr.shared_obs.inc.stool & !curr.shared.double_zeroes.plaque_stool)
    
    #How many personal SNVs are shared between the oral and gut samples?
    curr.personal_snv.shared.plaque_stool <- sum(curr.personal_snv.inc.plaque & curr.personal_snv.inc.stool)
    
    #Get some summary stats
    n.aSNV.obs.shared.plaque_stool <- length(curr.obs.shared.plaque_stool)
    n.aSNV.obs.total.plaque_stool <- length(curr.obs.total.plaque_stool)
  } else {
    curr.shared_obs.inc.plaque <- curr.shared.double_zeroes.plaque_stool <- curr.shared_alleles.plaque_stool <- curr.shared.double_zeroes.plaque_stool <- F 
    curr.obs.shared.plaque_stool <- curr.obs.total.plaque_stool <- curr.personal_snv.shared.plaque_stool <- n.aSNV.obs.shared.plaque_stool <- n.aSNV.obs.total.plaque_stool <- NA
  }
  #Quantify saliva <-> plaque
  if (!is.na(curr.sample.saliva) & !is.na(curr.sample.plaque)) {
    #Which alleles have observations in both oral and gut?
    curr.obs.shared.saliva_plaque <- intersect(which(curr.obs.saliva), which(curr.obs.plaque))
    curr.obs.total.saliva_plaque <- union(which(curr.obs.saliva), which(curr.obs.plaque))
    
    #Reduce to shared observations
    curr.shared_obs.inc.saliva <- curr.inc.saliva[curr.obs.shared.saliva_plaque]
    curr.shared_obs.inc.plaque <- curr.inc.plaque[curr.obs.shared.saliva_plaque]
    #Quantify agreement between oral and gut on shared pos
    curr.shared.double_zeroes.saliva_plaque <- !curr.shared_obs.inc.saliva & !curr.shared_obs.inc.plaque
    curr.shared_alleles.saliva_plaque <- which(curr.shared_obs.inc.saliva == curr.shared_obs.inc.plaque & !curr.shared.double_zeroes.saliva_plaque)
    
    #How many personal SNVs are shared between the oral and gut samples?
    curr.personal_snv.shared.saliva_plaque <- sum(curr.personal_snv.inc.saliva & curr.personal_snv.inc.plaque)
    
    #Get some summary stats
    n.aSNV.obs.shared.saliva_plaque <- length(curr.obs.shared.saliva_plaque)
    n.aSNV.obs.total.saliva_plaque <- length(curr.obs.total.saliva_plaque)
  } else {
    curr.shared_obs.inc.saliva <- curr.shared.double_zeroes.saliva_plaque <- curr.shared_alleles.saliva_plaque <- curr.shared.double_zeroes.saliva_plaque <- F 
    curr.obs.shared.saliva_plaque <- curr.obs.total.saliva_plaque <- curr.personal_snv.shared.saliva_plaque <- n.aSNV.obs.shared.saliva_plaque <- n.aSNV.obs.total.saliva_plaque <- NA
  }
  
  #Prepare data to pass to store results
  data.s_t.transmission <- rbind(data.s_t.transmission, data.frame(
    Tax_ID = tax.id,
    Scientific_Name = taxon,
    subject_timepoint = s_t,
    sample.saliva = curr.sample.saliva,
    sample.plaque = curr.sample.plaque,
    sample.stool = curr.sample.stool,
    #specI rel. abundances in oral and gut
    specI.saliva = curr.specI.saliva,
    specI.plaque = curr.specI.plaque,
    specI.stool = curr.specI.stool,
    #Horizontal and vertical coverage in oral and gut
    hor_cov.saliva = curr.hor_cov.saliva,
    hor_cov.plaque = curr.hor_cov.plaque,
    hor_cov.stool = curr.hor_cov.stool,
    ver_cov.saliva = curr.ver_cov.saliva,
    ver_cov.plaque = curr.ver_cov.plaque,
    ver_cov.stool = curr.ver_cov.stool,
    #Total observed alleles ("aSNV" for "all SNVs") in oral and gut
    n.aSNV.obs.saliva = sum(curr.obs.saliva),
    n.aSNV.obs.plaque = sum(curr.obs.plaque),
    n.aSNV.obs.stool = sum(curr.obs.stool),
    #Shared and observed alleles btw. oral and gut
    n.aSNV.obs.shared.saliva_stool = n.aSNV.obs.shared.saliva_stool,
    n.aSNV.obs.total.saliva_stool = n.aSNV.obs.total.saliva_stool,
    n.aSNV.obs.shared.plaque_stool = n.aSNV.obs.shared.plaque_stool,
    n.aSNV.obs.total.plaque_stool = n.aSNV.obs.total.plaque_stool,
    n.aSNV.obs.shared.saliva_plaque = n.aSNV.obs.shared.saliva_plaque,
    n.aSNV.obs.total.saliva_plaque = n.aSNV.obs.total.saliva_plaque,
    #Total allele count (incidence) in oral and gut
    n.aSNV.inc.saliva = sum(curr.inc.saliva),
    n.aSNV.inc.plaque = sum(curr.inc.plaque),
    n.aSNV.inc.stool = sum(curr.inc.stool),
    #Overlap of oral and gut allele profiles on positions with shared obs
    n.aSNV.shared.saliva_stool = length(curr.shared_alleles.saliva_stool),
    f.aSNV.shared.saliva_stool = length(curr.shared_alleles.saliva_stool) / (length(curr.obs.shared.saliva_stool) - sum(curr.shared.double_zeroes.saliva_stool)),
    n.aSNV.shared.plaque_stool = length(curr.shared_alleles.plaque_stool),
    f.aSNV.shared.plaque_stool = length(curr.shared_alleles.plaque_stool) / (length(curr.obs.shared.plaque_stool) - sum(curr.shared.double_zeroes.plaque_stool)),
    n.aSNV.shared.saliva_stool = length(curr.shared_alleles.saliva_stool),
    f.aSNV.shared.saliva_plaque = length(curr.shared_alleles.saliva_plaque) / (length(curr.obs.shared.saliva_plaque) - sum(curr.shared.double_zeroes.saliva_plaque)),
    #Personal SNV count
    n.personal_SNV = length(curr.personal_snv),
    #Personal SNV count in saliva
    n.pSNV.saliva = sum(curr.personal_snv.inc.saliva),
    #Personal SNV count in plaque
    n.pSNV.plaque = sum(curr.personal_snv.inc.plaque),
    #Personal SNV count in stool
    n.pSNV.stool = sum(curr.personal_snv.inc.stool),
    #Number of personal alleles observed in saliva
    n.pSNV.obs.saliva = sum(curr.personal_snv.obs.saliva),
    #Number of personal alleles observed in plaque
    n.pSNV.obs.plaque = sum(curr.personal_snv.obs.plaque),
    #Number of personal alleles observed in stool
    n.pSNV.obs.stool = sum(curr.personal_snv.obs.stool),
    #Number of shared personal SNVs
    n.pSNV.shared.saliva_stool = curr.personal_snv.shared.saliva_stool,
    n.pSNV.shared.plaque_stool = curr.personal_snv.shared.plaque_stool,
    n.pSNV.shared.saliva_plaque = curr.personal_snv.shared.saliva_plaque
  ))
}

#Store data
save(data.s_t.transmission, file=paste0(PARAM$folder.results, "quantify.transmission/", tax.id, ".", print.taxon,  ".transmission_per_subject.RData"))

#Report
writeLines(paste(date(), "=> Done with", tax.id, taxon))
################################################################################
################################################################################

q(save="no")




