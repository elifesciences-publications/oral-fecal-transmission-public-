#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Prepare allele data to detect putative transmission between two sites.
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load allele tables per taxon


#=> save and export individual datafiles in R format
#
#2017-10-17
#sebastian.schmidt@embl.de
################################################################################

#
#NOTE: This script was generated based on the "template" of a previous script to
#detect transmission between saliva and stool. Because of this, many variable
#names will be counter-intuitive: they were "x.oral" and "x.gut" before, where
#oral would be "saliva" and gut would be "stool". In the present script, names
#are remapped so that "*.gut" variables refer to saliva samples, and "*.oral" to
#plaque samples. Moreover, there is only one cohort in the dataset that contains
#plaque samples: Zhang-RA. Therefore, only within-cohort comparisons for this
#particular cohort are used.
#

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
PARAM$folder.base <- "###################"
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.logs <- paste0(PARAM$folder.base, "logs/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Set parameters
#=> minimum required shared alleles for comparisons
PARAM$min.shared_obs <- 20
#=> maximum number of background comparisons (beyond which distributions are sampled)
PARAM$max.background_size <- 10000

#Read current taxon of interest from command line
args <- commandArgs(TRUE)
tax.id <- args[2]
site_1 <- args[3]
site_2 <- args[4]

if (! site_1 %in% c("saliva", "stool", "dental_plaque")) {writeLines(paste(date(), "=> Unknown site_1:", tax.id, site_1)); q("no")}
if (! site_2 %in% c("saliva", "stool", "dental_plaque")) {writeLines(paste(date(), "=> Unknown site_2:", tax.id, site_2)); q("no")}

#Legacy renaming of dental_plaque -> plaque, saliva -> oral and stool -> gut
if (site_1 == "dental_plaque") {site_1.short <- "plaque"} else if (site_1 == "stool") {site_1.short <- "gut"} else {site_1.short <- site_1}
if (site_2 == "dental_plaque") {site_2.short <- "plaque"} else if (site_2 == "stool") {site_2.short <- "gut"} else {site_2.short <- site_2}

#Report
writeLines(paste(date(), "=> Processing", tax.id, site_1, site_2))
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
#Filter down to include only positions which are a SNV in ≥1 site_1 and ≥1 site_2 sample
if ("dental_plaque" %in% c(site_1, site_2)) {
  mat.allele <- mat.allele[data.allele[, paste0("covered.", site_1.short)] & data.allele[, paste0("covered.", site_2.short)], data.sample[colnames(mat.allele), "material"] %in% c(site_1, site_2) & data.sample[colnames(mat.allele), "cohort"] == "Zhang-RA"]
} else {
  mat.allele <- mat.allele[data.allele[, paste0("covered.", site_1.short)] & data.allele[, paste0("covered.", site_2.short)], data.sample[colnames(mat.allele), "material"] %in% c(site_1, site_2)]
}
rm(data.allele); invisible(gc())

#Get current habitat vector
habitats <- data.sample[colnames(mat.allele), "material"]
subjects <- data.sample[colnames(mat.allele), "subject"]
subjects.timepoints <- data.sample[colnames(mat.allele), "subject_timepoint"]

#Get current taxon name
taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
print.taxon <- gsub(" ", "_", taxon)
print.taxon <- gsub("/", "_", print.taxon)

###################################
#Identify positions of interest
#=> minimum number of observations in both saliva and gut (≥10 samples each)
#=> allele in at least one oral and one gut sample (filtered already above)
obs.allele <- mat.allele > 0

#Get an allele incidence matrix
#=> TRUE if allele was observed, FALSE if not (ref or "-1")
inc.allele <- mat.allele > 10^-12

#Exit if positions are distributed across habitats too unevenly
#if (nrow(obs.oral) < 20) {writeLines(paste(date(), "=> Too few oral positions:", tax.id)); q(save="no")}
#if (nrow(obs.gut) < 20) {writeLines(paste(date(), "=> Too few gut positions:", tax.id)); q(save="no")}

#Get raw observation frequencies per habitat
obs.site_1.sum <- rowSums(obs.allele[, habitats == site_1])
obs.site_2.sum <- rowSums(obs.allele[, habitats == site_2])

#Get raw incidence frequencies per habitat
inc.site_1.sum <- rowSums(inc.allele[, habitats == site_1])
inc.site_2.sum <- rowSums(inc.allele[, habitats == site_2])

#Require that an allele...
#=> is observed in at least one oral and one gut sample (otherwise, would be trivial)
#=> is not observed in *all* oral or *all* gut samples (otherwise, would be trivial, as well)
#=> has position coverage in at least 10 oral and 10 gut samples
keep.alleles <- obs.site_1.sum >= 10 & obs.site_2.sum >= 10 & inc.site_1.sum >= 1 & inc.site_2.sum >= 1 & inc.site_1.sum < sum(habitats == site_1) & inc.site_2.sum < sum(habitats == site_2)

#Exit if too few alleles remain
if (sum(keep.alleles) < 100) {writeLines(paste(date(), "=> Too few observed alleles:", tax.id)); q(save="no")}

#Filter down to relevant positions
mat.allele <- mat.allele[keep.alleles, ]
invisible(gc())
###################################

###################################
#Re-turn mat.allele into an incidence matrix (after filtering)
#=> TRUE if position was observed (as current allele or other), FALSE if not ("-1")
obs.allele <- mat.allele > 0
obs.site_1 <- obs.allele[, habitats == site_1]
obs.site_2 <- obs.allele[, habitats == site_2]

#Get an allele incidence matrix
#=> TRUE if allele was observed, FALSE if not (ref or "-1")
inc.allele <- mat.allele > 10^-12
inc.site_1 <- inc.allele[, habitats == site_1]
inc.site_2 <- inc.allele[, habitats == site_2]

#Tidy up
rm(mat.allele); invisible(gc())

#Get current total allele count
n.allele <- nrow(obs.allele)

#Get current habitat vector
habitats <- data.sample[colnames(obs.allele), "material"]
subjects <- data.sample[colnames(obs.allele), "subject"]
subjects.timepoints <- data.sample[colnames(obs.allele), "subject_timepoint"]

#Get the number of remaining eligible subject.timepoints with paired site_1-site_2 observations
curr.sample_pairs <- as.matrix(expand.grid(colnames(obs.site_1), colnames(obs.site_2)))
curr.pairs.st <- cbind(as.character(data.sample[curr.sample_pairs[,1], "subject_timepoint"]), as.character(data.sample[curr.sample_pairs[,2], "subject_timepoint"]))
curr.pairs.intra <- curr.sample_pairs[curr.pairs.st[,1] == curr.pairs.st[,2], ]

#Skip current taxon if there are no paired subject.timepoints left
if (nrow(curr.pairs.intra) < 2) {writeLines(paste(date(), tax.id, taxon, "=> Too few intra-individual paired observations.")); q("no")}

#Get list of subjects and subject_timepoints with paired intra observations
curr.pairs.subject <- as.character(data.sample[curr.pairs.intra[,1], "subject"])
curr.pairs.subject_timepoint <- as.character(data.sample[curr.pairs.intra[,1], "subject_timepoint"])
###################################

###################################
#Collapse timepoints in allele observation and incidence tables
#=> combine information from multiple timepoints of an individual
#=> retain any allele that was *ever* observed for a given subject
#=> the updated tables are alleles x subject, where subjects are collapsed longitudinally
#=> additionally, populate a data.subject frame for metadata
###################################
#Preallocate
data.subject <- data.frame()
subjects.with_site_1 <- unique(as.character(subjects[habitats == site_1]))
subjects.with_site_2 <- unique(as.character(subjects[habitats == site_2]))
obs.allele.subject.site_1 <- inc.allele.subject.site_1 <- matrix(F, nrow=n.allele, ncol=length(subjects.with_site_1), dimnames=list(rownames(obs.allele), subjects.with_site_1))
obs.allele.subject.site_2 <- inc.allele.subject.site_2 <- matrix(F, nrow=n.allele, ncol=length(subjects.with_site_2), dimnames=list(rownames(obs.allele), subjects.with_site_2))

#Iterate through subjects
for (subj in union(subjects.with_site_1, subjects.with_site_2)) {
  #Collapse oral samples
  if (subj %in% subjects.with_site_1) {
    obs.allele.subject.site_1[, subj] <- Reduce("|", as.data.frame(as.matrix(obs.allele[, subjects == subj & habitats == site_1])))
    inc.allele.subject.site_1[, subj] <- Reduce("|", as.data.frame(as.matrix(inc.allele[, subjects == subj & habitats == site_1])))
    available.site_1 <- TRUE
  } else {available.site_1 <- FALSE}
  
  #Collapse gut samples
  if (subj %in% subjects.with_site_2) {
    obs.allele.subject.site_2[, subj] <- Reduce("|", as.data.frame(as.matrix(obs.allele[, subjects == subj & habitats == site_2])))
    inc.allele.subject.site_2[, subj] <- Reduce("|", as.data.frame(as.matrix(inc.allele[, subjects == subj & habitats == site_2])))
    available.site_2 <- TRUE
  } else {available.site_2 <- FALSE}
  
  #Collect info on current subject
  data.subject <- rbind(data.subject, data.frame(
    subject = subj,
    available.site_1 = available.site_1,
    available.site_2 = available.site_2,
    n.tp = length(unique(data.sample$timepoint[data.sample$subject == subj])),
    cohort = data.sample$cohort[data.sample$subject == subj][1],
    family = data.sample$family[data.sample$subject == subj][1]
  ))
}

#Add an imputed cohort label for MM subjects (to pool them with FR-CRC)
data.subject$cohort.bg <- data.subject$cohort; data.subject$cohort.bg[data.subject$cohort == "MM"] <- "FR-CRC"
rownames(data.subject) <- data.subject$subject
###################################

###################################
#Further filter alleles to exlucde trivial cases
#=> do not consider alleles which were observed in *all* oral or *all* gut samples
keep.alleles <- (rowSums(inc.allele.subject.site_1) < rowSums(obs.allele.subject.site_1)) & (rowSums(inc.allele.subject.site_2) < rowSums(obs.allele.subject.site_2))
if (sum(keep.alleles) < 100) {writeLines(paste(date(), tax.id, taxon, "=> Too few alleles remaining after filtering.")); q(save="no")}

#Apply filter
obs.allele.subject.site_1 <- Matrix(obs.allele.subject.site_1[keep.alleles,], sparse=T)
inc.allele.subject.site_1 <- Matrix(inc.allele.subject.site_1[keep.alleles,], sparse=T)
obs.allele.subject.site_2 <- Matrix(obs.allele.subject.site_2[keep.alleles,], sparse=T)
inc.allele.subject.site_2 <- Matrix(inc.allele.subject.site_2[keep.alleles,], sparse=T)
obs.allele <- Matrix(obs.allele[keep.alleles,], sparse=T)
inc.allele <- Matrix(inc.allele[keep.alleles,], sparse=T)
###################################

###################################
#Store
save(
  obs.allele,
  inc.allele,
  obs.allele.subject.site_1,
  obs.allele.subject.site_2,
  inc.allele.subject.site_1,
  inc.allele.subject.site_2,
  curr.sample_pairs,
  curr.pairs.st,
  curr.pairs.intra,
  curr.pairs.subject,
  curr.pairs.subject_timepoint,
  data.subject,
  subjects.with_site_1,
  subjects.with_site_2,
  file=paste0(PARAM$folder.SNV_data, tax.id, ".allele_data.filtered.subset.", site_1.short, "_", site_2.short, ".RData")
)
################################################################################
################################################################################



q(save="no")


