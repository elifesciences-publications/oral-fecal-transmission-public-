#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Detect putative oral-gut transmission from SNV data.
#
#
#=> for one taxon at a time (read from command line input)
#=> load sample metadata
#=> load SNV tables per taxon
#=> calculate SNV frequencies in oral and gut samples, per position
#=> for each individual, calculate deviation of observed SNV patterns from bg
#=> save and export individual datafiles in R format
#
#2017-10-17
#sebastian.schmidt@embl.de
################################################################################

#
#NOTE: This script is an update to previous, site-specific scripts (where
#"saliva", "dental_plaque" and "stool" transmission detection was hard-coded).
#Moreover, this script previously included memory-intense early processing that
#has now been outsourced to the script
#"job.prepare.detect_transmission.alleles.R".
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
PARAM$folder.base <- "###########"
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.logs <- paste0(PARAM$folder.base, "logs/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Set parameters
#=> minimum required shared alleles for comparisons
PARAM$min.shared_obs <- 20
#=> maximum number of background comparisons (beyond which distributions are sampled)
PARAM$max.background_size <- 1000

#Read current taxon of interest from command line
args <- commandArgs(TRUE)
tax.id <- args[2]
site_1 <- args[3]
site_2 <- args[4]

if (! site_1 %in% c("saliva", "stool", "dental_plaque")) {writeLines(paste(date(), "=> Unknown site_1:", tax.id, site_1)); q("no")}
if (! site_2 %in% c("saliva", "stool", "dental_plaque")) {writeLines(paste(date(), "=> Unknown site_2:", tax.id, site_2)); q("no")}

#Legacy renaming of dental_plaque -> plaque and saliva -> oral
if (site_1 == "dental_plaque") {site_1.short <- "plaque"} else if (site_1 == "stool") {site_1.short <- "gut"} else {site_1.short <- site_1}
if (site_2 == "dental_plaque") {site_2.short <- "plaque"} else if (site_2 == "stool") {site_2.short <- "gut"} else {site_2.short <- site_2}

#Set switch to decide whether full backgrounds need to be calculated or not
#=> comparisons with plaque don't need this!
if ("dental_plaque" %in% c(site_1, site_2)) {calculate.full <- FALSE} else {calculate.full <- TRUE}

#Report
writeLines(paste(date(), "=> Processing", tax.id))
################################################################################
################################################################################


################################################################################
################################################################################
#Load allele data
file.allele_data <- paste0(PARAM$folder.SNV_data, tax.id, ".allele_data.filtered.subset.", site_1.short, "_", site_2.short, ".RData")
if (! file.exists(file.allele_data)) {writeLines(paste(tax.id, "file does not exist:", file.allele_data)); q("no")}
load(file.allele_data)

#Load sample data
load(paste0(PARAM$folder.parameters, "data.sample.RData"))
samples <- colnames(obs.allele)
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

#Get current habitat vector
habitats <- data.sample[colnames(obs.allele), "material"]
subjects <- data.sample[colnames(obs.allele), "subject"]
subjects.timepoints <- data.sample[colnames(obs.allele), "subject_timepoint"]

#Get current taxon name
taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
print.taxon <- gsub(" ", "_", taxon)
print.taxon <- gsub("/", "_", print.taxon)
################################################################################
################################################################################


################################################################################
################################################################################
###################################
#Calculate subject-specific allele frequency backgrounds
#=> per cohort and across cohorts
#=> block for family (exclude family/village members for subject-spec background)
###################################
#Preallocate
n.allele <- nrow(obs.allele)
if (calculate.full) {
  bg.allele.subject.site_1.full <- Matrix(matrix(0, nrow=n.allele, ncol=ncol(obs.allele.subject.site_1), dimnames=list(rownames(obs.allele), subjects.with_site_1)), sparse=T)
  bg.allele.subject.site_2.full <- Matrix(matrix(0, nrow=n.allele, ncol=ncol(obs.allele.subject.site_2), dimnames=list(rownames(obs.allele), subjects.with_site_2)), sparse=T)
}
bg.allele.subject.site_1.cohort <- Matrix(matrix(0, nrow=n.allele, ncol=ncol(obs.allele.subject.site_1), dimnames=list(rownames(obs.allele), subjects.with_site_1)), sparse=T)
bg.allele.subject.site_2.cohort <- Matrix(matrix(0, nrow=n.allele, ncol=ncol(obs.allele.subject.site_2), dimnames=list(rownames(obs.allele), subjects.with_site_2)), sparse=T)

#Iterate through subjects
for (subj in union(subjects.with_site_1, subjects.with_site_2)) {
  #Get site_1 background
  if (subj %in% subjects.with_site_1) {
    #Get indices of subjects to take into account as part of the global background (this includes the subject themselves, but not their families!)
    if (calculate.full) {
      idx.full <- sort(c(which(data.subject$family[data.subject$available.site_1] != data.subject[subj, "family"]), which(rownames(data.subject)[data.subject$available.site_1] == subj)))
      #Require a minimum background size
      if (length(idx.full) > 10) {
        #Calculate allele frequencies against this global background
        bg.allele.subject.site_1.full[, subj] <- rowSums(inc.allele.subject.site_1[, idx.full]) / rowSums(obs.allele.subject.site_1[, idx.full])
      }
    }
    
    #Get indices of subjects to take into account as part of the cohort-specific background (this includes the subject themselves, but not their families!)
    idx.cohort <- sort(c(which(data.subject$cohort.bg[data.subject$available.site_1] == data.subject[subj, "cohort.bg"] & data.subject$family[data.subject$available.site_1] != data.subject[subj, "family"]), which(rownames(data.subject)[data.subject$available.site_1] == subj)))
    if (length(idx.cohort) > 10) {
      #Calculate allele frequencies against cohort-specific background
      bg.allele.subject.site_1.cohort[, subj] <- rowSums(inc.allele.subject.site_1[, idx.cohort]) / rowSums(obs.allele.subject.site_1[, idx.cohort])
    }
  }
  
  #Get site_2 background
  if (subj %in% subjects.with_site_2) {
    #Get indices of subjects to take into account as part of the global background (this includes the subject themselves, but not their families!)
    if (calculate.full) {
      idx.full <- sort(c(which(data.subject$family[data.subject$available.site_2] != data.subject[subj, "family"]), which(rownames(data.subject)[data.subject$available.site_2] == subj)))
      #Require a minimum background size
      if (length(idx.full) > 10) {
        #Calculate allele frequencies against this global background
        bg.allele.subject.site_2.full[, subj] <- rowSums(inc.allele.subject.site_2[, idx.full]) / rowSums(obs.allele.subject.site_2[, idx.full])
      }
    }
    
    #Get indices of subjects to take into account as part of the cohort-specific background (this includes the subject themselves, but not their families!)
    idx.cohort <- sort(c(which(data.subject$cohort.bg[data.subject$available.site_2] == data.subject[subj, "cohort.bg"] & data.subject$family[data.subject$available.site_2] != data.subject[subj, "family"]), which(rownames(data.subject)[data.subject$available.site_2] == subj)))
    if (length(idx.cohort) > 10) {
      #Calculate allele frequencies against cohort-specific background
      bg.allele.subject.site_2.cohort[, subj] <- rowSums(inc.allele.subject.site_2[, idx.cohort]) / rowSums(obs.allele.subject.site_2[, idx.cohort])
    }
  }
  writeLines(paste(date(), "=> Done with subject-specific background calculations for", tax.id, taxon, site_1, site_2, subj))
}
#Report
writeLines(paste(date(), "=> Done with subject-specific background calculations for", tax.id, taxon, site_1, site_2))
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate background probabilities of observing alleles, per subject & per background (full vs cohort-specific)
#=> for each pair,
#=> sum log probabilities of observation for shared alleles (S_shared ≤ 0)
#=> sum -log probabilities of observation for non-shared alleles (S_nonshared ≥ 0)
#=> sum log probabilities for least probable shared obs per pos (S_max ≤ 0)
#=> calculate score by scaling S = (S_shared - S_obs) / S_max
###################################
#=> p(site_1 & site_2)
bg.log_p.cohort.1_1 <- log10(bg.allele.subject.site_1.cohort[, curr.pairs.subject] * bg.allele.subject.site_2.cohort[, curr.pairs.subject])
#=> p(site_1 & !site_2)
bg.log_p.cohort.1_0 <- log10(bg.allele.subject.site_1.cohort[, curr.pairs.subject] * (1 - bg.allele.subject.site_2.cohort[, curr.pairs.subject]))
#=> p(!site_1 & site_2)
bg.log_p.cohort.0_1 <- log10((1 -bg.allele.subject.site_1.cohort[, curr.pairs.subject]) * bg.allele.subject.site_2.cohort[, curr.pairs.subject])
#=> p(!site_1 & !site_2)
bg.log_p.cohort.0_0 <- log10((1 -bg.allele.subject.site_1.cohort[, curr.pairs.subject]) * (1 - bg.allele.subject.site_2.cohort[, curr.pairs.subject]))

#Add for full background if needed
if (calculate.full) {
  bg.log_p.full.1_1 <- log10(bg.allele.subject.site_1.full[, curr.pairs.subject] * bg.allele.subject.site_2.full[, curr.pairs.subject])
  bg.log_p.full.1_0 <- log10(bg.allele.subject.site_1.full[, curr.pairs.subject] * (1 - bg.allele.subject.site_2.full[, curr.pairs.subject]))
  bg.log_p.full.0_1 <- log10((1 - bg.allele.subject.site_1.full[, curr.pairs.subject]) * bg.allele.subject.site_2.full[, curr.pairs.subject])
  bg.log_p.full.0_0 <- log10((1 - bg.allele.subject.site_1.full[, curr.pairs.subject]) * (1 - bg.allele.subject.site_2.full[, curr.pairs.subject]))
}
###################################

###################################
#Calculate maximum theoretical scores per subject
#=> the maximum attainable score vs background is the min(log(p)) for the shared observations (1_1 and 0_0) per each allele
S.max_theo.bg.cohort <- pmin(bg.log_p.cohort.1_1, bg.log_p.cohort.0_0)
#Re-name columns
colnames(bg.log_p.cohort.1_1) <- colnames(bg.log_p.cohort.1_0) <- colnames(bg.log_p.cohort.0_1) <- colnames(bg.log_p.cohort.0_0) <- curr.pairs.subject_timepoint
colnames(S.max_theo.bg.cohort) <- curr.pairs.subject_timepoint
if (calculate.full) {
  S.max_theo.bg.full <- pmin(bg.log_p.full.1_1, bg.log_p.full.0_0)
  colnames(bg.log_p.full.1_1) <- colnames(bg.log_p.full.1_0) <- colnames(bg.log_p.full.0_1) <- colnames(bg.log_p.full.0_0) <- curr.pairs.subject_timepoint
  colnames(S.max_theo.bg.full) <- curr.pairs.subject_timepoint
}
################################################################################
################################################################################

################################################################################
################################################################################
#Define current subjects & timepoints of interest
curr.s_t.oi <- curr.pairs.subject_timepoint

#Reduce data.timepoint to only subjects of interest
data.s_t.oi <- data.timepoint[curr.s_t.oi, ]
#Define fields to add
data.s_t.fields <- c(
  "family",
  #Horizontal coverages
  "hor_cov.site_1", "hor_cov.site_2",
  #Vertical coverages
  "ver_cov.site_1", "ver_cov.site_2",
  #specI rel. abundances
  "abd.site_1", "abd.site_2",
  #Total number of alleles w/ observation in individual
  "n.obs", "n.allele",
  #Observed allele positions per body site
  "n.obs.site_1", "n.allele.site_1", "n.obs.site_2", "n.allele.site_2",
  #Number of shared observations and shared allele observations
  "n.obs.shared", "n.allele.shared",
  #Raw Jaccard index of shared allele observations
  "jac.shared",
  #Raw scores of current observation against different backgrounds
  "S.bg_full.agree", "S.bg_full.disagree", "S.bg_full.max_theo",
  "S.bg_cohort.agree", "S.bg_cohort.disagree", "S.bg_cohort.max_theo",
  #Probability scores against full and cohort-specific background, across all pairs or individual-specific pairs 
  "S.bg_full", "S.bg_cohort",
  #Z scores against full and cohort-specific background
  "S.z.bg_full.all", "S.z.bg_full.indiv", "S.z.bg_cohort.all", "S.z.bg_cohort.indiv"
)
#Preallocate additional columns
data.s_t.oi <- cbind(data.s_t.oi, as.data.frame(matrix(nrow=nrow(data.s_t.oi), ncol=length(data.s_t.fields), dimnames=list(rownames(data.s_t.oi), data.s_t.fields))))
###################################

###################################
#Get habitat, subject.timepoint, subject, family & background cohort for each pair of samples
#=> "cp" for "current pairs"
cp.habitat <- cbind(as.character(data.sample[curr.sample_pairs[,1], "material"]), as.character(data.sample[curr.sample_pairs[,2], "material"]))
cp.subject.timepoint <- cbind(as.character(data.sample[curr.sample_pairs[,1], "subject_timepoint"]), as.character(data.sample[curr.sample_pairs[,2], "subject_timepoint"]))
cp.subject <- cbind(as.character(data.sample[curr.sample_pairs[,1], "subject"]), as.character(data.sample[curr.sample_pairs[,2], "subject"]))
cp.family <- cbind(as.character(data.sample[curr.sample_pairs[,1], "family"]), as.character(data.sample[curr.sample_pairs[,2], "family"]))
cp.cohort.bg <- cbind(as.character(data.sample[curr.sample_pairs[,1], "cohort"]), as.character(data.sample[curr.sample_pairs[,2], "cohort"]))
cp.cohort.bg[cp.cohort.bg == "MM"] <- "FR-CRC"

#Get a list of "eligible" pairs
#=> pairs which satisfy the "minimum shared observations" criterion
#=> this can take a long time
cp.min_obs <- t(apply(curr.sample_pairs, 1, function(pair) {sum(rowSums(obs.allele[, pair]) == 2) > PARAM$min.shared_obs}))
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate pairwise scores and theoretical maximum scores
##=> per background
#=> for all site_1-site_2 pairs w/in background
#=> using subject-specific allele frequency profiles
################################################################################
################################################################################
#Score pairwise observations between all pairs of samples, for different backgrounds
S.bg.full <- S.bg.full.indiv <- S.bg.cohort <- S.bg.cohort.indiv <- S.bg.full.Z <- S.bg.full.indiv.Z <- S.bg.cohort.Z <- S.bg.cohort.indiv.Z <- list()
data.bg.full <- data.bg.full.indiv <- data.bg.cohort <- data.bg.cohort.indiv <- data.frame()

#Define function to score paired observations (OR-ST)
get.S <- function(pair, s_t, bg="full", keep=rep.int(T, nrow(S.max_theo.bg.full))) {
  #Get indices of alleles to consider
  idx.obs <- keep
  
  #Return NA if there are too few shared observations
  if (sum(idx.obs) < PARAM$min.shared_obs) {return(as.numeric(rep.int(NA, 4)))}
  
  #Subset data to current pair
  curr.inc.allele <- inc.allele[, pair]
  
  #Get indices for direct comparisons
  idx.1_1 <- idx.obs & (rowSums(curr.inc.allele) == 2)
  idx.1_0 <- idx.obs & curr.inc.allele[, 1] & (!curr.inc.allele[, 2])
  idx.0_1 <- idx.obs & curr.inc.allele[, 2] & (!curr.inc.allele[, 1])
  idx.0_0 <- idx.obs & (!curr.inc.allele[, 1]) & (!curr.inc.allele[, 2])
  
  #Compute raw scores
  if (bg == "full") {
    S.agree <- (sum(bg.log_p.full.1_1[idx.1_1, s_t], na.rm=T) + sum(bg.log_p.full.0_0[idx.0_0, s_t], na.rm=T))
    S.disagree <- (sum(bg.log_p.full.1_0[idx.1_0, s_t], na.rm=T) + sum(bg.log_p.full.0_1[idx.0_1, s_t], na.rm=T))
    S.max_theo <- sum(S.max_theo.bg.full[idx.obs, s_t])
    #S.max_oral <- sum(log_p[[cohort]][idx.1_0 | idx.1_1, "1_1"], na.rm=T) + sum(log_p[[cohort]][idx.0_1 | idx.0_0, "0_0"], na.rm=T)
  } else {
    S.agree <- (sum(bg.log_p.cohort.1_1[idx.1_1, s_t], na.rm=T) + sum(bg.log_p.cohort.0_0[idx.0_0, s_t], na.rm=T))
    S.disagree <- (sum(bg.log_p.cohort.1_0[idx.1_0, s_t], na.rm=T) + sum(bg.log_p.cohort.0_1[idx.0_1, s_t], na.rm=T))
    S.max_theo <- sum(S.max_theo.bg.cohort[idx.obs, s_t])
  }
  
  #Return
  c(S.agree,
    S.disagree,
    S.max_theo,
    #S.max_oral
    #Get "Sebastian"-style transmission score: scale by maximally unlikely agreement
    (S.agree - S.disagree) / S.max_theo
    #Get "Luis" score: average agreement, scaled by total "surprise" of observation
    #(S.agree - S.disagree) / (S.agree + S.disagree),
    #Get "Simone" score: like Luis, but w/out disagreement term
    #S.agree / (S.agree + S.disagree),
    #Get "Georg" score: scaled by theoretically possible agreement, conditioned on oral sample
    #(S.agree - S.disagree) / S.max_oral
  )
}

#Define function to Z transform a distribution
z.transform <- function(x) {(x - mean(x, na.rm=T)) / sd(x, na.rm=T)}

#Iterate through subject.timepoints and calculate transmission scores
#=> using the subject-specific background allele frequencies, both within and across cohorts
for (s_t in curr.s_t.oi) {
  ###################################
  #Get current subject ID
  subj <- data.s_t.oi[s_t, "subject"]
  fam <- as.character(data.subject[subj, "family"])
  cohort <- data.s_t.oi[s_t, "cohort"]
  if (cohort == "MM") {cohort.bg <- "FR-CRC"} else {cohort.bg <- cohort}
  ###################################
  
  ###################################
  #Select alleles which have non-trivial observations in current cohort
  #=> observed in ≥1 OR and ≥1 ST sample w/in cohort
  keep.alleles.cohort <- is.finite(S.max_theo.bg.cohort[, s_t])
  if (calculate.full) {
    keep.alleles.full <- is.finite(S.max_theo.bg.full[, s_t])
  }
  
  #Skip current subject.timepoint if too few alleles with finite background frequencies are left
  if (sum(keep.alleles.cohort) < 100) {next()}
  
  #Skip if there are too few shared observations for an intra-individual comparison
  n.shared_obs <- sum(rowSums(obs.allele[, pairs.cohort[pairs.cohort.st[,1] == s_t & pairs.cohort.st[,2] == s_t, ]]) == 2)
  if (n.shared_obs < PARAM$min.shared_obs) {writeLines(paste(date(), "=> Too few shared intra-individual observations:", n.shared_obs, tax.id, taxon, s_t)); next()}
  ###################################
  
  ###################################
  #Get current cohort-specific pairs
  pairs.cohort <- curr.sample_pairs[cp.cohort.bg[,1] == cohort.bg & cp.cohort.bg[,2] == cohort.bg, ]
  #Skip if there are too few pairs available
  if (nrow(pairs.cohort) < 2) {next()}
  #Inherit per-cohort annotations
  pairs.cohort.st <- cbind(as.character(data.sample[pairs.cohort[,1], "subject_timepoint"]), as.character(data.sample[pairs.cohort[,2], "subject_timepoint"]))
  pairs.cohort.subject <- cbind(as.character(data.sample[pairs.cohort[,1], "subject"]), as.character(data.sample[pairs.cohort[,2], "subject"]))
  pairs.cohort.family <- cbind(as.character(data.sample[pairs.cohort[,1], "family"]), as.character(data.sample[pairs.cohort[,2], "family"]))
  pairs.cohort.min_obs <- cp.min_obs[cp.cohort.bg[,1] == cohort.bg & cp.cohort.bg[,2] == cohort.bg]
  ###################################
  
  ###################################
  #Calculate scores across all pairs (THIS TAKES TIME!)
  #=> subsample to a (lower) maximum number of comparisons 
  #=> for both *.full and *.cohort
  #
  #Obtain background score distributions
  #=> pre-require minimum shared observations
  #=> ignore comparisons to subject's family
  #=> ignore comparisons within families
  #=> ignore intra-subject comparisons (implicit in family blocking)
  #=> but consider all comparisons to current subject.timepoint
  ###################################
  #=> against a full background
  if (calculate.full) {
    idx.pairs.full.bg <- cp.min_obs & (cp.family[,1] != cp.family[,2]) & cp.family[,1] != fam & cp.family[,2] != fam
    idx.pairs.full.bg[(cp.subject.timepoint[,1] == s_t | cp.subject.timepoint[,2] == s_t) & (cp.subject.timepoint[,1] != cp.subject.timepoint[,2])] <- TRUE
    if (sum(idx.pairs.full.bg) > PARAM$max.background_size) {
      curr.S.full <- t(apply(curr.sample_pairs[sample(which(idx.pairs.full.bg), size=PARAM$max.background_size, replace=F), ], 1, get.S, s_t=s_t, bg="full", keep=keep.alleles.full))
    } else {
      curr.S.full <- t(apply(curr.sample_pairs[idx.pairs.full.bg, ], 1, get.S, s_t=s_t, bg="full", keep=keep.alleles.full))
    }
    colnames(curr.S.full) <- c("S.agree", "S.disagree", "S.max_theo", "S")
    S.bg.full[[s_t]] <- na.omit(curr.S.full[, "S"])
  }
  ###################################
  #=> against the cohort-specific background
  idx.pairs.cohort.bg <- pairs.cohort.min_obs & (pairs.cohort.family[,1] != pairs.cohort.family[,2]) & pairs.cohort.family[,1] != fam & pairs.cohort.family[,2] != fam
  idx.pairs.cohort.bg[(pairs.cohort.st[,1] == s_t | pairs.cohort.st[,2] == s_t) & (pairs.cohort.st[,1] != pairs.cohort.st[,2])] <- TRUE
  if (sum(idx.pairs.cohort.bg) > PARAM$max.background_size) {
    curr.S.cohort <- t(apply(pairs.cohort[sample(which(idx.pairs.cohort.bg), size=PARAM$max.background_size, replace=F), ], 1, get.S, s_t=s_t, bg=cohort.bg, keep=keep.alleles.cohort))
  } else {
    curr.S.cohort <- t(apply(pairs.cohort[idx.pairs.cohort.bg, ], 1, get.S, s_t=s_t, bg=cohort.bg, keep=keep.alleles.cohort))
  }
  colnames(curr.S.cohort) <- c("S.agree", "S.disagree", "S.max_theo", "S")
  S.bg.cohort[[s_t]] <- na.omit(curr.S.cohort[, "S"])
  ###################################
  #Obtain individual-specific background score distributions
  #=> only comparisons to current subject.timepoint's oral or gut sample
  #=> include comparisons within current subjects family
  ###################################
  #=> with allele frequencies calculated against the full background
  if (calculate.full) {
    idx.pairs.full.bg.indiv <- (cp.subject.timepoint[,1] == s_t | cp.subject.timepoint[,2] == s_t) & (cp.subject[,1] != cp.subject[,2])
    curr.S.full.indiv <- t(apply(curr.sample_pairs[idx.pairs.full.bg.indiv, ], 1, get.S, s_t=s_t, bg="full", keep=keep.alleles.full))
    colnames(curr.S.full.indiv) <- c("S.agree", "S.disagree", "S.max_theo", "S")
    S.bg.full.indiv[[s_t]] <- na.omit(curr.S.full.indiv[, "S"])
  }
  #=> with allele frequencies calculated against the cohort-specific background
  idx.pairs.cohort.bg.indiv <- (pairs.cohort.st[,1] == s_t | pairs.cohort.st[,2] == s_t) & (pairs.cohort.subject[,1] != pairs.cohort.subject[,2])
  curr.S.cohort.indiv <- t(apply(pairs.cohort[idx.pairs.cohort.bg.indiv, ], 1, get.S, s_t=s_t, bg=cohort.bg, keep=keep.alleles.cohort))
  colnames(curr.S.cohort.indiv) <- c("S.agree", "S.disagree", "S.max_theo", "S")
  S.bg.cohort.indiv[[s_t]] <- na.omit(curr.S.cohort.indiv[, "S"])
  ###################################
  #Get current intra-individual (intra-timepoint) comparison
  ###################################
  collect.S.intra.cohort <- get.S(pairs.cohort[pairs.cohort.st[,1] == s_t & pairs.cohort.st[,2] == s_t, ], s_t=s_t, bg=cohort.bg, keep=keep.alleles.cohort)
  names(collect.S.intra.cohort) <- c("S.agree", "S.disagree", "S.max_theo", "S")
  if (calculate.full) {
    collect.S.intra.full <- get.S(curr.sample_pairs[cp.subject.timepoint[,1] == s_t & cp.subject.timepoint[,2] == s_t, ], s_t=s_t, bg="full", keep=keep.alleles.full)
    names(collect.S.intra.full) <- c("S.agree", "S.disagree", "S.max_theo", "S")
  }
  ###################################
  
  ###################################
  #Re-scale (Z-transform) background distributions
  ###################################
  if (calculate.full) {
    S.bg.full.Z[[s_t]] <- z.transform(S.bg.full[[s_t]])
    S.bg.full.indiv.Z[[s_t]] <- z.transform(S.bg.full.indiv[[s_t]])
  }
  S.bg.cohort.Z[[s_t]] <- z.transform(S.bg.cohort[[s_t]])
  S.bg.cohort.indiv.Z[[s_t]] <- z.transform(S.bg.cohort.indiv[[s_t]])
  ###################################
  #Z-transform current intra-individual comparisons
  ###################################
  if (calculate.full) {
    S.z.intra.full <- (collect.S.intra.full["S"] - mean(S.bg.full[[s_t]], na.rm=T)) / sd(S.bg.full[[s_t]], na.rm=T)
    S.z.intra.full.indiv <- (collect.S.intra.full["S"] - mean(S.bg.full.indiv[[s_t]], na.rm=T)) / sd(S.bg.full.indiv[[s_t]], na.rm=T)
  }
  S.z.intra.cohort <- (collect.S.intra.cohort["S"] - mean(S.bg.cohort[[s_t]], na.rm=T)) / sd(S.bg.cohort[[s_t]], na.rm=T)
  S.z.intra.cohort.indiv <- (collect.S.intra.cohort["S"] - mean(S.bg.cohort.indiv[[s_t]], na.rm=T)) / sd(S.bg.cohort.indiv[[s_t]], na.rm=T)
  ###################################
  
  ###################################
  #Store data on background distributions
  ###################################
  #Full background
  if (calculate.full) {
    data.bg.full <- rbind(data.bg.full, data.frame(
      s_t = s_t, n.comparisons = sum(idx.pairs.full.bg), S.mean = mean(S.bg.full[[s_t]], na.rm=T), S.sd = sd(S.bg.full[[s_t]], na.rm=T)
    ))
    #Full background, per current s_t
    data.bg.full.indiv <- rbind(data.bg.full.indiv, data.frame(
      s_t = s_t, n.comparisons = sum(idx.pairs.full.bg.indiv), S.mean = mean(S.bg.full.indiv[[s_t]], na.rm=T), S.sd = sd(S.bg.full.indiv[[s_t]], na.rm=T)
    ))
  }
  #Cohort-specific background
  data.bg.cohort <- rbind(data.bg.cohort, data.frame(
    s_t = s_t, n.comparisons = sum(idx.pairs.cohort.bg), S.mean = mean(S.bg.cohort[[s_t]], na.rm=T), S.sd = sd(S.bg.cohort[[s_t]], na.rm=T)
  ))
  #Cohort-specific background, per current s_t
  data.bg.cohort.indiv <- rbind(data.bg.cohort.indiv, data.frame(
    s_t = s_t, n.comparisons = sum(idx.pairs.cohort.bg.indiv), S.mean = mean(S.bg.cohort.indiv[[s_t]], na.rm=T), S.sd = sd(S.bg.cohort.indiv[[s_t]], na.rm=T)
  ))
  ###################################
  
  ###################################
  #Store data on intra-individual transmission
  #=> metadata (coverages, etc.)
  ###################################
  #Get current oral and gut sample id
  if (site_1 == "saliva") {intra.site_1 <- data.timepoint[s_t, "sample.oral"]}
  if (site_1 == "dental_plaque") {intra.site_1 <- data.timepoint[s_t, "sample.plaque"]}
  if (site_1 == "stool") {intra.site_1 <- data.timepoint[s_t, "sample.site_2"]}
  if (site_2 == "saliva") {intra.site_2 <- data.timepoint[s_t, "sample.oral"]}
  if (site_2 == "dental_plaque") {intra.site_2 <- data.timepoint[s_t, "sample.plaque"]}
  if (site_2 == "stool") {intra.site_2 <- data.timepoint[s_t, "sample.gut"]}
  ###################################
  #Update data.s_t.oi
  ###################################
  data.s_t.oi[s_t, "family"] <- fam
  #Coverage data
  data.s_t.oi[s_t, "hor_cov.site_1"] <- data.hor_cov[tax.id, intra.site_1]
  data.s_t.oi[s_t, "hor_cov.site_2"] <- data.hor_cov[tax.id, intra.site_2]
  data.s_t.oi[s_t, "ver_cov.site_1"] <- data.ver_cov[tax.id, intra.site_1]
  data.s_t.oi[s_t, "ver_cov.site_2"] <- data.ver_cov[tax.id, intra.site_2]
  #Abundance data
  data.s_t.oi[s_t, "abd.site_1"] <- data.specI.rel[tax.id, intra.site_1]
  data.s_t.oi[s_t, "abd.site_2"] <- data.specI.rel[tax.id, intra.site_2]
  #Allele observation and count statistics
  data.s_t.oi[s_t, "n.obs"] <- sum(rowSums(obs.allele[, c(intra.site_1, intra.site_2)]) >= 1)
  data.s_t.oi[s_t, "n.allele"] <- sum(rowSums(inc.allele[, c(intra.site_1, intra.site_2)]) >= 1)
  data.s_t.oi[s_t, "n.obs.site_1"] <- sum(obs.allele[, intra.site_1])
  data.s_t.oi[s_t, "n.allele.site_1"] <- sum(inc.allele[, intra.site_1])
  data.s_t.oi[s_t, "n.obs.site_2"] <- sum(obs.allele[, intra.site_2])
  data.s_t.oi[s_t, "n.allele.site_2"] <- sum(inc.allele[, intra.site_2])
  data.s_t.oi[s_t, "n.obs.shared"] <- sum(rowSums(obs.allele[, c(intra.site_1, intra.site_2)]) == 2)
  data.s_t.oi[s_t, "n.allele.shared"] <- sum(rowSums(inc.allele[, c(intra.site_1, intra.site_2)]) == 2)
  #Jaccard overlap of observed alleles
  data.s_t.oi[s_t, "jac.shared"] <- data.s_t.oi[s_t, "n.allele.shared"] / data.s_t.oi[s_t, "n.obs.shared"]
  #Scores on full background
  if (calculate.full) {
    data.s_t.oi[s_t, "S.bg_full.agree"] <- collect.S.intra.full["S.agree"]
    data.s_t.oi[s_t, "S.bg_full.disagree"] <- collect.S.intra.full["S.disagree"]
    data.s_t.oi[s_t, "S.bg_full.max_theo"] <- collect.S.intra.full["S.max_theo"]
    data.s_t.oi[s_t, "S.bg_full"] <- collect.S.intra.full["S"]
  } else {
    data.s_t.oi[s_t, "S.bg_full.agree"] <- data.s_t.oi[s_t, "S.bg_full.disagree"] <- data.s_t.oi[s_t, "S.bg_full.max_theo"] <-data.s_t.oi[s_t, "S.bg_full"] <- NA
  }
  #Scores on cohort-specific background
  data.s_t.oi[s_t, "S.bg_cohort.agree"] <- collect.S.intra.cohort["S.agree"]
  data.s_t.oi[s_t, "S.bg_cohort.disagree"] <- collect.S.intra.cohort["S.disagree"]
  data.s_t.oi[s_t, "S.bg_cohort.max_theo"] <- collect.S.intra.cohort["S.max_theo"]
  data.s_t.oi[s_t, "S.bg_cohort"] <- collect.S.intra.cohort["S"]
  #Z-transformed scores
  if (calculate.full) {
    data.s_t.oi[s_t, "S.z.bg_full.all"] <- S.z.intra.full
    data.s_t.oi[s_t, "S.z.bg_full.indiv"] <- S.z.intra.full.indiv
  } else {
    data.s_t.oi[s_t, "S.z.bg_full.all"] <- data.s_t.oi[s_t, "S.z.bg_full.indiv"] <- NA
  }
  data.s_t.oi[s_t, "S.z.bg_cohort.all"] <- S.z.intra.cohort
  data.s_t.oi[s_t, "S.z.bg_cohort.indiv"] <- S.z.intra.cohort.indiv
  ###################################
  
  ###################################
  writeLines(paste(date(), "=> Done with subject.timepoint", s_t, paste0("(", which(curr.s_t.oi == s_t), " out of ", length(curr.s_t.oi), ")"), "in", tax.id, taxon))
  ###################################
}

###################################
#Save current data
save(data.bg.full, data.bg.full.indiv, data.bg.cohort, data.bg.cohort.indiv, data.s_t.oi, file=paste0(PARAM$folder.results, "detect.transmission/", tax.id, ".", print.taxon,  ".detect_transmission.", site_1.short, "_", site_2.short, ".RData"))
save(S.bg.full, S.bg.full.indiv, S.bg.cohort, S.bg.cohort.indiv, file=paste0(PARAM$folder.results, "detect.transmission/", tax.id, ".", print.taxon,  ".background_scores.", site_1.short, "_", site_2.short, ".RData"))
#save(S.bg.full.Z, S.bg.full.indiv.Z, S.bg.cohort.Z, S.bg.cohort.indiv.Z, file=paste0(PARAM$folder.results, "detect.transmission/", tax.id, ".", print.taxon,  ".background_scores.z_transformed.", site_1.short, "_", site_2.short, ".RData"))
###################################

#Report
writeLines(paste(date(), "=> Done with", tax.id, taxon))
################################################################################
################################################################################



q(save="no")
