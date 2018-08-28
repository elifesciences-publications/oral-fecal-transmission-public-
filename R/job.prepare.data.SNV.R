#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Prepare SNV data for further analyses
#
#=> load sample metadata
#=> load SNV tables per taxon
#=> filter list of taxa-of-interest
#=> filter SNP tables for those taxa and samples
#=> save and export individual datafiles in R format
#
#2017-04-26
#sebastian.schmidt@embl.de
################################################################################


################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("gtools", warn.conflicts=F, quietly=T)
################################################################################
################################################################################


################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.base <- paste0("###################")
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNP_data <- paste0(PARAM$folder.data, "metaSNP_output/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")

#Read current taxon of interest from command line
args <- commandArgs(TRUE)
tax.id <- args[2]
################################################################################
################################################################################


################################################################################
################################################################################
#Load sample metadata
load(paste0(PARAM$folder.parameters, "data.sample.RData"))
load(paste0(PARAM$folder.parameters, "data.timepoint.RData"))
samples <- rownames(data.sample)
n.samples <- length(samples)

#Load specI abundance data
load(paste0(PARAM$folder.data, "data.specI.RData"))

#Load taxonomy data
load(paste0(PARAM$folder.data, "data.taxonomy.RData"))
################################################################################
################################################################################


################################################################################
################################################################################
#Load and process SNP count tables per taxon
#=> process population SNPs and individual SNPs separately
#=> skip all taxa which have not been observed in at least one oral and gut sample
################################################################################
################################################################################
#"Unfix" current taxonomy ID
#=> SNV tables are saved under a non-updated tax ID naming scheme, so to load tables, taxonomies have to made compatible
#=> with the exception of a handful of (mismatching) tax IDs, this step will not do anything
tax.id.unfixed <- unfix.tax_ids(tax.id)

#Get current taxon name 
curr.taxon <- data.taxonomy$Scientific_Name[data.taxonomy$Tax_ID == tax.id]
print.taxon <- gsub(" ", "_", curr.taxon)
print.taxon <- gsub("/", "_", print.taxon)

#Get current file names for current taxon
curr.file.pop <- paste0(PARAM$folder.SNP_data, "filtered.pop/", tax.id.unfixed, ".filtered.freq.gz")
curr.file.indiv <- paste0(PARAM$folder.SNP_data, "filtered.indiv/", tax.id.unfixed, ".filtered.freq.gz")

#Report
writeLines(paste(date(), "=> Processing", tax.id, tax.id.unfixed, curr.taxon))
###################################

###################################
#Read data
###################################
#Check if population-level SNV file exists
if (! file.exists(curr.file.pop)) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> No population-level SNV file available.\n")); q(save="no")}

#Read data
#=> first, read only the first few lines to check if data comes from ≥1 habitat
tmp.dat <- read.table(curr.file.pop, sep="\t", header=T, nrows=2, row.names=1)

#Process sample names
tmp.samples <- colnames(tmp.dat)
curr.samples <- curr.samples <- gsub("[.]", "_", tmp.samples)
curr.samples <- gsub("_sorted_unique_bam", "", curr.samples)

#Check if current taxon was observed in more than one habitat
curr.habitats <- as.character(data.sample[curr.samples, "body_site"])
n.hab <- length(unique(curr.habitats))

#Skip current taxon if only samples from one habitat are covered
if (n.hab == 1) {
  writeLines(paste0("Skip: only one habitat covered (", unique(curr.habitats), ")\n"))
  q(save="no")
}

#Report
writeLines(paste0(date(), " => ", tax.id, " ", curr.taxon, " => Both habitats covered: oral (", length(which(curr.habitats == "oral")), "), gut (", length(which(curr.habitats == "gut")), ")."))
###################################

###################################
#Read data
#=> for real, this time
tmp.dat <- read.table(curr.file.pop, sep="\t", header=T, row.names=1)

#Report
writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", nrow(tmp.dat), "raw pop-level SNP positions."))
###################################

###################################
#Check whether individual SNPs exist for current taxon
if (file.exists(curr.file.indiv)) {
  #Read individual SNP data
  tmp.dat.indiv <- read.table(curr.file.indiv, sep="\t", header=T, row.names=1)
  tmp.mat.pos.indiv <- as.matrix(tmp.dat.indiv)
  #Report
  writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", nrow(tmp.dat.indiv), "raw indiv-level SNP positions."))
  
  #Concatenate
  tmp.tmp.dat <- rbind(tmp.dat, tmp.dat.indiv)
  #Re-order
  new.order <- mixedorder(rownames(tmp.tmp.dat))
  tmp.dat <- tmp.tmp.dat[new.order, ]
  rm(tmp.tmp.dat); invisible(gc())
}

#Skip if too few positions were encountered
if (nrow(tmp.dat) < 10) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Skip: too few positions.\n")); q(save="no")}
###################################

###################################
#Process sample names
tmp.samples <- colnames(tmp.dat)
curr.samples <- curr.samples <- gsub("[.]", "_", tmp.samples)
curr.samples <- gsub("_sorted_unique_bam", "", curr.samples)
colnames(tmp.dat) <- curr.samples

#"Vertically" filter out samples which have no information (-1) across all positions
#=> moreover, remove samples that have previously been dropped based on abundance/coverage data
#=> this is also a sanity check
tmp.mat.pos <- as.matrix(tmp.dat)
keep.samples <- (colSums(tmp.mat.pos) > (-nrow(tmp.mat.pos))) & (curr.samples %in% rownames(data.sample))
curr.samples <- curr.samples[keep.samples]
tmp.mat.pos <- tmp.mat.pos[, keep.samples]

#Report
writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", nrow(tmp.dat), "SNP positions across", length(curr.samples), "samples."))

if (length(curr.samples) < 2) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Skip: too few samples covered appropriately.\n")); q(save="no")}

#Report if a sample ID does not match the original table
if (length(which(! curr.samples %in% rownames(data.sample))) > 0) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Unknown samples:", curr.samples[! curr.samples %in% rownames(data.sample)]))}

#Get current habitats (from fixed sample names)
curr.habitats <- as.character(data.sample[curr.samples, "body_site"])
n.hab <- length(unique(curr.habitats))
#Get current sampling materials
curr.materials <- as.character(data.sample[curr.samples, "material"])
n.materials <- length(unique(curr.materials))
###################################

###################################
#Extract current information per position
tmp.pos_info <- do.call(rbind, strsplit(rownames(tmp.dat), ":"))
tmp.ref_base <- substr(tmp.pos_info[,4], start=1, stop=1)
tmp.base <- substr(tmp.pos_info[,4], start=3, stop=3)
tmp.subst_effect <- substr(tmp.pos_info[,5], start=1, stop=1)
data.pos <- data.frame(Contig=tmp.pos_info[,1], Gene=tmp.pos_info[,2], Position=as.numeric(tmp.pos_info[,3]), Substitution=tmp.pos_info[,4], Ref_Base=tmp.ref_base, Base=tmp.base, Effect=tmp.subst_effect, Effect_Full=tmp.pos_info[,5])
rownames(data.pos) <- paste0(data.pos$Contig, ":", data.pos$Position, ":", data.pos$Base)
rm(tmp.dat); invisible(gc())
###################################

###################################
#Get clean positional frequency matrix
mat.pos <- tmp.mat.pos[, curr.samples]
rownames(mat.pos) <- rownames(data.pos)
rm(tmp.mat.pos); invisible(gc())

#Add a minute pseudo-count to ref positions and transform "-1"s to "0s"
mat.pos[mat.pos == 0] <- 10^-12
mat.pos[mat.pos == -1] <- 0

#Make matrix sparse
#=> but only if that won't create overflow trouble for the Matrix package
if (length(mat.pos) < (2^31 / 2)) {mat.pos <- Matrix(mat.pos, sparse=T)}

#"Horizontally" filter positions
#=> to select only those with positive non-reference frequencies
#=> this is a sanity check, as well
tmp.freq.sum <- rowSums(mat.pos > 0)
keep.pos <- tmp.freq.sum > 0
if (length(which(keep.pos)) <= 10) {writeLines(paste(date(), "=>", tax.id, curr.taxon, "=> Skip: too few positions after filtering for non-reference.")); q(save="no")}

#Report
writeLines(paste(date(), "=>", tax.id, curr.taxon, "=>", length(which(keep.pos)), "filtered non-reference SNP positions."))

#Filter
mat.pos <- mat.pos[keep.pos, ]
data.pos <- data.pos[keep.pos, ]
###################################

###################################
#Turn mat.pos into per-allele table
mat.tmp <- 1 - as.matrix(mat.pos)
#Account for observations of other alleles only ("1" in the original matrix => everything non-reference)
mat.tmp[mat.tmp == 0] <- 10^-12
#Account for missing observations (zeroes in the original allele freq matrix)
mat.tmp[mat.tmp == 1] <- 0
#Account for positions that were all ref (10^-12 in the original allele freq matrix)
mat.tmp[mat.tmp == 1 - 10^-12] <- 1
#Re-set rownames
rownames(mat.tmp) <- paste(data.pos$Contig, data.pos$Position, data.pos$Ref_Base, sep=":")
#Merge matrices and re-sparsify
if ((length(mat.pos) + length(mat.tmp)) < (2^31 / 2)) {mat.allele <- Matrix(rbind(mat.pos, mat.tmp), sparse = T)} else {mat.allele <- rbind(mat.pos, mat.tmp)}
keep.alleles <- mixedsort(unique(rownames(mat.allele)))     #use gtools::mixedsort to sort naturally
mat.allele <- mat.allele[keep.alleles, ]
rm(mat.tmp); invisible(gc())
#Get per-allele data
tmp.allele_info <- do.call(rbind, strsplit(rownames(mat.allele), ":"))
data.allele <- data.frame(Contig=tmp.allele_info[,1], Position=as.numeric(tmp.allele_info[,2]), Base=tmp.allele_info[,3])
rownames(data.allele) <- rownames(mat.allele)
###################################

###################################
#Filter positions
#=> keep only positions which have been called in at least one oral and at least one stool sample
###################################
#Get summed frequencies per position
if (length(which(curr.habitats == "oral")) == 1) {
  sums.oral.pos <- mat.pos[, curr.habitats == "oral"]
  sums.oral.allele <- mat.allele[, curr.habitats == "oral"]
} else {
  sums.oral.pos <- rowSums(mat.pos[, curr.habitats == "oral"], na.rm=T)
  sums.oral.allele <- rowSums(mat.allele[, curr.habitats == "oral"], na.rm=T)
}
if (length(which(curr.habitats == "gut")) == 1) {
  sums.gut.pos <- mat.pos[, curr.habitats == "gut"]
  sums.gut.allele <- mat.allele[, curr.habitats == "gut"]
} else {
  sums.gut.pos <- rowSums(mat.pos[, curr.habitats == "gut"], na.rm=T)
  sums.gut.allele <- rowSums(mat.allele[, curr.habitats == "gut"], na.rm=T)
}
#Get posision and allele sums per sampling material
if (length(which(curr.materials == "saliva")) == 1) {
  sums.saliva.pos <- mat.pos[, curr.materials == "saliva"]
  sums.saliva.allele <- mat.allele[, curr.materials == "saliva"]
} else {
  sums.saliva.pos <- rowSums(mat.pos[, curr.materials == "saliva"], na.rm=T)
  sums.saliva.allele <- rowSums(mat.allele[, curr.materials == "saliva"], na.rm=T)
}
if (length(which(curr.materials == "dental_plaque")) == 1) {
  sums.plaque.pos <- mat.pos[, curr.materials == "dental_plaque"]
  sums.plaque.allele <- mat.allele[, curr.materials == "dental_plaque"]
} else {
  sums.plaque.pos <- rowSums(mat.pos[, curr.materials == "dental_plaque"], na.rm=T)
  sums.plaque.allele <- rowSums(mat.allele[, curr.materials == "dental_plaque"], na.rm=T)
}

#Annotate positions and alleles covered in...
#=> at least one oral sample
data.pos$covered.oral <- sums.oral.pos > 0
data.allele$covered.oral <- sums.oral.allele > 0
#=> at least one saliva sample
data.pos$covered.saliva <- sums.saliva.pos > 0
data.allele$covered.saliva <- sums.saliva.allele > 0
#=> at least one plaque sample
data.pos$covered.plaque <- sums.plaque.pos > 0
data.allele$covered.plaque <- sums.plaque.allele > 0
#=> at least one stool sample
data.pos$covered.gut <- sums.gut.pos > 0
data.allele$covered.gut <- sums.gut.allele > 0

#Annotate positions that are an SNV in...
#=> at least one oral sample
data.pos$snv.oral <- sums.oral.pos > length(which(curr.habitats == "oral")) * 10^-12
#=> at least one saliva sample
data.pos$snv.saliva <- sums.saliva.pos > length(which(curr.materials == "saliva")) * 10^-12
#=> at least one plaque sample
data.pos$snv.plaque <- sums.plaque.pos > length(which(curr.materials == "dental_plaque")) * 10^-12
#=> at least one stool sample
data.pos$snv.gut <- sums.gut.pos > length(which(curr.habitats == "gut")) * 10^-12
###################################

###################################
#Establish intra-timepoint observed positions, SNVs and alleles
tp.pos.oral_gut <- tp.pos.saliva_gut <- tp.pos.plaque_gut <- tp.pos.saliva_plaque <- numeric()
tp.snv.oral_gut <- tp.snv.saliva_gut <- tp.snv.plaque_gut <- tp.snv.saliva_plaque <- numeric()
tp.allele.oral_gut <- tp.allele.saliva_gut <- tp.allele.plaque_gut <- tp.allele.saliva_plaque <- numeric()
n.st.oral_gut <- n.st.saliva_gut <- n.st.plaque_gut <- n.st_saliva_plaque <- 0
#Iterate throught subject.timepoints
for (st in rownames(data.timepoint)) {
  #Get current samples
  curr.sample.saliva <- data.timepoint[st, "sample.oral"]
  curr.sample.plaque <- data.timepoint[st, "sample.plaque"]
  curr.sample.gut <- data.timepoint[st, "sample.gut"]
  
  #Get data on saliva
  if ((!is.na(curr.sample.saliva)) & curr.sample.saliva %in% curr.samples) {
    curr.tp.pos.saliva <- which(mat.pos[, curr.sample.saliva] > 0)
    curr.tp.snv.saliva <- which(mat.pos[, curr.sample.saliva] > 10^-12)
    curr.tp.allele.saliva <- which(mat.allele[, curr.sample.saliva] > 0)
  } else {
    curr.tp.pos.saliva <- curr.tp.snv.saliva <- curr.tp.allele.saliva <- numeric()
  }
  
  #Get data on plaque
  if ((!is.na(curr.sample.plaque)) & curr.sample.plaque %in% curr.samples) {
    curr.tp.pos.plaque <- which(mat.pos[, curr.sample.plaque] > 0)
    curr.tp.snv.plaque <- which(mat.pos[, curr.sample.plaque] > 10^-12)
    curr.tp.allele.plaque <- which(mat.allele[, curr.sample.plaque] > 0)
  } else {
    curr.tp.pos.plaque <- curr.tp.snv.plaque <- curr.tp.allele.plaque <- numeric()
  }
  
  #Get data on gut
  if ((!is.na(curr.sample.gut)) & curr.sample.gut %in% curr.samples) {
    curr.tp.pos.gut <- which(mat.pos[, curr.sample.gut] > 0)
    curr.tp.snv.gut <- which(mat.pos[, curr.sample.gut] > 10^-12)
    curr.tp.allele.gut <- which(mat.allele[, curr.sample.gut] > 0)
  } else {
    curr.tp.pos.gut <- curr.tp.snv.gut <- curr.tp.allele.gut <- numeric()
  }
  
  #Process oral <-> gut
  if ((curr.sample.saliva %in% curr.samples | curr.sample.plaque %in% curr.samples) & curr.sample.gut %in% curr.samples) {
    tp.pos.oral_gut <- unique(c(tp.pos.oral_gut, intersect(union(curr.tp.pos.saliva, curr.tp.pos.plaque), curr.tp.pos.gut)))
    tp.snv.oral_gut <- unique(c(tp.snv.oral_gut, intersect(union(curr.tp.snv.saliva, curr.tp.snv.plaque), curr.tp.snv.gut)))
    tp.allele.oral_gut <- unique(c(tp.allele.oral_gut, intersect(union(curr.tp.allele.saliva, curr.tp.allele.plaque), curr.tp.allele.gut)))
    n.st.oral_gut <- n.st.oral_gut + 1
  }
  
  #Process saliva <-> gut
  if (curr.sample.saliva %in% curr.samples & curr.sample.gut %in% curr.samples) {
    tp.pos.saliva_gut <- unique(c(tp.pos.saliva_gut, intersect(curr.tp.pos.saliva, curr.tp.pos.gut)))
    tp.snv.saliva_gut <- unique(c(tp.snv.saliva_gut, intersect(curr.tp.snv.saliva, curr.tp.snv.gut)))
    tp.allele.saliva_gut <- unique(c(tp.allele.saliva_gut, intersect(curr.tp.allele.saliva, curr.tp.allele.gut)))
    n.st.saliva_gut <- n.st.saliva_gut + 1
  }
  
  #Process plaque <-> gut
  if (curr.sample.plaque %in% curr.samples & curr.sample.gut %in% curr.samples) {
    tp.pos.plaque_gut <- unique(c(tp.pos.plaque_gut, intersect(curr.tp.pos.plaque, curr.tp.pos.gut)))
    tp.snv.plaque_gut <- unique(c(tp.snv.plaque_gut, intersect(curr.tp.snv.plaque, curr.tp.snv.gut)))
    tp.allele.plaque_gut <- unique(c(tp.allele.plaque_gut, intersect(curr.tp.allele.plaque, curr.tp.allele.gut)))
    n.st.plaque_gut <- n.st.plaque_gut + 1
  }
  
  #Process saliva <-> plaque
  if (curr.sample.saliva %in% curr.samples & curr.sample.plaque %in% curr.samples) {
    tp.pos.saliva_plaque <- unique(c(tp.pos.saliva_plaque, intersect(curr.tp.pos.saliva, curr.tp.pos.plaque)))
    tp.snv.saliva_plaque <- unique(c(tp.snv.saliva_plaque, intersect(curr.tp.snv.saliva, curr.tp.snv.plaque)))
    tp.allele.saliva_plaque <- unique(c(tp.allele.saliva_plaque, intersect(curr.tp.allele.saliva, curr.tp.allele.plaque)))
    n.st_saliva_plaque <- n.st_saliva_plaque + 1
  }
}

###################################

###################################
#Store data on current taxon
curr.data.taxa <- data.frame(
  SNV_data.available = TRUE,
  SNV_data.oral_gut = TRUE,
  n.pos.total = nrow(data.pos),
  n.allele.total = nrow(data.allele),
  #Oral data
  n.pos.oral = length(which(data.pos$covered.oral)),
  n.snv.oral = length(which(data.pos$snv.oral)),
  n.allele.oral = length(which(data.allele$covered.oral)),
  #Saliva data
  n.pos.saliva = length(which(data.pos$covered.saliva)),
  n.snv.saliva = length(which(data.pos$snv.saliva)),
  n.allele.saliva = length(which(data.allele$covered.saliva)),
  #Plaque data
  n.pos.plaque = length(which(data.pos$covered.plaque)),
  n.snv.plaque = length(which(data.pos$snv.plaque)),
  n.allele.plaque = length(which(data.allele$covered.plaque)),
  #Gut data
  n.pos.gut = length(which(data.pos$covered.gut)),
  n.snv.gut = length(which(data.pos$snv.gut)),
  n.allele.gut = length(which(data.allele$covered.gut)),
  #Oral <-> gut data
  n.pos.shared.oral_gut = length(which(data.pos$covered.oral & data.pos$covered.gut)),
  n.snv.shared.oral_gut= length(which(data.pos$snv.oral & data.pos$snv.gut)),
  n.allele.shared.oral_gut = length(which(data.allele$covered.oral & data.allele$covered.gut)),
  #Saliva <-> gut data
  n.pos.shared.saliva_gut = length(which(data.pos$covered.saliva & data.pos$covered.gut)),
  n.snv.shared.saliva_gut= length(which(data.pos$snv.saliva & data.pos$snv.gut)),
  n.allele.shared.saliva_gut = length(which(data.allele$covered.saliva & data.allele$covered.gut)),
  #Plaque <-> gut data
  n.pos.shared.plaque_gut = length(which(data.pos$covered.plaque & data.pos$covered.gut)),
  n.snv.shared.plaque_gut= length(which(data.pos$snv.plaque & data.pos$snv.gut)),
  n.allele.shared.plaque_gut = length(which(data.allele$covered.plaque & data.allele$covered.gut)),
  #Saliva <-> plaque data
  n.pos.shared.saliva_plaque = length(which(data.pos$covered.plaque & data.pos$covered.saliva)),
  n.snv.shared.saliva_plaque= length(which(data.pos$snv.plaque & data.pos$snv.saliva)),
  n.allele.shared.saliva_plaque = length(which(data.allele$covered.plaque & data.allele$covered.saliva)),
  #Stats on intra-timepoint observed positions, SNVs and alleles
  #=> oral <-> gut
  n.st.oral_gut = n.st.oral_gut,
  n.pos.intra_tp.oral_gut = length(tp.pos.oral_gut),
  n.snv.intra_tp.oral_gut = length(tp.snv.oral_gut),
  n.allele.intra_tp.oral_gut = length(tp.allele.oral_gut),
  #=> saliva <-> gut
  n.st.saliva_gut = n.st.saliva_gut,
  n.pos.intra_tp.saliva_gut = length(tp.pos.saliva_gut),
  n.snv.intra_tp.saliva_gut = length(tp.snv.saliva_gut),
  n.allele.intra_tp.saliva_gut = length(tp.allele.saliva_gut),
  #=> plaque <-> gut
  n.st.plaque_gut = n.st.plaque_gut,
  n.pos.intra_tp.plaque_gut = length(tp.pos.plaque_gut),
  n.snv.intra_tp.plaque_gut = length(tp.snv.plaque_gut),
  n.allele.intra_tp.plaque_gut = length(tp.allele.plaque_gut),
  #=> saliva <-> plaque
  n.st_saliva_plaque = n.st_saliva_plaque,
  n.pos.intra_tp.saliva_plaque = length(tp.pos.saliva_plaque),
  n.snv.intra_tp.saliva_plaque = length(tp.snv.saliva_plaque),
  n.allele.intra_tp.saliva_plaque = length(tp.allele.saliva_plaque)
)
###################################

###################################
#Store
save(data.allele, mat.allele, file=paste0(PARAM$folder.data, "SNV.data/", tax.id, ".allele_data.filtered.RData"))
save(data.pos, mat.pos, file=paste0(PARAM$folder.data, "SNV.data/", tax.id, ".SNV_data.filtered.RData"))
save(curr.data.taxa, file=paste0(PARAM$folder.data, "SNV.data/", tax.id, ".taxa_data.RData"))
#Beautify output
writeLines("\n")
################################################################################
################################################################################







q(save="no")


