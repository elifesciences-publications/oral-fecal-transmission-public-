#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Consolidate SNV data across taxa for further analyses
#
#=> load preliminary taxa data
#=> load taxa data per taxon (from parallel processing)
#=> consolidate and store
#
#2017-04-26
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

#Set relevant folder names
PARAM$folder.base <- paste0(getwd(), "/")
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
################################################################################
################################################################################


################################################################################
################################################################################
#Load preliminary taxon data
load(paste0(PARAM$folder.data, "data.taxa.abd_corr.RData"))

#Preallocate
data.taxa.add <- data.frame()

#Iterate through taxa and append additional taxon data
for (tax.id in as.character(data.taxa$Tax_ID)) {
  load(paste0(PARAM$folder.SNV_data, tax.id, ".taxa_data.RData"))
  data.taxa.add <- rbind(data.taxa.add, curr.data.taxa)
}

#Rename novel columns (fixing a bug from earlier runs)
# colnames(data.taxa.add) <- c(
#   "SNV_data.available",
#   "SNV_data.oral_gut",
#   "n.pos.total",
#   "n.pos.oral",
#   "n.snv.oral",
#   "n.pos.gut",
#   "n.snv.gut",
#   "n.pos.shared_habitat",
#   "n.snv.shared_habitat",
#   "n.subjects_timepoints.oi",
#   "n.pos.shared_subject_timepoint",
#   "n.snv.shared_subject_timepoint"
# )

#Append novel taxa data
data.taxa <- cbind(data.taxa, data.taxa.add)

#Export to tsv
write.table(data.taxa, file=paste0(PARAM$folder.data, "data.taxa.SNV.tsv"), row.names=F, sep="\t", quote=F)
#Save in R format
save(data.taxa, file=paste0(PARAM$folder.data, "data.taxa.SNV.RData"))
################################################################################
################################################################################





q()

