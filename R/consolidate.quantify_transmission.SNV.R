#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Consolidate per-taxon jobs on oral-gut transmission quantification from SNV data.
#
#
#=> iterate through taxa
#=> correlate 
#=> collect data in large frame for further analyses
#=> save and export individual datafiles in R format
#
#2017-07-27
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("plyr", warn.conflicts=F, quietly=T)
library("ggplot2", warn.conflicts=F, quietly=T)
library("RColorBrewer", warn.conflicts=F, quietly=T)
library("reshape2", warn.conflicts=F, quietly=T)
################################################################################
################################################################################

################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.base <- gsub("R/", "", PARAM$folder.R)
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")
PARAM$folder.results_detection <- paste0(PARAM$folder.results, "detect.transmission/")
PARAM$folder.results_quantification <- paste0(PARAM$folder.results, "quantify.transmission/")
################################################################################
################################################################################

################################################################################
################################################################################
#Load sample data
load(paste0(PARAM$folder.parameters, "data.sample.RData"))

#Load timepoint data
load(paste0(PARAM$folder.parameters, "data.timepoint.RData"))

#Load taxa data
load(paste0(PARAM$folder.data, "data.taxa.abd_corr.RData"))
data.taxa.abd <- data.taxa
load(paste0(PARAM$folder.data, "data.taxa.SNV.RData"))
data.taxa.SNV <- data.taxa

#Load specI abundance data
load(paste0(PARAM$folder.data, "data.specI.RData"))

#Load coverage data
load(paste0(PARAM$folder.data, "data.coverage.RData"))

#Load taxonomy data
load(paste0(PARAM$folder.data, "data.taxonomy.RData"))
################################################################################
################################################################################

################################################################################
################################################################################
#Preallocate data.frame to collect results
data.quantify_transmission <- data.frame()

#Iterate through taxa and process/collect results
for (tax.id in as.character(data.taxa.abd$Tax_ID)){
  #Get current taxon name
  taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
  print.taxon <- gsub(" ", "_", taxon)
  print.taxon <- gsub("/", "_", print.taxon)
  
  #Get current file names
  file.data_quantify <- paste0(PARAM$folder.results_quantification, tax.id, ".", print.taxon, ".transmission_per_subject.RData")
  #file.data_detect <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".detect_transmission.RData")
  
  #Skip if files for taxon don't exist
  if (! file.exists(file.data_quantify)) {next()}
  #if (! file.exists(file.data_detect)) {next()}
  
  #Load data
  load(file.data_quantify)
  #load(file.data_detect)
  
  if (ncol(data.s_t.transmission) != 43) {next()}
  
  #Get transmission score from "detect.transmission" data
  #data.s_t.transmission$z.transmission <- data.s_t.oi[match(rownames(data.s_t.oi), data.s_t.transmission$subject_timepoint), "S.Seb.z.bg_full.all"]
  
  #Concatenate
  data.quantify_transmission <- rbind(data.quantify_transmission, data.s_t.transmission)
}

#Save
save(data.quantify_transmission, file=paste0(PARAM$folder.results, "data.quantify.transmission.RData"))
################################################################################
################################################################################



q(save="no")

