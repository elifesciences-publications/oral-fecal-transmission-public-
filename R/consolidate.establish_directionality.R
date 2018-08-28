#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Consolidate per-taxon jobs on directionality detection of oral <-> gut transmission.
#
#=> iterate through taxa
#=> collect data in large frame for further analyses
#=> save and export individual datafiles in R format
#
#2017-10-25
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("plyr", warn.conflicts=F, quietly=T)
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
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")
PARAM$folder.results_directionality <- paste0(PARAM$folder.results, "establish.directionality/")
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
data.directionality <- data.frame()

#Iterate through taxa and process/collect results
for (tax.id in as.character(data.taxa.abd$Tax_ID)){
  #Get current taxon name
  taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
  print.taxon <- gsub(" ", "_", taxon)
  print.taxon <- gsub("/", "_", print.taxon)
  
  #Get current file names
  file.data_directionality <- paste0(PARAM$folder.results_directionality, tax.id, ".", print.taxon, ".directionality.fisher_test.RData")
  
  #Skip if files for taxon don't exist
  if (! file.exists(file.data_directionality)) {next()}
  
  #Load data
  load(file.data_directionality)
  
  #Preallocate current results collection structure
  curr.data <- data.frame()
  
  #Iterate through candidate "sink" sites
  for (sink.site in c("saliva", "dental_plaque", "stool")) {
    curr.sub.delta_t <- sub.delta_t[sub.delta_t$site.1 == sink.site, ]
    
    #Iterate through candidate "source" sites
    for (source.site in c("saliva", "dental_plaque", "stool")) {
      if (source.site == sink.site) {next()}
      if (! paste(source.site, "p", sep=".") %in% colnames(curr.sub.delta_t)) {next()}
      curr.sub.delta_t.sub <- curr.sub.delta_t[!is.na(curr.sub.delta_t[, paste(source.site, "p", sep=".")]), ]
      if (nrow(curr.sub.delta_t.sub) == 0) {next()}
      
      curr.data <- rbind(curr.data, cbind(
        data.frame(Tax_ID=rep.int(tax.id, nrow(curr.sub.delta_t.sub)), Taxon=rep.int(taxon, nrow(curr.sub.delta_t.sub)), stringsAsFactors=F),
        curr.sub.delta_t.sub[, 1:3],
        group = data.timepoint[paste(curr.sub.delta_t.sub$subject, curr.sub.delta_t.sub$timepoint.1, sep="."), "group"],
        data.frame(source.site = rep.int(source.site, nrow(curr.sub.delta_t.sub)), sink.site=sink.site, stringsAsFactors=F),
        source_sink = paste(source.site, sink.site, sep="_to_"),
        curr.sub.delta_t.sub[, 6:15],
        source_sink.n.shared_obs = curr.sub.delta_t.sub[, paste0(source.site, ".n.shared_obs")],
        source_sink.n.gained.coupled = curr.sub.delta_t.sub[, paste0(source.site, ".n.gain.coupled")],
        source_sink.n.gained.not_coupled = curr.sub.delta_t.sub[, paste0(source.site, ".n.gain.not_coupled")],
        source_sink.n.not_gained.coupled = curr.sub.delta_t.sub[, paste0(source.site, ".n.not_gained.coupled")],
        source_sink.odds_ratio = curr.sub.delta_t.sub[, paste0(source.site, ".odds_ratio")],
        source_sink.log2_or = log2(curr.sub.delta_t.sub[, paste0(source.site, ".odds_ratio")]),
        source_sink.p.raw = curr.sub.delta_t.sub[, paste0(source.site, ".p")]
      ))
    }
  }
  
  #Adjust p-values for current taxon
  curr.data$source_sink.p.adjusted <- p.adjust(curr.data$source_sink.p.raw, method = "BH")
  
  #Append to global results collector
  data.directionality <- rbind(data.directionality, curr.data)
}

#Re-order taxa by phylogeny
data.directionality$Taxon <- factor(data.directionality$Taxon, levels=levels(data.taxa$Scientific_Name)[levels(data.taxa$Scientific_Name) %in% data.directionality$Taxon])

#Re-set -Inf and +Inf odds ratios to global minimum/maximum
data.directionality$source_sink.log2_or[data.directionality$source_sink.log2_or == Inf] <- max(data.directionality$source_sink.log2_or[data.directionality$source_sink.log2_or != Inf])
data.directionality$source_sink.log2_or[data.directionality$source_sink.log2_or == -Inf] <- min(data.directionality$source_sink.log2_or[data.directionality$source_sink.log2_or != -Inf])

#Annotate significance levels
data.directionality$significance <- "FALSE"
data.directionality$significance[data.directionality$source_sink.p.adjusted < 0.05] <- "<0.05"
data.directionality$significance[data.directionality$source_sink.p.adjusted < 0.01] <- "<0.01"
data.directionality$significance[data.directionality$source_sink.p.adjusted < 0.001] <- "<0.001"
data.directionality$significance <- factor(data.directionality$significance, levels=c("FALSE", "<0.05", "<0.01", "<0.001"))

#Re-order subjects, nested by
#=> cohort
#=> family/village
#=> group
subjects.ordered <- order(as.character(data.directionality$cohort), as.character(data.directionality$family), as.character(data.directionality$group))
subjects.sorted <- unique(as.character(data.directionality$subject)[subjects.ordered])
data.directionality$subject <- factor(as.character(data.directionality$subject), levels=subjects.sorted)

#Fix levels for source_sink pairs
data.directionality$source_sink <- factor(data.directionality$source_sink, levels=c("saliva_to_stool", "stool_to_saliva", "dental_plaque_to_stool", "stool_to_dental_plaque", "saliva_to_dental_plaque", "dental_plaque_to_saliva"))

#Save data
save(data.directionality, file=paste0(PARAM$folder.results, "data.directionality.RData"))
################################################################################
################################################################################



