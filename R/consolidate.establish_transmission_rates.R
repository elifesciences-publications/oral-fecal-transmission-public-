#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Consolidate per-taxon jobs on rates of oral <-> gut transmission.
#
#=> iterate through taxa
#=> collect data in large frame for further analyses
#=> save and export individual datafiles in R format
#
#2018-06-28
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
PARAM$folder.results_transmission_rates <- paste0(PARAM$folder.results, "establish.transmission_rates/")
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
data.rates <- data.rates.bg <- data.frame()

#Iterate through taxa and process/collect results
for (tax.id in as.character(data.taxa.abd$Tax_ID)){
  #Get current taxon name
  taxon <- as.character(data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"])
  print.taxon <- gsub(" ", "_", taxon)
  print.taxon <- gsub("/", "_", print.taxon)
  
  print(paste(date(), "=> Processing", tax.id, taxon))
  
  #Get current file names
  file.data_transmission_rates <- paste0(PARAM$folder.results_transmission_rates, tax.id, ".", print.taxon, ".transmission_rates.RData")
  
  #Skip if files for taxon don't exist
  if (! file.exists(file.data_transmission_rates)) {next()}
  
  #Load data
  load(file.data_transmission_rates)
  
  if (nrow(sub.delta_t) == 0) {next()}
  if (is.null(sub.delta_t.bg)) {next()}
  if (nrow(sub.delta_t.bg) == 0) {next()}
  if (! "bg.coupled" %in% sub.delta_t.bg$delta_t.bin) {print(paste(tax.id, taxon, "is still in old version. Skipping.")); next()}
  
  #Add current sample pair/triplet identifier
  sub.delta_t$triplet <- paste(as.character(sub.delta_t$subject), as.character(sub.delta_t$timepoint.1), as.character(sub.delta_t$timepoint.2), as.character(sub.delta_t$site.1), sep=".")
  
  #Add taxon information
  sub.delta_t$Tax_ID <- sub.delta_t.bg$Tax_ID <- tax.id
  sub.delta_t$Taxon <- sub.delta_t.bg$Taxon <- taxon
  
  #Drop all triplets with insufficient information
  sub.delta_t <- sub.delta_t[!is.na(sub.delta_t$source_sink) & sub.delta_t$cross.n.shared_obs.coupled >= 100 & !is.nan(sub.delta_t$turnover.transmission.raw) & !is.nan(sub.delta_t$turnover.transmission), ]
  if (nrow(sub.delta_t) == 0) {next()}
  
  #Preallocate
  sub.delta_t$size.bg_longit <- sub.delta_t$Z.bg_longit <- sub.delta_t$size.bg_coupled <- sub.delta_t$Z.bg_coupled <- NA
  
  #Loop through triplets to sort out the appropriate backgrounds
  curr.bg <- data.frame()
  for (triplet in as.character(sub.delta_t$triplet)) {
    curr.sub.delta_t <- sub.delta_t[sub.delta_t$triplet == triplet, ]
    
    #Get current longitudinal background
    bg.longit <- sub.delta_t.bg[
      sub.delta_t.bg$delta_t.bin == "bg.longitudinal" &
        sub.delta_t.bg$cross.n.shared_obs.coupled >= 100 &
        !is.nan(sub.delta_t.bg$turnover.transmission.raw) &
        !is.nan(sub.delta_t.bg$turnover.transmission) &
        as.character(sub.delta_t.bg$subject) == curr.sub.delta_t$subject[1] &
        as.character(sub.delta_t.bg$sample.1) == curr.sub.delta_t$sample.1[1],
      ]
    
    if (nrow(bg.longit) > 0) {
      bg.longit$triplet <- triplet
      curr.bg <- rbind(curr.bg, bg.longit)
      sub.delta_t[sub.delta_t$triplet == triplet, "size.bg_longit"] <- nrow(bg.longit)
      if (sum(!is.na(bg.longit$turnover.non_transmission.raw)) > 5) {
        sub.delta_t[sub.delta_t$triplet == triplet, "Z.bg_longit"] <- (curr.sub.delta_t$turnover.transmission.raw - mean(bg.longit$turnover.transmission.raw, na.rm=T)) / sd(bg.longit$turnover.transmission.raw, na.rm=T)
      }
    }
    
    #Get current coupled background
    bg.coupled <- sub.delta_t.bg[
      sub.delta_t.bg$delta_t.bin == "bg.coupled" &
        sub.delta_t.bg$cross.n.shared_obs.coupled >= 100 &
        !is.nan(sub.delta_t.bg$turnover.transmission.raw) &
        !is.nan(sub.delta_t.bg$turnover.transmission) &
        as.character(sub.delta_t.bg$subject) == curr.sub.delta_t$subject[1] &
        as.character(sub.delta_t.bg$sample.1) == curr.sub.delta_t$sample.1[1] &
        as.character(sub.delta_t.bg$sample.2) == curr.sub.delta_t$sample.2[1],
      ]
    
    if (nrow(bg.coupled) > 0) {
      bg.coupled$triplet <- triplet
      curr.bg <- rbind(curr.bg, bg.coupled)
      sub.delta_t[sub.delta_t$triplet == triplet, "size.bg_coupled"] <- nrow(bg.longit)
      if (sum(!is.na(bg.coupled$turnover.non_transmission.raw)) > 5) {
        sub.delta_t[sub.delta_t$triplet == triplet, "Z.bg_coupled"] <- (curr.sub.delta_t$turnover.transmission.raw - mean(bg.coupled$turnover.transmission.raw, na.rm=T)) / sd(bg.coupled$turnover.transmission.raw, na.rm=T)
      }
    }
  }
  
  #Append to global results collector
  data.rates <- rbind(data.rates, sub.delta_t)
  data.rates.bg <- rbind(data.rates.bg, curr.bg)
}

#Re-order taxa by phylogeny
data.rates$Taxon <- factor(data.rates$Taxon, levels=levels(data.taxa$Scientific_Name)[levels(data.taxa$Scientific_Name) %in% data.rates$Taxon])
data.rates.bg$Taxon <- factor(data.rates.bg$Taxon, levels=levels(data.taxa$Scientific_Name)[levels(data.taxa$Scientific_Name) %in% data.rates$Taxon])

#Re-order subjects, nested by
#=> cohort
#=> family/village
#=> group
subjects.ordered <- order(as.character(data.rates$cohort), as.character(data.rates$family))
subjects.sorted <- unique(as.character(data.rates$subject)[subjects.ordered])
data.rates$subject <- factor(as.character(data.rates$subject), levels=subjects.sorted)
data.rates.bg$subject <- factor(as.character(data.rates.bg$subject), levels=subjects.sorted)

#Fix levels for source_sink pairs
data.rates$source_sink <- factor(data.rates$source_sink, levels=c("saliva_to_stool", "stool_to_saliva"))
data.rates.bg$source_sink <- factor(data.rates.bg$source_sink, levels=c("saliva_to_stool", "stool_to_saliva"))

#Fix levels for background types.
data.rates$delta_t.bin <- "intra.indiv"
data.rates.bg$delta_t.bin <- factor(data.rates.bg$delta_t.bin, levels=c("bg.coupled", "bg.longitudinal"))

#Sub-set to complete cases
data.rates <- data.rates[complete.cases(data.rates[, c("Taxon", "source_sink")]), ]

#Save data
save(data.rates, data.rates.bg, file=paste0(PARAM$folder.results, "data.transmission_rates.RData"))
################################################################################
################################################################################



