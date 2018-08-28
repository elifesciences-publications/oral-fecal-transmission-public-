#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Consolidate per-taxon jobs on longitudinal stability of allele profiles in the various sites.
#
#=> iterate through taxa
#=> collect data in large frame for further analyses
#=> save and export individual datafiles in R format
#
#2017-11-08
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("plyr", warn.conflicts=F, quietly=T)
library("reshape2", warn.conflicts=F, quietly=T)
library("ggplot2", warn.conflicts=F, quietly=T)
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
PARAM$folder.results_stability <- paste0(PARAM$folder.results, "longitudinal.stability/")
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
data.stability <- data.frame()

#Iterate through taxa and process/collect results
for (tax.id in as.character(data.taxa.abd$Tax_ID)){
  #Get current taxon name
  taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
  print.taxon <- gsub(" ", "_", taxon)
  print.taxon <- gsub("/", "_", print.taxon)
  
  #Get current file names
  file.data_stability <- paste0(PARAM$folder.results_stability, tax.id, ".", print.taxon, ".longitudinal_allele_profiles.RData")
  
  #Skip if files for taxon don't exist
  if (! file.exists(file.data_stability)) {next()}
  
  #Load data
  load(file.data_stability)
  
  #Sub-select current results (only within-site comparisons are considered, as only these were indeed computed previously)
  curr.data <- results.allele_overlap[results.allele_overlap$site.1 == results.allele_overlap$site.2, ]
  
  #Iterate through candidate "sink" sites
  for (site in c("saliva", "dental_plaque", "stool")) {
    curr.data.sub <- curr.data[curr.data$site.1 == site, ]
    
    #Append to global results collection frame
    data.stability <- rbind(data.stability, cbind(
      data.frame(Tax_ID=rep.int(tax.id, nrow(curr.data.sub)), Taxon=rep.int(taxon, nrow(curr.data.sub))),
      curr.data.sub[, c(1:3,6:11)],
      data.frame(
        site = site,
        allele.overlap = curr.data.sub[, paste0(site, ".allele.overlap")],
        allele.overlap.bg_generic.mean = curr.data.sub[, paste0(site, ".allele.overlap.bg_generic.mean")],
        allele.overlap.bg_generic.sd = curr.data.sub[, paste0(site, ".allele.overlap.bg_generic.sd")]
      )
    ))
  }
  
  #Plot longitudinal profiles for current taxon
  #=> if bg data is available
  if (!exists("collect.allele_overlap")) {next()}
  plot.data <- cbind(data.stability[data.stability$Tax_ID == tax.id, c("cohort", "delta_t.bin", "site", "allele.overlap")], data.frame(type=rep.int("obs", sum(data.stability$Tax_ID == tax.id))))
  plot.data$cohort <- as.character(plot.data$cohort)
  for (coh in names(collect.allele_overlap)) {
    for (site in names(collect.allele_overlap[[coh]])) {
      plot.data <- rbind(plot.data, data.frame(
        cohort=coh,
        delta_t.bin=runif(length(collect.allele_overlap[[coh]][[site]]), min=0, max=max(data.delta_t$delta_t.bin, na.rm=T)),
        site=site,
        allele.overlap=collect.allele_overlap[[coh]][[site]],
        type="bg",
        stringsAsFactors = F
      ))
    }
  }
  plot.data$cohort[plot.data$cohort %in% c("MM", "FR-CRC")] <- "MM.FR-CRC"
  #Generate plot
  curr.plot <- ggplot(plot.data, aes(x=delta_t.bin, y=allele.overlap, colour=interaction(type, site), fill=interaction(type, site), alpha=type)) +
    geom_point() +
    geom_smooth(method="lm", alpha=0.6) +
    scale_color_manual(values=c("#1f78b4", "#bdbdbd", "#33a02c", "#bdbdbd", "#f1a340", "#bdbdbd")) +
    scale_fill_manual(values=c("#1f78b4", "#bdbdbd", "#33a02c", "#bdbdbd", "#f1a340", "#bdbdbd")) +
    scale_alpha_discrete(range=c(1,0.4)) +
    scale_y_continuous(limits=c(0.5,1)) +
    scale_x_continuous(trans = "sqrt") +
    ggtitle(taxon) +
    facet_grid(cohort ~ site) +
    theme_bw()
  ggsave(curr.plot, width=15, height=15, filename=paste0(PARAM$folder.results_stability, print.taxon, ".longitudinal_stability.pdf"), useDingbats=F)
  #Tidy up
  rm(collect.allele_overlap)
}

#Re-order taxa by phylogeny
data.stability$Taxon <- factor(data.stability$Taxon, levels=levels(data.taxa$Scientific_Name)[levels(data.taxa$Scientific_Name) %in% data.stability$Taxon])

#Z transform allele overlap
data.stability$allele.overlap.Z <- (data.stability$allele.overlap - data.stability$allele.overlap.bg_generic.mean) / (data.stability$allele.overlap.bg_generic.sd)

#Generate global overview plot, by taxon
curr.plot <- ggplot(data.stability, aes(x=Taxon, y=allele.overlap.Z, colour=delta_t.bin)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="black") +
  geom_jitter(alpha=0.6, width=0.2) +
  scale_color_distiller(type="div", palette="PRGn") +
  scale_x_discrete(limits = rev(levels(data.stability$Taxon))) +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1.645) + geom_hline(yintercept=2.33) + geom_hline(yintercept=3.0902) +
  geom_hline(yintercept=-1.645) + geom_hline(yintercept=-2.33) + geom_hline(yintercept=-3.0902) +
  ylim(c(-5, 5)) +
  ylab("Allele Overlap within Subject\n(Z transformed)") +
  coord_flip() +
  facet_grid(. ~ site) +
  theme_bw()
ggsave(curr.plot, width=25, height=40, filename=paste0(PARAM$folder.results, "longitudinal_stability.per_taxon.pdf"), useDingbats=F)

#Re-order subjects, nested by
#=> cohort
#=> family/village
subjects.ordered <- order(as.character(data.stability$cohort), as.character(data.stability$family))
subjects.sorted <- unique(as.character(data.stability$subject)[subjects.ordered])
data.stability$subject <- factor(as.character(data.stability$subject), levels=subjects.sorted)

#Generate overview plot, by subject
plot.data <- data.stability
plot.data$subject.delta_t.bin <- droplevels(interaction(plot.data$subject, plot.data$delta_t.bin))
subj.dt.ordered <- order(as.character(plot.data$cohort), as.character(plot.data$family), as.character(plot.data$subject))
subj.dt.sorted <- unique(as.character(plot.data$subject.delta_t.bin)[subj.dt.ordered])
plot.data$subject.delta_t.bin <- factor(as.character(plot.data$subject.delta_t.bin), levels=subj.dt.sorted)
curr.plot <- ggplot(plot.data, aes(x=subject.delta_t.bin, y=allele.overlap.Z, colour=delta_t.bin)) +
  geom_boxplot(outlier.colour=NA, fill=NA, colour="black") +
  geom_jitter(alpha=0.3, width=0.2, size=1) +
  scale_color_distiller(type="div", palette="PRGn") +
  scale_x_discrete(limits = rev(levels(plot.data$subject.delta_t.bin))) +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1.645) + geom_hline(yintercept=2.33) + geom_hline(yintercept=3.0902) +
  geom_hline(yintercept=-1.645) + geom_hline(yintercept=-2.33) + geom_hline(yintercept=-3.0902) +
  ylim(c(-5, 5)) +
  ylab("Allele Overlap within Subject\n(Z transformed)") +
  coord_flip() +
  facet_grid(. ~ site) +
  theme_bw()
ggsave(curr.plot, width=25, height=30, filename=paste0(PARAM$folder.results, "longitudinal_stability.per_subject.pdf"), useDingbats=F)

#Save data
save(data.stability, file=paste0(PARAM$folder.results, "data.stability.RData"))
################################################################################
################################################################################



