#!~/.nix-profile/bin/Rscript
################################################################################
#Project: Microbial Transmission on the Oral-Gut Axis
#
#Consolidate per-taxon jobs on oral-gut transmission detection from SNV data.
#
#
#=> iterate through taxa
#=> plot background distributions of scores (per cohort)
#=> plot Z scores
#=> collect data in large frame for further analyses
#=> save and export individual datafiles in R format
#
#2017-05-02
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
# Load Packages
library("Matrix", warn.conflicts=F, quietly=T)
library("plyr", warn.conflicts=F, quietly=T)
library("ggplot2", warn.conflicts=F, quietly=T)
library("RColorBrewer", warn.conflicts=F, quietly=T)
library("corrplot", warn.conflicts=F, quietly=T)
library("reshape2", warn.conflicts=F, quietly=T)
################################################################################
################################################################################

################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.base <- gsub("R/", "", PARAM$folder.R)
PARAM$folder.data <- paste0(PARAM$folder.base, "data/")
PARAM$folder.SNV_data <- paste0(PARAM$folder.data, "SNV.data/")
PARAM$folder.parameters <- paste0(PARAM$folder.base, "parameters/")
PARAM$folder.results <- paste0(PARAM$folder.base, "results/")
PARAM$folder.results_detection <- paste0(PARAM$folder.results, "detect.transmission/")

#Thresholds to filter detection results
PARAM$thresh.min_shared_obs <- 20

#Predefine list of taxa with per-subject computations
PARAM$taxa.per_subject <- c(
  866776,
  879309,
  388919,
  1051074,
  596315,
  706433,
  592028,
  944560,
  411466,
  680646,
  246201,
  585202,
  585204,
  365659,
  1035189,
  997830,
  1005705,
  871237,
  864570,
  888049,
  546273,
  655813,
  888833,
  467705,
  889201,
  904306,
  435842,
  638301,
  562982,
  562983,
  521095,
  469621,
  1048332,
  862965,
  479436,
  888052,
  760570,
  904296,
  585501,
  575593,
  887325,
  944565,
  411465,
  767100,
  706434,
  546271,
  563038,
  1054460,
  487214,
  888746,
  927666,
  768726,
  1000588,
  563037,
  864567,
  904294,
  706437,
  742820,
  1035184,
  210007,
  712961,
  324831,
  401473,
  762963,
  633147,
  1028805,
  935897,
  1035188,
  888057,
  888743,
  767031,
  908937,
  575611,
  717959,
  626523,
  592010,
  888060,
  684738,
  626369,
  562981,
  546270,
  649764,
  431947,
  203275,
  553174,
  553171,
  537011,
  866771,
  997352,
  649761,
  585543,
  619693,
  862515,
  997353,
  888832,
  634176,
  563008,
  626522,
  575615,
  688246,
  585502,
  457405,
  1028803,
  679190,
  873533,
  679189,
  702438,
  575614,
  264731,
  449673,
  596327,
  553175,
  908612
)

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
data.detect_transmission <- data.frame()
data.taxa.transmission <- list()

#Define combinations of body sites to explore
site.combinations <- c("saliva_gut", "plaque_gut", "plaque_saliva")
for (coupled.sites in site.combinations) {data.taxa.transmission[[coupled.sites]] <- data.frame()}

#Define function to Z transform a distribution
z.transform <- function(x) {(x - mean(x, na.rm=T)) / sd(x, na.rm=T)}

#Get list of unique subjects
subjects <- as.character(unique(data.sample$subject))

#Iterate through taxa and process/collect results
for (tax.id in as.character(data.taxa.abd$Tax_ID)){
  #Get current taxon name
  taxon <- data.taxonomy[data.taxonomy$Tax_ID == tax.id, "Scientific_Name"]
  print.taxon <- gsub(" ", "_", taxon)
  print.taxon <- gsub("/", "_", print.taxon)
  
  #Iterate through habitat pairs
  for (coupled.sites in site.combinations) {
    #Check whether we have to consolidate per-subject files first
    #=> otherwise load and process global dataset
    if ((coupled.sites == "saliva_gut" & tax.id %in% as.character(PARAM$taxa.per_subject)) | (coupled.sites == "plaque_gut" & tax.id %in% as.character(PARAM$taxa.per_subject.plaque_gut))) {
      #Preallocate per-subject results collectors
      ps.data.s_t.oi <- data.frame()
      ps.data.bg.full <- ps.data.bg.full.indiv <- ps.data.bg.cohort <- ps.data.bg.cohort.indiv <- data.frame()
      ps.bg.full <- ps.bg.full.indiv <- ps.bg.cohort <- ps.bg.cohort.indiv <- list()
      
      #Iterate through subjects to load and process file
      for (subj in subjects) {
        curr.file.detect <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".detect_transmission.", coupled.sites, ".", subj, ".RData")
        curr.file.bg <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".background_scores.", coupled.sites, ".", subj, ".RData")
        
        #Skip if data is not available for current subject
        if (!file.exists(curr.file.detect) | !file.exists(curr.file.bg)) {next()}
        
        #Load current data
        load(curr.file.detect)
        load(curr.file.bg)
        
        if (nrow(data.s_t.oi) == 0) {writeLines(paste(subj, "is empty.")); next()}
        
        #Append
        ps.data.s_t.oi <- rbind(ps.data.s_t.oi, data.s_t.oi)
        #Per-background data
        ps.data.bg.full <- rbind(ps.data.bg.full, data.bg.full)
        ps.data.bg.full.indiv <- rbind(ps.data.bg.full.indiv, data.bg.full.indiv)
        ps.data.bg.cohort <- rbind(ps.data.bg.cohort, data.bg.cohort)
        ps.data.bg.cohort.indiv <- rbind(ps.data.bg.cohort.indiv, data.bg.cohort.indiv)
        #Background scores
        for (i in 1:length(S.bg.cohort)) {
          if (coupled.sites == "saliva_gut") {
            ps.bg.full[[names(S.bg.cohort)[i]]] <- S.bg.full[[i]]
            ps.bg.full.indiv[[names(S.bg.cohort)[i]]] <- S.bg.full.indiv[[i]]
          }
          ps.bg.cohort[[names(S.bg.cohort)[i]]] <- S.bg.cohort[[i]]
          ps.bg.cohort.indiv[[names(S.bg.cohort)[i]]] <- S.bg.cohort.indiv[[i]]
        }
      }
      
      #Pass for further use
      data.s_t.oi <- ps.data.s_t.oi
      data.bg.full <- ps.data.bg.full
      data.bg.full.indiv <- ps.data.bg.full.indiv
      data.bg.cohort <- ps.data.bg.cohort
      data.bg.cohort.indiv <- ps.data.bg.cohort.indiv
      S.bg.full <- ps.bg.full
      S.bg.full.indiv <- ps.bg.full.indiv
      S.bg.cohort <- ps.bg.cohort
      S.bg.cohort.indiv <- ps.bg.cohort.indiv
    } else {
      #Get current file names
      file.data_transmission <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".detect_transmission.", coupled.sites, ".RData")
      file.data_bg <- paste0(PARAM$folder.results_detection, tax.id, ".", print.taxon, ".background_scores.", coupled.sites, ".RData")
      
      #Skip if files for taxon don't exist
      if (! file.exists(file.data_bg)) {next()}
      if (file.info(file.data_bg)$size == 0) {next()}
      
      #Load data
      load(file.data_transmission)
      load(file.data_bg)
    }
    
    #Remove all intra-individual pairs with too few shared observations
    data.s_t.oi <- data.s_t.oi[data.s_t.oi$n.obs.shared > PARAM$thresh.min_shared_obs & !is.na(data.s_t.oi$n.obs.shared), ]
    #Skip if no data remains
    if (nrow(data.s_t.oi) == 0) {next()}
    
    #Skip if too few points exist (in full bg) to properly estimate background
    if (nrow(data.bg.cohort) == 0) {next()}
    if (median(data.bg.cohort$n.comparisons) < 1000) {next()}
    
    ###################################
    #Plot background score distributions
    #=> plot densities for non-transformed scores, for full and per-cohort-backgrounds
    ###################################
    plot.bg <- data.frame()
    #Extract scores
    if (coupled.sites == "saliva_gut") {
      bg.full <- unlist(S.bg.full)
      bg.full.indiv <- unlist(S.bg.full.indiv)
      plot.bg <- rbind(
        plot.bg,
        data.frame(Background="full", Score=sample(bg.full, min(10000, length(bg.full)))),
        data.frame(Background="full.indiv", Score=sample(bg.full.indiv, min(10000, length(bg.full.indiv))))
      )
    }
    #Extract per-cohort scores
    data.coh.bg <- data.timepoint
    data.coh.bg$cohort[data.coh.bg$cohort == "MM"] <- "FR-CRC"
    for (coh in unique(as.character(data.coh.bg$cohort))) {
      bg.cohort <- unlist(S.bg.cohort[rownames(data.coh.bg)[data.coh.bg$cohort == coh]])
      bg.cohort.indiv <- unlist(S.bg.cohort.indiv[rownames(data.coh.bg)[data.coh.bg$cohort == coh]])
      if (length(bg.cohort.indiv) == 0) {next()}
      if (length(bg.cohort) == 0) {next()}
      plot.bg <- rbind(
        plot.bg,
        data.frame(Background=coh, Score=sample(bg.cohort, min(10000, length(bg.cohort)))),
        data.frame(Background=paste0(coh, ".indiv"), Score=sample(bg.cohort.indiv, min(10000, length(bg.cohort.indiv))))
      )
    }
    #Generate plot
    curr.plot <- ggplot(plot.bg, aes(x=Score, colour=Background)) +
      stat_density(geom="line", position="identity") +
      ggtitle(paste("Raw Background Score Distribution for", taxon)) +
      scale_color_brewer(type="qual", palette="Paired") +
      theme_bw()
    ggsave(curr.plot, width=10, height=10, filename=paste0(PARAM$folder.results_detection, "bg_scores/", tax.id, ".", print.taxon,  ".bg_scores.", coupled.sites, ".pdf"), useDingbats=F)
    ###################################
    
    ###################################
    #Skip if too few subjects were score-able
    if (sum(!is.na(data.s_t.oi$S.z.bg_cohort.all)) <= 5) {next()}
    ###################################
    
    ###################################
    #Fix column names in data.s_t.oi
    colnames(data.s_t.oi) <- c("subject", "timepoint", "cohort", "extraction", "group", "sample.oral", 
                               "sample.gut", "available.oral", "available.gut", "available.plaque", 
                               "sample.plaque", "paired.oral_gut", "paired.saliva_gut", "paired.plaque_gut", 
                               "paired.saliva_plaque", "family", "hor_cov.site_1", "hor_cov.site_2", 
                               "ver_cov.site_1", "ver_cov.site_2", "abd.site_1", "abd.site_2", 
                               "n.obs", "n.allele", "n.obs.site_1", "n.allele.site_1", "n.obs.site_2", 
                               "n.allele.site_2", "n.obs.shared", "n.allele.shared", "jac.shared", 
                               "S.bg_full.agree", "S.bg_full.disagree", "S.bg_full.max_theo", 
                               "S.bg_cohort.agree", "S.bg_cohort.disagree", "S.bg_cohort.max_theo", 
                               "S.bg_full", "S.bg_cohort", "S.z.bg_full.all", "S.z.bg_full.indiv", 
                               "S.z.bg_cohort.all", "S.z.bg_cohort.indiv")
    
    ###################################
    
    ###################################
    #Process data per subject.timepoint
    #=> z-transform
    #=> per-subject one-sample t.test p values (intra-individual score vs subject-specific background)
    #=> relative rank vs subject-specific background
    ###################################
    for (s_t in rownames(data.s_t.oi)) {
      if (coupled.sites == "saliva_gut") {
        #Z-transform scores
        data.s_t.oi[s_t, "S.z.bg_full.all"] <- (data.s_t.oi[s_t, "S.bg_full"] - mean(S.bg.full[[s_t]], na.rm=T)) / sd(S.bg.full[[s_t]], na.rm=T)
        data.s_t.oi[s_t, "S.z.bg_full.indiv"] <- (data.s_t.oi[s_t, "S.bg_full"] - mean(S.bg.full.indiv[[s_t]], na.rm=T)) / sd(S.bg.full.indiv[[s_t]], na.rm=T)
        
        #Get current observation's relative rank against the subject-specific background
        data.s_t.oi[s_t, "S.rank.bg_full.all"] <- sum(S.bg.full[[s_t]] > data.s_t.oi[s_t, "S.bg_full"]) / length(na.omit(S.bg.full[[s_t]]))
        data.s_t.oi[s_t, "S.rank.bg_full.indiv"] <- sum(S.bg.full.indiv[[s_t]] > data.s_t.oi[s_t, "S.bg_full"]) / length(na.omit(S.bg.full.indiv[[s_t]]))
        
        #Get per-subject t-test p values
        if (is.na(data.s_t.oi[s_t, "S.bg_full"])) {
          data.s_t.oi[s_t, "S.t_test.bg_full.all.p.raw"] <- data.s_t.oi[s_t, "S.t_test.bg_full.indiv.p.raw"] <- NA
        } else {
          data.s_t.oi[s_t, "S.t_test.bg_full.all.p.raw"] <- t.test(S.bg.full[[s_t]], mu=data.s_t.oi[s_t, "S.bg_full"], alternative = "less")$p.value
          data.s_t.oi[s_t, "S.t_test.bg_full.indiv.p.raw"] <- t.test(S.bg.full.indiv[[s_t]], mu=data.s_t.oi[s_t, "S.bg_full"], alternative = "less")$p.value
        }
      } else {
        data.s_t.oi[s_t, "S.rank.bg_full.all"] <- data.s_t.oi[s_t, "S.rank.bg_full.indiv"] <- data.s_t.oi[s_t, "S.t_test.bg_full.all.p.raw"] <- data.s_t.oi[s_t, "S.t_test.bg_full.indiv.p.raw"] <- NA
      }
      #Repeat, for per-cohort backgrounds
      data.s_t.oi[s_t, "S.z.bg_cohort.all"] <- (data.s_t.oi[s_t, "S.bg_cohort"] - mean(S.bg.cohort[[s_t]], na.rm=T)) / sd(S.bg.cohort[[s_t]], na.rm=T)
      data.s_t.oi[s_t, "S.z.bg_cohort.indiv"] <- (data.s_t.oi[s_t, "S.bg_cohort"] - mean(S.bg.cohort.indiv[[s_t]], na.rm=T)) / sd(S.bg.cohort.indiv[[s_t]], na.rm=T)
      data.s_t.oi[s_t, "S.rank.bg_cohort.all"] <- sum(S.bg.cohort[[s_t]] > data.s_t.oi[s_t, "S.bg_cohort"]) / length(na.omit(S.bg.cohort[[s_t]]))
      data.s_t.oi[s_t, "S.rank.bg_cohort.indiv"] <- sum(S.bg.cohort.indiv[[s_t]] > data.s_t.oi[s_t, "S.bg_cohort"]) / length(na.omit(S.bg.cohort.indiv[[s_t]]))
      if (is.na(data.s_t.oi[s_t, "S.bg_cohort"])) {
        data.s_t.oi[s_t, "S.t_test.bg_cohort.all.p.raw"] <- data.s_t.oi[s_t, "S.t_test.bg_cohort.indiv.p.raw"] <- NA
      } else {
        data.s_t.oi[s_t, "S.t_test.bg_cohort.all.p.raw"] <- t.test(S.bg.cohort[[s_t]], mu=data.s_t.oi[s_t, "S.bg_cohort"], alternative = "less")$p.value
        data.s_t.oi[s_t, "S.t_test.bg_cohort.indiv.p.raw"] <- t.test(S.bg.cohort.indiv[[s_t]], mu=data.s_t.oi[s_t, "S.bg_cohort"], alternative = "less")$p.value
      }
      
    }
    ###################################
    
    ###################################
    #Z transform background data
    if (coupled.sites == "saliva_gut") {
      #Z transform background data
      S.bg.full.Z <- unlist(lapply(S.bg.full[rownames(data.s_t.oi)], z.transform))
      S.bg.full.indiv.Z <- unlist(lapply(S.bg.full.indiv[rownames(data.s_t.oi)], z.transform))
    }
    S.bg.cohort.Z <- unlist(lapply(S.bg.cohort[rownames(data.s_t.oi)], z.transform))
    S.bg.cohort.indiv.Z <- unlist(lapply(S.bg.cohort.indiv[rownames(data.s_t.oi)], z.transform))
    ###################################
    
    ###################################
    #Generate correlogram for relevant parameters
    if (coupled.sites == "saliva_gut") {
      curr.cor <- cor(data.s_t.oi[, c(
        "hor_cov.site_1", "hor_cov.site_2",
        "ver_cov.site_1", "ver_cov.site_2",
        "abd.site_1", "abd.site_2",
        "n.obs.shared", "n.allele.shared",
        "jac.shared",
        "S.z.bg_full.indiv", "S.z.bg_full.all",
        "S.z.bg_cohort.indiv", "S.z.bg_cohort.all"
      )], method="spearman", use="na.or.complete")
    } else {
      curr.cor <- cor(data.s_t.oi[, c(
        "hor_cov.site_1", "hor_cov.site_2",
        "ver_cov.site_1", "ver_cov.site_2",
        "abd.site_1", "abd.site_2",
        "n.obs.shared", "n.allele.shared",
        "jac.shared",
        "S.z.bg_cohort.indiv", "S.z.bg_cohort.all"
      )], method="spearman", use="na.or.complete")
    }
    if (! all(is.na(curr.cor))) {
      pdf(file=paste0(PARAM$folder.results_detection, "diagnostic_corrplots/", tax.id, ".", print.taxon, ".", coupled.sites, ".corrplot.pdf"), width=10, height=10, useDingbats=F)
      curr.corrplot <- corrplot(curr.cor, type="lower", method="circle", diag=F, title=paste(taxon, coupled.sites), mar=c(0,0,2,0))
      dev.off()
    }
    ###################################
    
    ###################################
    #Plot Z score distributions per subject.timepoint
    #=> gridded per background
    ###################################
    #Concatenate data for plotting
    if (coupled.sites == "saliva_gut") {
      df.melt <- melt(data.s_t.oi, id.vars="sample.oral", measure.vars=c("S.z.bg_full.all", "S.z.bg_full.indiv", "S.z.bg_cohort.all", "S.z.bg_cohort.indiv"), value.name="Z.score", variable.name="Test")
      df.melt$Type <- "within subject"
      df.melt <- rbind(
        df.melt,
        data.frame(sample.oral="dummy", Test="S.z.bg_full.all", Z.score=unlist(S.bg.full.Z), Type="bg between subjects"),
        data.frame(sample.oral="dummy", Test="S.z.bg_full.indiv", Z.score=unlist(S.bg.full.indiv.Z), Type="bg between subjects"),
        data.frame(sample.oral="dummy", Test="S.z.bg_cohort.all", Z.score=unlist(S.bg.cohort.Z), Type="bg between subjects"),
        data.frame(sample.oral="dummy", Test="S.z.bg_cohort.indiv", Z.score=unlist(S.bg.cohort.indiv.Z), Type="bg between subjects")
      )
    } else {
      df.melt <- melt(data.s_t.oi, id.vars="sample.oral", measure.vars=c("S.z.bg_cohort.all", "S.z.bg_cohort.indiv"), value.name="Z.score", variable.name="Test")
      df.melt$Type <- "within subject"
      df.melt <- rbind(
        df.melt,
        data.frame(sample.oral="dummy", Test="S.z.bg_cohort.all", Z.score=unlist(S.bg.cohort.Z), Type="bg between subjects"),
        data.frame(sample.oral="dummy", Test="S.z.bg_cohort.indiv", Z.score=unlist(S.bg.cohort.indiv.Z), Type="bg between subjects")
      )
    }
    #Generate plot
    curr.plot <- ggplot(droplevels(df.melt), aes(x=Z.score, colour=interaction(Type, Test))) +
      geom_density() +
      ggtitle(paste("Z Score Distribution Across Individuals for", coupled.sites, "in", taxon)) +
      scale_color_brewer(type = "qual", palette = "Paired") +
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = 1.645) + geom_vline(xintercept = 2.33) + geom_vline(xintercept = 3.09) +
      geom_vline(xintercept = -1.645) + geom_vline(xintercept = -2.33) + geom_vline(xintercept = -3.09) +
      facet_grid(Test ~ .) +
      theme_bw()
    ggsave(curr.plot, width=15, height=15, filename=paste0(PARAM$folder.results_detection, "z_distributions/", tax.id, ".", print.taxon, ".", coupled.sites, ".z_score_distribution.pdf"), useDingbats=F)
    ###################################
    
    ###################################
    #Get current data
    curr.frame <- cbind(
      data.frame(
        Tax_ID=rep.int(tax.id, nrow(data.s_t.oi)),
        Taxon=rep.int(taxon, nrow(data.s_t.oi)),
        coupled.sites=coupled.sites,
        stringsAsFactors=F
      ),
      data.s_t.oi
    )
    curr.frame <- curr.frame[!is.na(curr.frame$S.z.bg_cohort.all), ]
    ###################################
    
    ###################################
    #Calculate per-subject p values
    #=> is a subject's transmission score significant relative to the subject-specific background?
    curr.frame[, "S.z.bg_full.all.p.adj"] <- p.adjust(1 - pnorm(curr.frame$S.z.bg_full.all), method="BH")
    curr.frame[, "S.z.bg_full.indiv.p.adj"] <- p.adjust(1 - pnorm(curr.frame$S.z.bg_full.indiv), method="BH")
    curr.frame[, "S.z.bg_cohort.all.p.adj"] <- p.adjust(1 - pnorm(curr.frame$S.z.bg_cohort.all), method="BH")
    curr.frame[, "S.z.bg_cohort.indiv.p.adj"] <- p.adjust(1 - pnorm(curr.frame$S.z.bg_cohort.indiv), method="BH")
    ###################################
    
    ###################################
    #Test significance of transmission across subject.timepoints
    #curr.bg_full.all.p.raw <- wilcox.test(curr.frame$S.z.bg_full.all, unlist(S.bg.full.Z), paired=F)$p.value
    ###################################
    
    ###################################
    #Calculate "transmissibility"
    #=> number of subjects with abd >10^-6 and Z > median(bg(Z))
    ###################################
    #Get number of subjects with abd.oral, abd.gut or abd.oral&abd.gut >=10^-6 in >=1tp
    #subj.abd.site_1 <- ddply(curr.frame, "subject", function(x) {any(x$abd.site_1 > 10^-6)})$V1
    #subj.abd.site_2 <- ddply(curr.frame, "subject", function(x) {any(x$abd.site_2 > 10^-6)})$V1
    #subj.abd.paired <- subj.abd.site_1 & subj.abd.site_2
    subj.abd.site_1 <- curr.frame$abd.site_1 > 10^-6
    subj.abd.site_2 <- curr.frame$abd.site_2 > 10^-6
    subj.abd.paired <- subj.abd.site_1 & subj.abd.site_2
    #Get overall "transmissibility" per subject
    #tmp.transmissibility <- sapply(names(S.bg.cohort.Z), function(s_t) {curr.frame[s_t, "S.z.bg_cohort.all"] > median(S.bg.cohort.Z[[s_t]], na.rm=T)})
    #tmp.transmissibility <- sum(curr.frame[, "S.z.bg_cohort.all"] > median(S.bg.cohort.Z, na.rm=T)) / sum(curr.frame[, "S.z.bg_cohort.all"] < median(S.bg.cohort.Z, na.rm=T))
    #tmp.df <- data.frame(subject=curr.frame$subject, transmissibility=tmp.transmissibility)
    #subj.transmissibility <- ddply(tmp.df, "subject", function(x) {sum(x$transmissibility) > 0})$V1
    subj.above_median <- curr.frame[, "S.rank.bg_cohort.all"] < 0.5
    subj.below_median <- curr.frame[, "S.rank.bg_cohort.all"] > 0.5
    #Calculate transmissibility (as fraction of subjects)
    #transmissibility.site_1.abs <- sum(subj.abd.site_1 & subj.above_median) / 
    #transmissibility.site_2.abs <- sum(subj.abd.site_2 & subj.transmissibility)
    #transmissibility.paired.abs <- sum(subj.abd.paired & subj.transmissibility)
    transmissibility.raw <- sum(subj.above_median) / length(subj.above_median)
    transmissibility.rank_mean <- mean(curr.frame[, "S.rank.bg_cohort.all"], na.rm=T)
    transmissibility.site_1 <- log2(sum(subj.abd.site_1 & subj.above_median) / sum(subj.abd.site_1 & subj.below_median))
    transmissibility.site_2 <- log2(sum(subj.abd.site_2 & subj.above_median) / sum(subj.abd.site_2 & subj.below_median))
    transmissibility.paired <- log2(sum(subj.abd.paired & subj.above_median) / sum(subj.abd.paired & subj.below_median))
    ###################################
    
    ###################################
    #Plot transmission scores by cohort
    if (coupled.sites == "saliva_gut") {
      plot.melt <- melt(data.s_t.oi, id.vars=c("cohort", "group"), measure.vars=c("S.z.bg_full.all", "S.z.bg_full.indiv", "S.z.bg_cohort.all", "S.z.bg_cohort.indiv"), value.name="Z.score", variable.name="Test")
    } else {
      plot.melt <- melt(data.s_t.oi, id.vars=c("cohort", "group"), measure.vars=c("S.z.bg_cohort.all", "S.z.bg_cohort.indiv"), value.name="Z.score", variable.name="Test")
    }
    #Remove "small adenoma" category
    plot.melt <- plot.melt[! plot.melt$group %in% "small_adenoma", ]
    #Subset to include only controls
    plot.controls <- plot.melt[plot.melt$group == "control", ]
    if (length(unique(plot.melt$group)) > 1 & length(unique(plot.melt$Test)) > 1) {
      curr.plot <- ggplot(plot.controls, aes(x=cohort, y=Z.score, fill=cohort)) +
        geom_boxplot(alpha=0.7, outlier.colour=NA) +
        geom_point(position=position_jitterdodge(jitter.width=0.4)) +
        geom_hline(yintercept=0) + geom_hline(yintercept=1.645) + geom_hline(yintercept=2.326) + geom_hline(yintercept=3.090) +
        ggtitle(paste("Z Score Distributions per Cohort for", coupled.sites, "in", taxon)) +
        facet_grid(. ~ Test) +
        theme_bw() +
        theme(axis.text.x=element_text(angle = 90, hjust = 0))
      ggsave(curr.plot, width=15, height=10, filename=paste0(PARAM$folder.results_detection, "cohort/", tax.id, ".", print.taxon, ".", coupled.sites, ".by_cohort.pdf"), useDingbats=F)
      
      ###################################
      #Plot transmission scores per group
      #=> within cohort
      curr.plot <- ggplot(plot.melt, aes(x=group, y=Z.score, fill=group)) +
        geom_boxplot(alpha=0.7, outlier.colour=NA) +
        geom_point(position=position_jitterdodge(jitter.width=0.4), alpha=0.5) +
        geom_hline(yintercept=0) + geom_hline(yintercept=1.645) + geom_hline(yintercept=2.326) + geom_hline(yintercept=3.090) +
        scale_fill_brewer(type="qual", palette="Accent") +
        ggtitle(paste("Z Score Distributions per Cohort for", coupled.sites, "in", taxon)) +
        facet_grid(cohort ~ Test) +
        theme_bw() +
        theme(axis.text.x=element_text(angle = 90, hjust = 0))
      ggsave(curr.plot, width=15, height=20, filename=paste0(PARAM$folder.results_detection, "subject_group/", tax.id, ".", print.taxon, ".", coupled.sites, ".by_group.pdf"), useDingbats=F)
    }
    ###################################
    
    ###################################
    #Annotate current taxon
    curr.data.taxa <- cbind(
      data.taxa.abd[data.taxa.abd$Tax_ID == tax.id, ],
      data.frame(
        coupled.sites=coupled.sites,
        n.indiv.scored=nrow(curr.frame),
        #Full background, all comparisons
        Z.bg_full.all.mean = mean(curr.frame$S.z.bg_full.all, na.rm=T),
        Z.bg_full.all.sd = sd(curr.frame$S.z.bg_full.all, na.rm=T),
        Z.bg_full.all.median = median(curr.frame$S.z.bg_full.all, na.rm=T),
        #cor.Z_score.abd_oral = curr.cor["abd.oral", "S.z.bg_full.all"],
        #cor.Z_score.abd_gut = curr.cor["abd.gut", "S.z.bg_full.all"],
        Z.bg_full.all.wilcox_p.raw = ifelse(all(is.na(curr.frame$S.z.bg_full.all)), NA, wilcox.test(curr.frame$S.z.bg_full.all, unlist(S.bg.full.Z), paired=F, alternative="greater")$p.value),
        #Full background, comparisons against focal individual
        Z.bg_full.indiv.mean = mean(curr.frame$S.z.bg_full.indiv, na.rm=T),
        Z.bg_full.indiv.sd = sd(curr.frame$S.z.bg_full.indiv, na.rm=T),
        Z.bg_full.indiv.median = median(curr.frame$S.z.bg_full.indiv, na.rm=T),
        Z.bg_full.indiv.wilcox_p.raw = ifelse(all(is.na(curr.frame$S.z.bg_full.indiv)), NA, wilcox.test(curr.frame$S.z.bg_full.indiv, unlist(S.bg.full.indiv.Z), paired=F, alternative="greater")$p.value),
        #Full background, comparisons against focal individual
        Z.bg_cohort.all.mean = mean(curr.frame$S.z.bg_cohort.all, na.rm=T),
        Z.bg_cohort.all.sd = sd(curr.frame$S.z.bg_cohort.all, na.rm=T),
        Z.bg_cohort.all.median = median(curr.frame$S.z.bg_cohort.all, na.rm=T),
        Z.bg_cohort.all.wilcox_p.raw = ifelse(all(is.na(curr.frame$S.z.bg_cohort.all)), NA, wilcox.test(curr.frame$S.z.bg_cohort.all, unlist(S.bg.cohort.Z), paired=F, alternative="greater")$p.value),
        Z.bg_cohort.all.transmissibility = transmissibility.paired,
        Z.bg_cohort.all.transmissibility.raw = transmissibility.raw,
        Z.bg_cohort.all.transmissibility.rank_mean = transmissibility.rank_mean,
        #Z.bg_cohort.all.transmissibility.abs = transmissibility.paired.abs,
        Z.bg_cohort.all.transmissibility.site_1 = transmissibility.site_1,
        #Z.bg_cohort.all.transmissibility.site_1.abs = transmissibility.site_1.abs,
        Z.bg_cohort.all.transmissibility.site_2 = transmissibility.site_2,
        #Z.bg_cohort.all.transmissibility.site_2.abs = transmissibility.site_2.abs,
        #Full background, comparisons against focal individual
        Z.bg_cohort.indiv.mean = mean(curr.frame$S.z.bg_cohort.indiv, na.rm=T),
        Z.bg_cohort.indiv.sd = sd(curr.frame$S.z.bg_cohort.indiv, na.rm=T),
        Z.bg_cohort.indiv.median = median(curr.frame$S.z.bg_cohort.indiv, na.rm=T),
        Z.bg_cohort.indiv.wilcox_p.raw = ifelse(all(is.na(curr.frame$S.z.bg_cohort.indiv)), NA, wilcox.test(curr.frame$S.z.bg_cohort.indiv, unlist(S.bg.cohort.indiv.Z), paired=F)$p.value)
      )
    )
    data.taxa.transmission[[coupled.sites]] <- rbind(data.taxa.transmission[[coupled.sites]], curr.data.taxa)
    ###################################
    
    ###################################
    #Append current data to filtered results collector frame
    data.detect_transmission <- rbind(data.detect_transmission, curr.frame)
    ###################################
    
    ###################################
    #Report
    writeLines(paste(date(), "=> done with", taxon, coupled.sites))
    ###################################
  }
}
################################################################################
################################################################################


################################################################################
################################################################################
#Post-process transmission detection data
################################################################################
################################################################################
#Remove all intra-individual comparisons with too few shared observations
data.detect_transmission <- data.detect_transmission[data.detect_transmission$n.obs.shared >= PARAM$thresh.min_shared_obs, ]

###################################
#Correct p.values on transmission Z scores and generate genus-level plots
for (coupled.sites in site.combinations[1:2]) {
  if (nrow(data.taxa.transmission[[coupled.sites]]) == 0) {next}
  data.taxa.transmission[[coupled.sites]]$Z.bg_full.all.wilcox.p <- p.adjust(data.taxa.transmission[[coupled.sites]]$Z.bg_full.all.wilcox_p.raw, method="BH")
  data.taxa.transmission[[coupled.sites]]$Z.bg_full.indiv.wilcox.p <- p.adjust(data.taxa.transmission[[coupled.sites]]$Z.bg_full.indiv.wilcox_p.raw, method="BH")
  data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.p <- p.adjust(data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox_p.raw, method="BH")
  data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.indiv.wilcox.p <- p.adjust(data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.indiv.wilcox_p.raw, method="BH")
  
  #Set significance levels as factor
  data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.sig <- "FALSE"
  curr.no_na <- !is.na(data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.p)
  data.taxa.transmission[[coupled.sites]][data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.p < 0.05 & curr.no_na, "Z.bg_cohort.all.wilcox.sig"] <- "<0.05"
  data.taxa.transmission[[coupled.sites]][data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.p < 0.01 & curr.no_na, "Z.bg_cohort.all.wilcox.sig"] <- "<0.01"
  data.taxa.transmission[[coupled.sites]][data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.p < 0.001 & curr.no_na, "Z.bg_cohort.all.wilcox.sig"] <- "<0.001"
  
  #Assign logical vector "transmitter"
  data.taxa.transmission[[coupled.sites]]$transmitter <- F
  data.taxa.transmission[[coupled.sites]][data.taxa.transmission[[coupled.sites]]$Z.bg_cohort.all.wilcox.sig != "FALSE", "transmitter"] <- T
}
###################################


###################################
#Reorder taxa factors in results frames
taxa.order.sci_name <- as.character(data.taxa$Scientific_Name[data.taxa$Scientific_Name %in% data.detect_transmission$Taxon])
data.detect_transmission$Taxon <- factor(data.detect_transmission$Taxon, levels=taxa.order.sci_name)
###################################


#Re-format data (drop superfluous taxa)
plot.detect_transmission <- data.detect_transmission
plot.detect_transmission$Taxon <- droplevels(plot.detect_transmission$Taxon)
plot.taxa <- data.taxa
plot.taxa <- plot.taxa[plot.taxa$Tax_ID %in% plot.detect_transmission$Tax_ID, ]
plot.taxa$Scientific_Name <- droplevels(plot.taxa$Scientific_Name)

#Generate per-taxon Z score boxplots for bg_cohort.all
p.boxplot.transmission <- ggplot(plot.detect_transmission, aes(x=Taxon, y=S.z.bg_cohort.all)) +
  #Add boxplot
  geom_boxplot(alpha=0.8) +
  #Add lines to indicate Z score significances
  geom_hline(yintercept=0) + geom_hline(yintercept=1.645) + geom_hline(yintercept=2.326) + geom_hline(yintercept=3.090) +
  #Invert order of taxa axis back to natural
  scale_x_discrete(limits=rev(levels(plot.detect_transmission$Taxon))) +
  ylab("Transmission Score per Individual\n(Z Score vs Background)") +
  facet_grid(. ~ coupled.sites) +
  #Put taxa on y axis
  coord_flip() +
  theme_bw()

p.boxplot.transmission <- ggplot(data.taxa.transmission[["saliva_gut"]], aes(x=Scientific_Name, y=transmissibility.rank_mean)) +
  geom_point() +
  #Add lines to indicate Z score significances
  #geom_hline(yintercept=0) + geom_hline(yintercept=1.645) + geom_hline(yintercept=2.326) + geom_hline(yintercept=3.090) +
  #Invert order of taxa axis back to natural
  scale_x_discrete(limits=rev(levels(plot.detect_transmission$Taxon))) +
  ylab("Transmission Score per Individual\n(Z Score vs Background)") +
  facet_grid(. ~ coupled.sites) +
  #Put taxa on y axis
  coord_flip() +
  theme_bw()

################################################################################
################################################################################


#Save data
save(data.detect_transmission, data.taxa.transmission, file=paste0(PARAM$folder.results, "data.detect.transmission.RData"))
################################################################################
################################################################################



#q(save="no")





