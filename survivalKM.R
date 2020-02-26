# script to generate survival (KM) plots using expression
#setwd("/Volumes/target_nbl_ngs/KP/RShiny")
library(reshape2)
library(data.table)
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)

# fit <- survfit(Surv(time, status) ~ sex, data = lung) 
# In the code, time is overall/event-free survival time, status is 0 (alive) or 1 (dead), 
# instead of sex you will use gene expression binary values for the gene of interest: 0 (low expression) or 1 (high expression)
# ggsurvplot(fit, data = lung)


load('./data/GSE3960_data.RData')
load('./data/GSE3960_mData.RData')

load('./data/TARGET_NBL_FPKM_PST_data.RData')
load('./data/TARGET_NBL_FPKM_PST_mData.RData')


exprData = TARGET_NBL_FPKM_PST_data
metadata = TARGET_NBL_FPKM_PST_mData
gene = 'PHOX2B'
Risk = 'All'
endpoint = 'os'

survivalKM(exprData = TARGET_NBL_FPKM_PST_data, metadata = TARGET_NBL_FPKM_PST_mData, endpoint = 'os', Risk = 'All', gene = 'NME1')

survivalKM <- function(exprData, metadata, gene, endpoint, Risk) {

  if(endpoint == "os") {
    time <- 'nti_surv_overall'
    status <- 'nti_event_overall_num'
  } else {
    time <- 'nti_surv_progrfree'
    status <- 'nti_event_progrfree_num'
  }

  # Step1: Prepare data ----------------------------------
  # condition for specific RISK group
  if (Risk != 'All') {
    metadata <- metadata[metadata$RISK == Risk,]
    exprData <- exprData[, which(colnames(exprData) %in% rownames(metadata))]
  } else {
    metadata <- metadata
    exprData <- exprData
  }
  
  # subset metadata
  metadata <- metadata[,c("RISK",time,status)]
  
  # subset expr data
  exprData <- exprData[gene,]
  
  
  # Step2: Normalizing expression data
  # z-score > 1 = high expression
  # z-score < -1 = low expression
  exprData <- data.frame("gene.Z" = apply(X = exprData, MARGIN = 1, FUN = scale))
  

  # adding expr data to metadata
  metadata <- cbind(metadata, exprData)
  metadata <- metadata %>%
      gather(key = 'gene', value = 'zScore', -c("RISK",time,status))

  metadata$gene <- gsub('gene.Z.','',metadata$gene)
  metadata$gene <- as.factor(metadata$gene)
  # assigning exprStatus
  # 1 = high
  # 0 = low
  # 2 = mid expression
  # ignore mid expressions
  metadata$exprStatus <- ifelse(metadata$zScore > 0, 1,
                                ifelse(metadata$zScore < 0, 0, 2))


  # remove the ones with mid expressions
  metadata <- metadata[metadata$exprStatus != 2,]

  # Step3: Fitting survival curves ----------------------------------
  group = "exprStatus"

  fit <- survfit(Surv(get(time), get(status)) ~ exprStatus, data = metadata)


  # Step4: Comparing survival times between groups ----------------------------------
  # conduct between-group significance tests using a log-rank test
  # to be conducted when more than 1 gene is selected
  if(length(gene) > 1){
    diff <- survdiff(formula = Surv(get(time), get(status)) ~ gene,
             data = metadata)
  
    pval <- pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
  
    adjpval <- p.adjust(pval, method = "BH", n = length(gene))
  
    plotTitle <- paste0(paste(gene, collapse = " | "), " P-val(Adj) :", format(pval, scientific=T, digits=3), "(", format(adjpval, scientific=T, digits=3), ")")
  } else {
   plotTitle <- paste0(gene)
 }


  # Step5: Generate plot ----------------------------------
  # change strata names
  low <- paste0("Low : n = ", fit$n[1])
  high <- paste0("High : n = ", fit$n[2])
  names(fit$strata) <- c(low,high)

  plot <- ggsurvplot(fit,
             data = metadata,
             pval = TRUE,
             pval.size = 5,
             conf.int = TRUE,
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             surv.median.line = "hv", # Specify median survival
             #ggtheme = theme_Publication_scatter(base_size = 14), # Change ggplot2 theme
             # legend = "right",
             # font.x = 14, font.y = 14, font.main = 14, font.legend = 14, font.tickslab = 14,
             risk.table.fontsize = 5,
             palette = c("#E7B800", "#2E9FDF"),
             title = plotTitle) + xlab('Survival Time')

  return(plot)
} # function ends here

