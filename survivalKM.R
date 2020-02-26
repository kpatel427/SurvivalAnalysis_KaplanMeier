# script to generate survival (KM) plots using expression
library(reshape2)
library(data.table)
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)

# function call = returns KM plot
survivalKM(exprData = expression_data, metadata = associated_metadata, endpoint = 'overall survival', Risk = 'All', gene = 'ALK')

# function definition
survivalKM <- function(exprData, metadata, gene, endpoint, Risk) {

  if(endpoint == "overall survival") {
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
  
  
  # Step2: Normalizing expression data ----------------------------------
  # z-score > 0 = high expression
  # z-score < 0 = low expression
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
             ggtheme = theme_bw(),
             risk.table.fontsize = 5,
             palette = c("#E7B800", "#2E9FDF"),
             title = plotTitle) + xlab('Survival Time')

  return(plot)
} # function ends here

