# script to scan a dataset to get genes with significant difference in survival, based on expression
# also outputs which group has poor survival
# plots KM plot if survival probability of any one group is 0
# does not plot for those genes with only 1 sample in either of the group with 0 survival probability
# setwd("~/KP/significantSurvival")

library(reshape2)
library(data.table)
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)


# Read Data ------------------
load('~/data/expr_data.RData')
load('~/data/expr_mData.RData')

# empty df to store pvals and survival probabilities for all genes
result <- data.frame()


# loop to call function -------------------------
for(i in row.names(expr_data)){
  print(i)
  
  # function call = returns KM plot
  df <- survivalKM(exprData = expr_data, 
                   metadata = expr_mData, 
                   endpoint = 'overall survival', 
                   Risk = 'All', 
                   gene = as.character(i))
  
  result <- rbind(result,df)
}

# result formatting -----------------------
# omit NAs
result <- na.omit(result) 
# OS/EFS =243 genes did not have both expr groups

# filter for significant diffs
result <- result[result$pval < 0.05,] 
# OS = only 3309 genes have significant difference in survival
# EFS = only 1580 genes have significant difference in survival


# to annotate worseOutcomes
result <- result %>%
  group_by(gene) %>%
  mutate('worseOutcome'=min(survProb))

result <- result[result$survProb == result$worseOutcome,]

result$strata <- gsub('exprStatus=0','Low', result$strata)
result$strata <- gsub('exprStatus=1','High', result$strata)

result <- result[,c(-5,-6)]

# write to a table
write.table(result[,c(1,4)], file = 'allGenes_Target_signfOS.txt', col.names = T, row.names = F, sep = '\t', quote = F)


# function definition ------------------------------------------------
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
  
  # if no two expression groups present, set NAs for that gene
  if(length(unique(metadata$exprStatus)) != 2){
    
    merged <- data.frame(gene = gene, pval = NA, adjpval = NA, strata = NA, survProb = NA)
    
  } else {
    
    # Step3: Fitting survival curves ----------------------------------
    group = "exprStatus"
    
    fit <- survfit(Surv(get(time), get(status)) ~ exprStatus, data = metadata)
    # Median survival is the time corresponding to a survival probability of 0.5 (50% survival)
    
    # to get survival probabilty at a specific time
    summ <- summary(fit, times = 4500, extend = TRUE)
    #df <- data.frame('strata'=summ$strata, 'survProb' = summ$surv)

    # Step4: Comparing survival times between groups ----------------------------------
    # conduct between-group significance tests using a log-rank test
    # to be conducted between high-low expressing groups
    
    diff <- survdiff(formula = Surv(get(time), get(status)) ~ exprStatus, data = metadata)
    pval <- pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
    adjpval <- p.adjust(pval, method = "fdr",  n = (fit$n[1] + fit$n[2]))
    
    # plot for those genes whose survival probability is zero
    if(0 %in% summ$surv){
        
        plotTitle <- paste0(paste(gene, collapse = " | "), " P-val(Adj) :", format(pval, scientific=T, digits=3), "(", format(adjpval, scientific=T, digits=3), ")")
    
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
    
        pdf(paste0('plots/OS/',gene,"_OS_survplot.pdf"), width = 12, height = 9)
        print(plot, newpage = FALSE)
        dev.off()
    } # if block ends here

    merged <- data.frame(gene = gene, pval = pval, adjpval = adjpval, strata = summ$strata, survProb = summ$surv)

  } # else block ends here
  
  return(merged)
  
  
} # function ends here
