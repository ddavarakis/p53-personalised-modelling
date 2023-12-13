
# Assignment 2: Modelling Project
# Personalised modelling of the p53 DNA damage response

# Q1 : Can the p53 model predict patient survival (use TCGA data)?

# References for citated genes
# https://www.mdpi.com/2072-6694/13/9/2125
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316127/
# RAF family : https://www.frontiersin.org/articles/10.3389/fcell.2020.00217/full
# MYC family: https://pubmed.ncbi.nlm.nih.gov/8960139/#:~:text=Overexpression%20and%20amplification%20of%20the,the%20progression%20of%20colorectal%20cancer.
# RAS family: https://oncologypro.esmo.org/education-library/factsheets-on-biomarkers/ras-in-colorectal-cancer

# Dimitrios – Tilemachos Davarakis
# UCD student ID: 22200956


suppressPackageStartupMessages({
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(ggsurvfit)
  library(gtsummary)
  library("glmnet")
  library("glmpath")
  library("glmnet")
  library(survival)
  library(survminer)
  library(lubridate)
  library(tidycmprsk)
  library(preprocessCore)
})


#set directory
setwd("C:\\Users\\mtheo\\Desktop\\assignment 2")
setwd("C:\\Users\\davarakis\\Desktop\\Anaconda\\PrecisionMedicine")

p53_folder = "p53 model python\\"
# Dataset Used: coadread tcga pan can atlas 2018

# ------------------------------

# READ p53_model predictions!
filename = paste0(p53_folder,"data\\count_predictions\\comp_data_new_surv.csv")
filename
all_predictions <- read.csv(filename)

# Can p53_s15 or p53_46 states predict progression free survival?

# Create survival object
os_surv_obj <- Surv(time = all_predictions$OS_MONTHS, 
                     event = all_predictions$OS_STATUS)


# ----------------------------
# Run UNIVARIATE cox regression for each p53_s15 state or p53_s46 state
# Univariate cox regression analysis identify of individual factors 
# related to patient survival

p53_cox_univariate <- function(p53_state, range) {
  vec <- paste0(p53_state, range)
  vec
  vec <- append(vec,c("OS_STATUS", "OS_MONTHS"),after=0)
  print(vec)
  data <- all_predictions[,vec]
  tbl_uvregression(
    data,
    method=coxph,
    y = Surv(time = OS_MONTHS, event = OS_STATUS),
    exponentiate = TRUE,
    #include = -ID
  )
}
# p53_s15
# DDR range log[-3,1]
p53_cox_univariate("s15_", c(1:10))
# DDR range log[1,5]
p53_cox_univariate("s15_", c(11:20))
# DDR range log[5,10]
p53_cox_univariate("s15_", c(21:30))
# p53_s46
# DDR range log[-3,1]
p53_cox_univariate("s46_", c(1:10))
# DDR range log[1,5]
p53_cox_univariate("s46_", c(11:20))
# DDR range log[5,10]
p53_cox_univariate("s46_", c(21:30))

# Statistical significant p53 states are:
p53_sign_states <- c("s46_8")


# get p-values
p53_sign_states_pvalues <- c()
for (p53_state in p53_sign_states) {
  #print(p53_state)
  fS <- os_surv_obj ~ .
  fs = reformulate(p53_state, fS[[2]])
  fs
  fit.coxph <- coxph(fs, data = all_predictions)
  sum <- summary(fit.coxph)
  p_value <- sum$coefficients[1,5]
  print(paste0("for ", p53_state, " pvalue is: ", p_value))
  p53_sign_states_pvalues <- append(p53_sign_states_pvalues, p_value)
}
p53_sign_states_pvalues


# ------------------------------------
# Analysis with citated genes
### Does the genes, reported by the literature as those genes 
# that might regulate p53 activate in 
# colorectal cancer, also influence the PFS ?  

citated_genes <- c("MYC",	"KRAS", "NRAS",	"ARF1",	
                   "MDM4", "ARL4C", "ARL5A", "SAR1B", 
                   "AGAP2", "ASAP1", "SMAP1", 
                   # best 5 scored DNA_PK genes
                   # refer to https://www.genecards.org/Search/Keyword?queryString=dna-pk
                   "PRKDC", "XRCC5", "XRCC6", "LIG4", "XRCC4")  

fS <- os_surv_obj ~ . 
citated_genes
fs = reformulate(citated_genes, fS[[2]])
fs
fit.coxph <- coxph(fs, data = all_predictions)
summary(fit.coxph)
fit.coxph %>% 
  tbl_regression(exp = TRUE)
ggforest(fit.coxph, data = all_predictions)
#!!! good !!!

# Multivariant Cox Regression for each p53_significant_state + "SAR1B", "XRCC4", "PRKDC"
for (p53_state in p53_sign_states) {
  print(p53_state)
  fS <- os_surv_obj ~ .
  genes <- c("SAR1B", "XRCC4", "PRKDC")
  formula <- append(genes, p53_state, after=0)
  fs = reformulate(formula, fS[[2]])
  fs
  fit.coxph <- coxph(fs, data = all_predictions)
  summary(fit.coxph)
  a<-fit.coxph %>% 
    tbl_regression(exp = TRUE)
  print(a)
}
p53_sign_states

#save citated genes
library(xlsx)
df <- data.frame(citated_genes)
write.xlsx(df,"os_significant_genes.xlsx",sheetName = "os_citated_genes")

# ---------------------------------
# Multivariant Cox regression analysis 
# for: significant_p53_state + "SAR1B", "XRCC4", "PRKDC" + stage2 + stage3 + stage4

stages <- unique(all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE)
stages

# calculate cancer stages!
all_predictions$stage2 <- ifelse(all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE II" |
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIA" | 
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIB" |
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIC" , "Yes", "No")
all_predictions$stage3 <- ifelse(all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE III" |
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIA" | 
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIB" |
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIC" , "Yes", "No")
all_predictions$stage4 <- ifelse(all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IV" |
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IVA" | 
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IVB" |
                                   all_predictions$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IVC" , "Yes", "No")


for (p53_state in p53_sign_states) {
  print(p53_state)
  fS <- os_surv_obj ~ .
  others <- c("SAR1B", "XRCC4", "PRKDC", "stage2", "stage3", "stage4")
  formula <- append(others, p53_state, after=0)
  fs = reformulate(formula, fS[[2]])
  fs
  fit.coxph <- coxph(fs, data = all_predictions)
  summary(fit.coxph)
  #fit.coxph %>% 
  #  tbl_regression(exp = TRUE)
  print(ggforest(fit.coxph, data = all_predictions))
}


# -------------------------------------
# Fit Regularized Cox regression (by using glmnet R package) to find out the 
# genes that statistically influence the PFS survival. 
# Prepare the data for glmnet function

# get only averaged TCGA gene data
my_x <- all_predictions[,63:20580]
# convert RNA data from dataframe to matrix
my_x <- as.matrix(my_x)
# get the  time
my_time <- all_predictions$OS_MONTHS
# get the status
my_status <- all_predictions$OS_STATUS
# create a matrix with time and status
my_y <- as.matrix(data.frame(time=my_time,status=my_status))

# create a dataframe that contain time, status and RNA data
my_data <- data.frame(time=my_time,status=my_status,my_x)

# Create a survival object (for overall survival)
my_z <- Surv(my_time, my_status)

# Fit the glmnet - Regularized Cox regression
# It takes some time !
my_fit <- glmnet(my_x, my_z, family = "cox")
print(my_fit)
plot(my_fit)
title("Regularized Cox Regression on PFS Survival", line = 2.5, cex.main = 1)

# Analysis
my_cfs = coef(my_fit)
summary(my_cfs)
# extract the coefficients at a certain value of λ 
# λ = 0.077290 returns the 21 most significant genes
# λ = 0.067230 returns the 32 most significant genes  (given to Cillian)
my_cfs = coef(my_fit, s = 0.067230) # gets 32 genes

# get the relevant genes
my_meaning_coefs = rownames(my_cfs)[my_cfs[,1]!= 0]
my_meaning_coefs


# dynamically built of formula
fS <- my_z ~ . 
fs = reformulate(my_meaning_coefs, fS[[2]])
fs

# use the built formula to fit the Cox regression
fit.coxph <- coxph(fs, data = all_predictions)

ggforest(fit.coxph, data = my_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)

# find genes that are statistically significant
a <- summary(fit.coxph)
genes <- rownames(a$coefficients)
length(genes)
stat_genes <- c()
for (i in 1:length(genes)) {
  if (a$coefficients[i,5] <= 0.05) {
    print(genes[i])
    stat_genes <- append(stat_genes,genes[i])
  }
}
stat_genes
length(stat_genes)

#save significant genes
df <- data.frame(stat_genes)
write.xlsx(df,"os_significant_genes.xlsx",sheetName = "os_significant_genes",append = TRUE,row.names = FALSE)


for (p53_state in p53_sign_states) {
  print(p53_state)
  fS <- os_surv_obj ~ .
  formula <- append(stat_genes, p53_state, after=0)
  fs = reformulate(formula, fS[[2]])
  fs
  fit.coxph <- coxph(fs, data = all_predictions)
  summary(fit.coxph)
  a<-fit.coxph %>% 
    tbl_regression(exp = TRUE)
  print(a)
  print(ggforest(fit.coxph, data = all_predictions))
}


# ----------------------------------
# Kaplan Meier Plots for each statistically significant p53 state
# Run a cycle by changing a threshold for p53state-separation into low and
# high groups, from min to max, and find a least p-value, 
# meaning the best statistical significance. 

km_data <- all_predictions
k = 1
for (p53state_name in p53_sign_states) {
  #print(p53state_name)
  p53threshold <- p53_sign_states_pvalues[k]
  #print(p53threshold)
  
  # starts
  text <- paste0("Computing & Ploting KM for: ", p53state_name)
  print(text)
  # start
  p53state <- km_data[p53state_name]
  p53state_strata <- paste0(p53state_name, "_strata")
  p53state_sorted <- sort(p53state[,1])
  #length(p53state_sorted)
  for (i in 2:(length(p53state_sorted)-1)) {
    cur_cutoff = p53state_sorted[i]
    if (cur_cutoff>0) {
      km_data[p53state_strata] <- ifelse(p53state < cur_cutoff, "LOW", "HIGH")
      
      form <- paste0("os_surv_obj ~ ",p53state_strata)
      
      fit <- survfit(as.formula(form), data=km_data)
      w <- surv_pvalue(fit)
      if (w$pval <= p53threshold){
        #print(paste0("Found pvalue < threshold: ", w$pval))
        minimum_pvalue<- w$pval
        icutoff <- i
        break
      }
    }
  }
  #print(minimum_pvalue)
  #print(icutoff)
  best_cut_off <- p53state_sorted[icutoff]
  #print(best_cut_off)
  km_data[p53state_strata] <- ifelse(p53state < best_cut_off, "LOW", "HIGH")
  
  form <- paste0("os_surv_obj ~ ",p53state_strata)
  form
  
  fit <- survfit(as.formula(form), data=km_data)
  fit
  print(ggsurvplot(fit,
                   data = km_data,
                   pval = T,
                   risk.table = T,
                   xlab = "Months", ylab = "PFS survival"))
  
  # ends
  
  k <- k + 1
}





