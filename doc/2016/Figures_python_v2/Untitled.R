#library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
library(psych)


setwd('/Users/kanaan/Google Drive/TS-EUROTRAIN/Papers/2016_QSM_paper/Figures_python_v2')
source("/Users/kanaan/SCR/R/ggbiplot2.R")

datadir = '/Users/kanaan/Google Drive/TS-EUROTRAIN/RESULTS_QSMv3/AUG5/phenotypic'
dfp = read.table(file.path(datadir,'df_patients_qc.csv'), sep = ",", header=TRUE, row.names = 1)

# Remove outlier HHQP (row#12) 
#dfp <- dfp[-c(3),]

data  = dfp[,c('CLN_YGTSS_Motoric_Score',
               'CLN_YGTSS_Vocal_Score', 
               #'CLN_YGTSS_Total_Score_incl_Impairment',
               #'CLN_YGTSS_Total_Tic_Score',
               'CLN_RVTRS',
               #'CLN_qol_score',
               'CLN_puts',
               'CLN_DSM4_ADHD_Score',
               #'CLN_CAARS_Score_ADHS_Symptoms_Total_Cat_G_T_Score',
               #'CLN_YBOCS_Totalscore_Items_1to10', 
               'CLN_OCIR_total_score',
               'CLN_BAI',
               'CLN_BDI12',
               'CLN_MADRS'
)] 

colnames(data) = c('YGTSS-MOTOR', 'YGTSS-Vocal', 'RVTRS', 'PUTS', 'DSM4-ADHD',#'CAARS', 'YBOCS', 
                   'OCI-R', 'BAI', 'BDI-II', 'MADRS')



# Pariwise Scatter Plots
#pairs(data)


PCA2 = principal(data, nfactors=9)
biplot(PCA2)
