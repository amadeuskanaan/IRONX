#library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
library(psych)
#devtools::install_github("kassambara/factoextra")
#library("factoextra")

setwd('/Users/kanaan/Google Drive/TS-EUROTRAIN/Papers/2016_QSM_paper/Figures_python_v2')
datadir = '/Users/kanaan/Google Drive/TS-EUROTRAIN/RESULTS_QSMv3/SEPT10/AHBA'

geneset = 'STR_IRON_T2'
  
df = read.table(file.path(datadir,paste(geneset, '.csv', sep = "")), sep = ",", header=TRUE, row.names = 1)
drops <- c('coords_native', 'donor_names', 'struct_id', 'struct_name', 'top_struct', 'Mean', 'Median', 
           'corrected_mni_x', 'corrected_mni_y', 'corrected_mni_z', 'SVD1g', 'SVD2g', 'SVD3g', 'SVD1p', 'SVD2p', 'SVD3p')
data = df[ , !(names(df) %in% drops)]

# Perform PCA
#PCA <- princomp(data, scores = TRUE, cor = TRUE) 
#summary(PCA)
#PCA = principal(data, nfactors=5)
PCA = principal(data, nfactors=3)
PCA

# Scree plot
#plot(PCA, type='l')

# biplot
#biplot(PCA )
#ggbiplot(PCA, obs.scale = 0, var.scale = 0, ellipse = FALSE, circle = TRUE)

# Print Loadings
PCA$loadings

# Print Scores (projections)
#PCA$scores

df_scores   = PCA$scores
df_loadings = PCA$loadings
#colnames(df_scores) = c('PC1', 'PC2', 'PC3') 
#colnames(df_scores) = c('PC1', 'PC2', 'PC3') 
dfmerge = merge(df, df_scores, by=0, all=TRUE)
row.names(dfmerge) <- dfmerge$Row.names
dfmerge = dfmerge[ , -which(names(dfmerge) %in% c("Row.names"))]
write.csv(dfmerge, file = file.path(datadir, paste('AHBA_', geneset, '_PCA.csv', sep = "")))
write.csv(df_loadings, file = file.path(datadir, paste('loadings_', geneset, '.csv', sep = "")))

