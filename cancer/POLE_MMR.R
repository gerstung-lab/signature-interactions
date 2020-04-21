##########################################################################
## MCMC sampling model for interaction effects in POLE and MMR cancers  ##
## N. Volkova, EBI, November 2018                                       ##
##########################################################################

# Filtering and running criteria:
# Uterine cancers (TCGA-UCEC)
# Missense mutations in the proofreading domains of POLE or POLD1
# Bethesda MSI status
# 6 signatures, 1 factor (POLE+MMR)


library(VariantAnnotation)
library(MASS)
source('../useful_functions.R')
library(ggplot2)
library(reshape2)
source('../plotting_functions.R')
library(greta)
library(openxlsx)

cosine <- function(x,y) {
  x %*% y / sqrt(sum(x**2)) / sqrt(sum(y**2))
}
types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")

# TCGA data
PATH_TO_TCGA_MATRIX='path_to_table_with_mutation_counts_for_TCGA'
bigmat <- read.table(PATH_TO_TCGA_MATRIX,sep='\t')
bigmat <- bigmat[rowSums(bigmat) > 50 & rowSums(bigmat) < 50000,]

PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)
alternative_rownames_bigmat <- paste0('TCGA',substr(rownames(bigmat),5,nchar(row.names(bigmat))))
donor.mut.mat <- bigmat[metadata$project[match(alternative_rownames_bigmat, metadata$tumour)] == 'UCEC',]
rownames(donor.mut.mat) <- sapply(rownames(donor.mut.mat), function(x) gsub(pattern = 'H_LR', replacement = 'TCGA', x = x))
donor.mut.mat[is.na(donor.mut.mat)] <- 0
donor.mut.mat <- donor.mut.mat[rowSums(donor.mut.mat[,1:96])>0,]
ucec <- rownames(donor.mut.mat)

# clinical data from TCGA
PATH_TO_UCEC='path_to_UCEC_data_from_ICGC/'
donor <- read.delim(paste0(PATH_TO_UCEC,'donor.tsv'), sep='\t', header=T)
clin.sample <- read.delim(paste0(PATH_TO_UCEC,'sample.tsv'), sep='\t', header=T)
clin.specimen <- read.delim(paste0(PATH_TO_UCEC,'specimen.tsv'), sep='\t', header=T)
ucec <- ucec[substr(ucec,1,16) %in% clin.specimen$submitted_specimen_id]
data <- data.frame(row.names = ucec, age = donor$donor_age_at_diagnosis[match(substr(ucec,1,12), donor$submitted_donor_id)])

# POLE mutations in TCGA
# get mutations using the mutations.sh bash script
PATH_TO_POLE_MUTATIONS='path_to_table_with_POLE_mutations_per_sample'
tmp <- read.delim(PATH_TO_POLE_MUTATIONS, sep='\t', header = T)
tmp <- tmp[as.character(tmp$Effect)=='missense' & tmp$Gene=='POLE' & nchar(as.character(tmp$Protein))==7,]
tmp <- tmp[as.numeric(substr(as.character(tmp$Protein),4,6)) < 472 & 
             as.numeric(substr(as.character(tmp$Protein),4,6)) > 267,]
data$POLE <- gsub(ucec,pattern = 'TCGA',replacement = 'H_LR') %in% as.character(tmp$Sample)
tmp <- tmp[as.character(tmp$Sample) %in% gsub(ucec,pattern = 'TCGA',replacement = 'H_LR'),]
data$whichPOLE <- tmp$Protein[match(gsub(ucec,pattern = 'TCGA',replacement = 'H_LR'), as.character(tmp$Sample))]

data$POLE_1 <- data$whichPOLE == 'p.P286R'
data$POLE_2 <- data$whichPOLE == 'p.V411L'

# MSI data from TCGA Clinical Explorer
# go to http://genomeportal.stanford.edu/pan-tcga/data_download
# select UCEC
# press Go
msi <- read.table('UCEC/UCEC_2015-04-02_ClinicalParameters.txt', sep='\t', header=T)
data$MMR <- msi$MSIstatus[match(substr(ucec,1,12),msi$SampleCode)]
data$MMR <- data$MMR=='MSI-H'

ucec <- ucec[-which(is.na(data$MMR))]
data <- data[-which(is.na(data$MMR)),]
u <- donor.mut.mat[ucec,]
data <- data[,-match('whichPOLE',colnames(data))]
data[is.na(data)] <- F

# Run signature extraction with MCMC
ucec.data <- list(
  y = as.matrix(donor.mut.mat[ucec,]),
  R = 104,
  N = length(ucec),
  M = max(rowSums(donor.mut.mat[ucec,])),
  S = 5,
  X = cbind(as.numeric(data$POLE_1), as.numeric(data$POLE_2), as.numeric(data$POLE * data$MMR > 0)),
  K = 3
)

S <- variable(lower = 0, upper = 1, dim = c(ucec.data[['S']],ucec.data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = ucec.data[['R']]))
E <- variable(lower = 0, upper = ucec.data[['M']], dim = c(ucec.data[['N']],ucec.data[['S']]))

sigma = 0.5
beta = normal(mean = 0, sd = sigma, dim = c(ucec.data[['K']],ucec.data[['R']]))

size = 50

mu1 = E[,-(ucec.data[['S']]),drop=F] %*% S[-(ucec.data[['S']]),,drop=F] + 
  (E[,ucec.data[['S']],drop=F] %*% S[(ucec.data[['S']]),,drop=F]) * exp(ucec.data[['X']] %*% beta) / 
  ((exp(ucec.data[['X']] %*% beta) %*% t(S[ucec.data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))
prob = size/(size + mu1)
distribution(ucec.data[['y']]) = negative_binomial(size = matrix(1, nrow = ucec.data[['N']], ncol = ucec.data[['R']]) * size, prob = prob)
m <- model(S,E,beta)

# sampling
draws <- mcmc(m,n_samples = 500, warmup = 500,chains=4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(ucec.data[['S']]*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[1:9]])

# inspect the chains - may swap POLE and MMR
draws_all <- draws[[4]]

S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = ucec.data[['S']], ncol = ucec.data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = ucec.data[['S']], ncol = ucec.data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = ucec.data[['S']], ncol = ucec.data[['R']])
S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = ucec.data[['S']], ncol = ucec.data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:ucec.data[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = ucec.data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = ucec.data[['K']], ncol = ucec.data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = ucec.data[['K']], ncol = ucec.data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = ucec.data[['K']], ncol = ucec.data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = ucec.data[['K']], ncol = ucec.data[['R']])

plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_subindel_wb(cbind(t(S_est), S_est[ucec.data[['S']],] * exp(beta_est[1,]),S_est[ucec.data[['S']],] * exp(beta_est[2,]),S_est[ucec.data[['S']],] * exp(beta_est[3,])))
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1))
interaction_effect_plot_human(beta_est[2,], lwd = 2, CI = T, low = beta_low[2,], high = beta_high[2,], at = c(-1,0,1))
interaction_effect_plot_human(beta_est[3,], lwd = 2, CI = T, low = beta_low[3,], high = beta_high[3,], at = c(-1,0,1))

fit <- glm(E_est[data$POLE,6] ~ (data$POLE*data$MMR)[data$POLE], family = stats::poisson(link = 'identity'))

pdf('POLE_boxplot.pdf',3,4)
par(mar = c(4,4,2,2))
f <- beeswarm(log10(E_est[data$POLE,6]) ~ (data$POLE*data$MMR)[data$POLE], pch = 16, col = 'gray74', 
              bty = 'n', ylab = 'Total burden', xaxt = 'n', xlab = '',
              corral = 'wrap',
              method='hex', yaxt = 'n',
              corralWidth = 0.3)
axis(side = 1, at = c(1,2), labels = c('POLE','POLE+MMR'), col = 'white')
axis(side = 2, at = c(1.3,2,3,4,4.6), labels = c(20,100,1000,10000,40000), las = 2)
axis(side = 2, at = c(log10(c(c(3:9)*10,c(2:9)*100,c(2:9)*1000,20000,30000))), labels = rep('',25), las = 2)
boxplot(log10(E_est[data$POLE,6]) ~ (data$POLE*data$MMR)[data$POLE], frame = F, add = T,
        col=rgb(1,1,1,alpha=0.2), outline = F, xaxt = 'n', yaxt = 'n', staplewex=0, lty = 1,boxwex=.5)
axis(side=3, tcl=0.25, at=1:2, labels=NA)
mtext(side=3, paste0("P = ",signif(summary(fit)$coef[2,4],1)))
dev.off()


mu <- E_est %*% S_est

plot(as.vector(ucec.data[['y']]),
     as.vector(mu), pch = 16)
abline(a=0,b=1,col='red',lty=2)

# ids - clinical data on all samples from ICGC (donot.tsv)
AGE = ids$donor_age_at_diagnosis[match(substr(rownames(ucec.data[['y']]),1,12), ids$submitted_donor_id)]
S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = ucec.data[['S']], ncol = ucec.data[['R']])
tmp <- data.frame(S_est[6,],
                  S_est[6,] * exp(beta_est[3,]))
var_tmp <- exp(2*beta_est[3,]) * S_var[6,] + S_est[6,]**2 * exp(2*beta_est[3,]) * beta_var[3,] # variance of S_est[2,] * exp(beta_est[1,])
tmp_low <- data.frame(S_low[6,],
                      S_est[6,] * exp(beta_est[3,]) - 1.96 * sqrt(var_tmp))
tmp_up <- data.frame(S_high[6,],
                     S_est[6,] * exp(beta_est[3,]) + 1.96 * sqrt(var_tmp))
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('POLE','POLE + MMR')

plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())

no_muts_average <- mean(rowSums(ucec.data[['y']])[data$POLE>0 & ucec.data[['X']][,3]==0] / AGE[data$POLE>0 & ucec.data[['X']][,3]==0], na.rm = T)
tmp$POLE <- tmp$POLE * no_muts_average
tmp_low$POLE <- tmp_low$POLE * no_muts_average
tmp_up$POLE <- tmp_up$POLE * no_muts_average

muts_average <- mean(rowSums(ucec.data[['y']])[data$POLE>0 & ucec.data[['X']][,3]==1]/AGE[data$POLE>0 & ucec.data[['X']][,3]==1], na.rm = T)
tmp_low$`POLE + MMR` <- tmp_low$`POLE + MMR` * muts_average / sum(tmp$`POLE + MMR`)
tmp_up$`POLE + MMR` <- tmp_up$`POLE + MMR` * muts_average / sum(tmp$`POLE + MMR`)
tmp$`POLE + MMR` <- tmp$`POLE + MMR` * muts_average / sum(tmp$`POLE + MMR`)

interaction_effect_plot_human(log10(tmp$`POLE + MMR` / tmp$POLE),
                              CI = F, log = T, at = c(-1,0,1))


pdf('POLE_signatures.pdf',12,5)
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())
dev.off()

sd_pole_lfc <- sqrt(1/tmp$`POLE + MMR`**2 * var_tmp * muts_average**2 + 
                    1/tmp$POLE**2  * S_var[6,] * no_muts_average**2)

POLE_LFC <- mean(tmp$`POLE + MMR` / tmp$POLE) # 1.057215
POLE_LFC_sd <- sqrt(sum(sd_pole_lfc**2) / 104**2) # 0.06272744

pdf('POLE_LFC.pdf',12,4)
interaction_effect_plot_human(log(tmp$`POLE + MMR` / tmp$POLE), lwd = 2, CI = T,
                              low = log(tmp$`POLE + MMR` / tmp$POLE) - 1.96*sd_pole_lfc,
                              high = log(tmp$`POLE + MMR` / tmp$POLE) + 1.96*sd_pole_lfc, at = c(-1,log10(0.5),0,log10(2),log10(5),1,2),
                              labels = c('<0.1',0.5,1,2,5,10,100),
                              plot_main = 'POLE change w.r.t. MSI status')
dev.off()