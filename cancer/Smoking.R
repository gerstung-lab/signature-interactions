##################################################################
## MCMC sampling model for interaction effects in lung cancers  ##
## N. Volkova, EBI, November 2018                               ##
##################################################################

library(VariantAnnotation)
library(MASS)
source('../useful_functions.R')
library(ggplot2)
library(reshape2)
source('../plotting_functions.R')
library(greta)
library(openxlsx)

NER.core <- c('XPC','XPA','CUL5','ERCC1','ERCC2','ERCC4','ERCC5','ERCC6','CUL3')
TLS.core <- c('POLN','POLQ','REV3L','REV1','SHPRH','POLK','POLH','MAD2L2','POLI')

types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")

# Get the new signature set from here: https://www.synapse.org/#!Synapse:syn11967914
new_cancer_signatures <- read.csv('sigProfiler_exome_SBS_signatures.csv')

mutations_in_samples <- read.xlsx("Supplementary_Tables/Supplement/Supplementary Table 4.xlsx", sheet = 1, startRow = 2)
samples <- lapply(c(NER.core, TLS.core), function(x) unique(mutations_in_samples$Sample[grep(x, mutations_in_samples$Mutated_genes)]))

PATH_TO_TCGA_MATRIX='path_to_table_with_mutation_counts_for_TCGA'
bigmat <- read.table(PATH_TO_TCGA_MATRIX,sep='\t')
rownames(bigmat) <- paste0('TCGA',substr(rownames(bigmat),5,nchar(rownames(bigmat))))
PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)

lung <- metadata$tumour[metadata$project %in% c('LUAD','LUSC')]
lung <- lung[which(sapply(match(lung,alternative_rownames_bigmat), function(x) cosine(as.numeric(bigmat[x,1:96]), new_cancer_signatures[,'SBS4'])) > 0.8)]
lung.mut.mat <- bigmat[match(lung,rownames(bigmat)),]

# Pull together
X <- data.frame(NER = as.numeric(substr(lung,1,12) %in% unlist(samples[NER.core])),
                TLS = as.numeric(substr(lung,1,12) %in% unlist(samples[TLS.core])))

lung_data <- list(
  N = nrow(lung.mut.mat),
  R = ncol(lung.mut.mat),
  K = ncol(X)-1,
  X = as.matrix(X)[,1],
  y = as.matrix(lung.mut.mat),
  M = max(rowSums(lung.mut.mat)),
  S = 2
)

# Sampling
library(greta)
S <- variable(lower = 0, upper = 1, dim = c(lung_data[['S']],lung_data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = lung_data[['R']]))
E <- variable(lower = 0, upper = lung_data[['M']], dim = c(lung_data[['N']],lung_data[['S']]))

sigma = 0.5
beta = normal(mean = 0, sd = sigma, dim = c(lung_data[['K']],lung_data[['R']]))

size = 50

mu1 = E[,-lung_data[['S']],drop=F] %*% S[-lung_data[['S']],,drop=F] + 
  (E[,lung_data[['S']],drop=F] %*% S[lung_data[['S']],,drop=F]) * exp(lung_data[['X']] %*% beta) / 
  ((exp(lung_data[['X']] %*% beta) %*% t(S[lung_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(lung_data[['y']]) = negative_binomial(size = matrix(1, nrow = lung_data[['N']], ncol = lung_data[['R']]) * size, prob = prob)
# defining the model
m <- model(S,E,beta)

# sampling
draws <- mcmc(m, n_samples = 500, warmup = 500, chains = 4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(lung_data[['S']]*c(1:104),12)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(lung_data[['K']]*c(1:104),12)]])

# chains where effect was assigned to smoking signature 
draws_all <- draws[[3]]
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = lung_data[['S']], ncol = lung_data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = lung_data[['S']], ncol = lung_data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = lung_data[['S']], ncol = lung_data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',c(1:lung_data[['S']]))
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = lung_data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = lung_data[['K']], ncol = lung_data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = lung_data[['K']], ncol = lung_data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = lung_data[['K']], ncol = lung_data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = lung_data[['K']], ncol = lung_data[['R']])

plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1))


AGE <- ids$donor_age_at_diagnosis[match(metadata$patient_id[match(rownames(lung_data[['y']]),
                                                                  metadata$tumour)], ids$submitted_donor_id)]

S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = lung_data[['S']], ncol = lung_data[['R']])
tmp <- data.frame(S_est[2,],
                  S_est[2,] * exp(beta_est[1,]))
var_tmp <- exp(2*beta_est[1,]) * S_var[2,] + S_est[2,]**2 * exp(2*beta_est[1,]) * beta_var[1,] # variance of S_est[1,] * exp(beta_est[1,])
tmp_low <- data.frame(S_low[2,],
                      S_est[2,] * exp(beta_est[1,]) - 1.96 * sqrt(var_tmp))
tmp_up <- data.frame(S_high[2,],
                     S_est[2,] * exp(beta_est[1,]) + 1.96 * sqrt(var_tmp))
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('TOBACCO','TOBACCO + muts')
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, norm = F, ymax = 0.1) + theme(panel.grid = element_blank())

no_muts_average <- mean(rowSums(lung_data[['y']])[lung_data[['X']][,1]==0] / AGE[lung_data[['X']][,1]==0], na.rm = T)
tmp$TOBACCO <- tmp$TOBACCO * no_muts_average
tmp_low$TOBACCO <- tmp_low$TOBACCO * no_muts_average
tmp_up$TOBACCO <- tmp_up$TOBACCO * no_muts_average

muts_average <- mean(rowSums(lung_data[['y']])[lung_data[['X']][,1]==1]/AGE[lung_data[['X']][,1]==1], na.rm = T)
tmp_low$`TOBACCO + muts` <- tmp_low$`TOBACCO + muts` * muts_average / sum(tmp$`TOBACCO + muts`)
tmp_up$`TOBACCO + muts` <- tmp_up$`TOBACCO + muts` * muts_average / sum(tmp$`TOBACCO + muts`)
tmp$`TOBACCO + muts` <- tmp$`TOBACCO + muts` * muts_average / sum(tmp$`TOBACCO + muts`)

interaction_effect_plot_human(log10(tmp$`TOBACCO + muts` / tmp$TOBACCO),
                              CI = F, log = T, at = c(-1,0,1))

pdf('smoking_signatures.pdf',12,5)
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = 1, norm = F) + theme(panel.grid = element_blank())
dev.off()

sd_sm_lfc <- sqrt(1/tmp$`TOBACCO + muts`**2 * var_tmp * muts_average**2 + 
                    1/tmp$TOBACCO**2  * S_var[1,] * no_muts_average**2)

SMOKING_LFC <- mean(tmp$`TOBACCO + muts` / tmp$TOBACCO) # 0.8
SMOKING_LFC_sd <- sqrt(sum(sd_sm_lfc**2) / 104**2) # 0.05

pdf('smoking_LFC.pdf',12,4)
interaction_effect_plot_human(log10(tmp$`TOBACCO + muts` / tmp$TOBACCO), lwd = 2, CI = T,log=T,
                              low = log10(tmp$`TOBACCO + muts` / tmp$TOBACCO) - 2*sd_sm_lfc,
                              high = log10(tmp$`TOBACCO + muts` / tmp$TOBACCO) + 2*sd_sm_lfc, at = c(-1,0,log10(2),1),
                              labels = c('<0.1',1,2,10),
                              plot_main = 'TOBACCO change w.r.t. NER mutations')
dev.off()



