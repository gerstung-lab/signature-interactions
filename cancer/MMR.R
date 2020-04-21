#############################################################################
## MCMC sampling model for interaction effects in MMR-def cancers          ##
## N. Volkova, EBI, November 2018                                          ##
#############################################################################

library(VariantAnnotation)
library(MASS)
source('../useful_functions.R')
library(ggplot2)
library(reshape2)
source('../plotting_functions.R')
library(greta)
library(openxlsx)


types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")

# Get the new signature set from here: https://www.synapse.org/#!Synapse:syn11967914
new_cancer_signatures <- read.csv('sigProfiler_exome_SBS_signatures.csv')

PATH_TO_TCGA_MATRIX='path_to_table_with_mutation_counts_for_TCGA'
bigmat <- read.table(PATH_TO_TCGA_MATRIX,sep='\t')

PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)
rownames(bigmat) <- paste0('TCGA', substr(rownames(bigmat), 5,nchar(rownames(bigmat))))

# MSI status for TCGA from TCGA Clinical Explorer
files <- list.files('~/Downloads/MSI_TCGA/')
cancer_types <- as.character(sapply(files, function(x) unlist(strsplit(x, split='[_]'))[1]))
status <- lapply(files, function(x) read.table(paste0('~/Downloads/MSI_TCGA/',x), sep='\t', header=T))
msi_status = data.frame(sample = do.call('c', lapply(status, function(x) as.character(x$SampleCode))),
                        cancer = rep(cancer_types, times = sapply(status, function(x) nrow(x))),
                        MSI = do.call('c', lapply(status, function(x) as.character(x$MSIstatus))))

# Labels
mmr.genes <- c('MLH1','PMS2','MSH2','MSH3','MSH6')
#mmr.genes <- MMR.core[-6]

mutations_in_samples <- read.xlsx("Supplementary_Tables/Supplement/Supplementary Table 4.xlsx", sheet = 1, startRow = 2)
samples <- lapply(ner.genes, function(x) unique(mutations_in_samples$Sample[grep(x, mutations_in_samples$Mutated_genes)]))

new.mmr <- unique(c(as.character(unlist(samples[mmr.genes])),
             as.character(msi_status$sample[msi_status$MSI == 'MSI-H'])))
new.mmr <- alternative_rownames_bigmat[substr(alternative_rownames_bigmat,1,12) %in% new.mmr |
                                         alternative_rownames_bigmat %in% new.mmr]

new.mmr <- new.mmr[rowSums(bigmat[match(new.mmr,alternative_rownames_bigmat),1:96]) < 20000 & rowSums(bigmat[match(new.mmr,alternative_rownames_bigmat),1:96]) > 100]

# Remove what is clearly not MMR
new.mmr <- new.mmr[which(sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), new_cancer_signatures[,'SBS2']+
                                                              new_cancer_signatures[,'SBS13'])) < 0.7 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), new_cancer_signatures[,'SBS7a']+
                                                                new_cancer_signatures[,'SBS7b'])) < 0.7 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), new_cancer_signatures[,'SBS10a']+
                                                                new_cancer_signatures[,'SBS10b'])) < 0.7 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), new_cancer_signatures[,'SBS4'])) < 0.7)] 
table(metadata$project[match(new.mmr, metadata$tumour)])


mmr_data <- list()
mmr_data[['y']] <- bigmat[match(new.mmr,alternative_rownames_bigmat),]
mmr_data[['X']] <- model.matrix( ~ .,
    data = data.frame(v1 = as.character(metadata$project[match(new.mmr, metadata$tumour)])))[,-1]
colSums(mmr_data[['X']])
mmr_data[['X']] <- data.frame(BREAST = as.numeric(mmr_data[['X']][,'v1BRCA']>0),
                              CERVIX = as.numeric(mmr_data[['X']][,'v1CESC']>0),
                              HNSC = as.numeric(mmr_data[['X']][,'v1HNSC']>0),
                              STOMACH = as.numeric(mmr_data[['X']][,'v1STAD']>0),
                              COLORECT = as.numeric((mmr_data[['X']][,'v1COAD'] + mmr_data[['X']][,'v1READ']) > 0),
                              LIVER = as.numeric(mmr_data[['X']][,'v1LIHC'] > 0),
                              LUNG = as.numeric((mmr_data[['X']][,'v1LUSC']+mmr_data[['X']][,'v1LUAD']) > 0),
                              PROSTATE = as.numeric(mmr_data[['X']][,'v1PRAD']>0),
                              UTERUS = as.numeric((mmr_data[['X']][,'v1UCEC'] + mmr_data[['X']][,'v1UCS'])>0))
ind <- which(rowSums(mmr_data[['X']])>0)
mmr_data[['X']] <- mmr_data[['X']][ind,]
mmr_data[['y']] <- mmr_data[['y']][ind,]

mmr_data[['X']] <- mmr_data[['X']][,-match('COLORECT',colnames(mmr_data[['X']]))] # COLORECT - default

mmr_data[['N']] <- nrow(mmr_data[['y']]) 
mmr_data[['R']] <- ncol(mmr_data[['y']]) 
mmr_data[['K']] <- ncol(mmr_data[['X']])
mmr_data[['S']] <- 3
mmr_data[['M']] <- max(rowSums(mmr_data[['y']]))

# Greta model
S <- variable(lower = 0, upper = 1, dim = c(mmr_data[['S']],mmr_data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = mmr_data[['R']]))
E <- variable(lower = 0, upper = mmr_data[['M']], dim = c(mmr_data[['N']],mmr_data[['S']]))

sigma = 0.5
beta = normal(mean = 0, sd = sigma, dim = c(mmr_data[['K']],mmr_data[['R']]))

size = 50

mu1 = E[,-(mmr_data[['S']]),drop=F] %*% S[-(mmr_data[['S']]),,drop=F] + 
  (E[,mmr_data[['S']],drop=F] %*% S[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta) / 
((exp(mmr_data[['X']] %*% beta) %*% t(S[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(mmr_data[['y']]) = negative_binomial(size = matrix(1, nrow = mmr_data[['N']], ncol = mmr_data[['R']]) * size, prob = prob)

#m <- model(S,E)
m <- model(S,E,beta)

# sampling
draws <- mcmc(m,n_samples = 1000, warmup = 2000, chains = 4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(mmr_data[['S']]*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(1:(mmr_data[['K']]*104),9)]])

draws_all <- draws[[4]] # 1,3
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:mmr_data[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = mmr_data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data[['K']], ncol = mmr_data[['R']])

rownames(beta_est) <- colnames(mmr_data[['X']])
plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_all_changes(plot_subindel_wb, S_est, beta_est)

mu <- E_est[,-(mmr_data[['S']]),drop=F] %*% S_est[-(mmr_data[['S']]),,drop=F] + 
  (E_est[,mmr_data[['S']],drop=F] %*% S_est[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta_est)

mu <- E_est %*% S_est * exp(mmr_data[['X']] %*% beta_est) / ((exp(mmr_data[['X']] %*% beta_est) %*% t(S_est[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

pdf('~/MMR_observed_predicted.pdf',5,5)
plot(rowSums(mmr_data[['y']]),rowSums(mu), pch = 16,
     xlab = 'Observed', ylab= 'Predicted')
abline(a=0,b=1,col='red',lty=2)
dev.off()

rownames(beta_est) <- colnames(mmr_data[['X']])
altered_signatures <- data.frame(S_est[mmr_data[['S']],],
                                 S_est[mmr_data[['S']],] * exp(beta_est[1,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[2,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[3,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[4,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[5,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[6,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[7,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[8,]))
colnames(altered_signatures) <- c('COLORECT',rownames(beta_est))
var_tmp <- t(rbind(S_var[mmr_data[['S']],],exp(2*beta_est) * S_var[mmr_data[['S']],] + S_est[mmr_data[['S']],]**2 * exp(2*beta_est) * beta_var)) # variance of S_est[1,] * exp(beta_est[1,])

altered_low <- altered_signatures - 1.96*sqrt(var_tmp)
altered_high <- altered_signatures + 1.96*sqrt(var_tmp)

plot_subindel_wb(altered_signatures, CI = T, low = altered_low, high = altered_high)

pdf('MMR_signature.pdf',12,7)
plot_subindel_wb(altered_signatures[,c(1,5,7,9)], CI = T, low = altered_low[,c(1,5,7,9)], high = altered_high[,c(1,5,7,9)]) + theme(panel.grid = element_blank())
dev.off()

altered_sums <- colSums(altered_signatures)
altered_sums_var <- colSums(var_tmp)

new_var_tmp <- sqrt(1/altered_sums**2 * var_tmp + altered_signatures**2 / altered_sums**4 * altered_sums_var)
new_altered_signatures <- apply(altered_signatures,2,function(x) x/sum(x))

plot_subindel_wb(new_altered_signatures, CI = T, low = new_altered_signatures - 1.96*new_var_tmp, high = new_altered_signatures + 1.96 * new_var_tmp)



bigmat <- bigmat[rowSums(bigmat) < 20000 & rowSums(bigmat) > 100,]
alternative_rownames_bigmat <- paste0('TCGA',substr(rownames(bigmat),5,nchar(rownames(bigmat))))

# mean and SE of the mean
mean_c_to_t <- list()
mean_c_to_t_sd <- list()
mmr_mean_c_to_t <- list()
mmr_mean_c_to_t_sd <- list()
mean_del <- list()
mean_del_sd <- list()
mmr_mean_del <- list()
mmr_mean_del_sd <- list()
mean_ins <- list()
mean_ins_sd <- list()
mmr_mean_ins <- list()
mmr_mean_ins_sd <- list()
ids <- read.table('~/Downloads/ICGCtoTCGA.tsv', sep = '\t', header = T)
age <- ids$donor_age_at_diagnosis[match(substr(alternative_rownames_bigmat,1,12),ids$submitted_donor_id)]
ctype <- sapply(alternative_rownames_bigmat, function(x) metadata$project[match(x,metadata$tumour)])
for (proj in c('COAD','BRCA','CESC','HNSC','STAD','LIHC','LUSC','PRAD','UCEC')) {
  
  rest_of_samples <- setdiff(alternative_rownames_bigmat[ctype == proj],new.mmr)
  

  mean_c_to_t[[proj]] <- mean(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mean_c_to_t_sd[[proj]] <- sd(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) / 
    sum(!is.na(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))

  mean_del[[proj]] <- mean(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mean_del_sd[[proj]] <- sd(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) / 
    sum(!is.na(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)]))

  
  mean_ins[[proj]] <- mean((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mean_ins_sd[[proj]] <- sd((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) / 
    sum(!is.na((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))

  rest_of_samples <- intersect(alternative_rownames_bigmat[ctype == proj],new.mmr)
  mmr_mean_c_to_t[[proj]] <- mean(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mmr_mean_c_to_t_sd[[proj]] <- sd(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) / 
    sum(!is.na(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))

  mmr_mean_del[[proj]] <- mean(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mmr_mean_del_sd[[proj]] <- sd(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) / 
    sum(!is.na(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)]))

  mmr_mean_ins[[proj]] <- mean((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mmr_mean_ins_sd[[proj]] <- sd((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) / 
    sum(!is.na((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))
}

pdf('MMR_C_T_scatterplot.pdf',4,4)
plot(y = unlist(mmr_mean_c_to_t), 
     x = unlist(mean_c_to_t),pch = 16,col='red',
     xlab = 'Non-MMR', ylab = 'MMR', bty = 'n', main = 'CpG > TpG',
     ylim = c(0,15), xlim = c(0,15))
for (j in 1:length(mean_c_to_t)) {
  lines(y = c(mmr_mean_c_to_t[[j]] - 1.96*mmr_mean_c_to_t_sd[[j]],mmr_mean_c_to_t[[j]] + 1.96*mmr_mean_c_to_t_sd[[j]]),
        x = c(mean_c_to_t[[j]],mean_c_to_t[[j]]), col = 'red', lwd = 1)
  lines(x = c(mean_c_to_t[[j]] - 1.96*mean_c_to_t_sd[[j]],mean_c_to_t[[j]] + 1.96*mean_c_to_t_sd[[j]]),
        y = c(mmr_mean_c_to_t[[j]],mmr_mean_c_to_t[[j]]), col = 'red', lwd = 1)
}
lm(unlist(mmr_mean_c_to_t) ~ unlist(mean_c_to_t))
abline(a = 0, b = 3, lty = 2)
#text(x = 6, y = 12, expression(y == -0.76 + 3.2 * x))
text(x = unlist(mean_c_to_t), y = unlist(mmr_mean_c_to_t), labels = names(unlist(mean_c_to_t)))
dev.off()

pdf('MMR_DEL_scatterplot.pdf',4,4)
plot(y = unlist(mmr_mean_del), 
     x = unlist(mean_del), pch = 16,col='orange',
     xlab = 'Non-MMR', ylab = 'MMR', bty = 'n', main = '1bp deletions',
     xlim = c(0,4), ylim = c(0,4))
for (j in 1:length(mean_del)) {
  lines(y = c(mmr_mean_del[[j]] - 1.96*mmr_mean_del_sd[[j]],mmr_mean_del[[j]] + 1.96*mmr_mean_del_sd[[j]]),
        x = c(mean_del[[j]],mean_del[[j]]), col = 'orange', lwd = 2)
  lines(x = c(mean_del[[j]] - 1.96*mean_del_sd[[j]],mean_del[[j]] + 1.96*mean_del_sd[[j]]),
        y = c(mmr_mean_del[[j]],mmr_mean_del[[j]]), col = 'orange', lwd = 1.5)
}
lm(unlist(mmr_mean_del) ~ 0 + unlist(mean_del))
abline(a = 0, b = 8.6, lty = 2)
text(x = unlist(mean_del), y = unlist(mmr_mean_del), labels = names(unlist(mean_del)))
#text(x = 1.2, y = 4, expression(y == -0.5 + 7.3 * x))
dev.off()

pdf('MMR_INS_scatterplot.pdf',4,4)
plot(y = unlist(mmr_mean_ins), 
     x = unlist(mean_ins), pch = 16,col='purple',
     xlab = 'Non-MMR', ylab = 'MMR', bty = 'n', main = '1 bp insertions',
     xlim = c(0,1), ylim = c(0,1))
for (j in 1:length(mean_ins)) {
  lines(y = c(mmr_mean_ins[[j]] - 1.96*mmr_mean_ins_sd[[j]],mmr_mean_ins[[j]] + 1.96*mmr_mean_ins_sd[[j]]),
        x = c(mean_ins[[j]],mean_ins[[j]]), col = 'purple', lwd = 2)
  lines(x = c(mean_ins[[j]] - 1.96*mean_ins_sd[[j]],mean_ins[[j]] + 1.96*mean_ins_sd[[j]]),
        y = c(mmr_mean_ins[[j]],mmr_mean_ins[[j]]), col = 'purple', lwd = 1.5)
}
lm(unlist(mmr_mean_ins) ~ 0 + unlist(mean_ins))
abline(a = 0, b = 5, lty = 2)
text(x = unlist(mean_ins), y = unlist(mmr_mean_ins), labels = names(unlist(mean_ins)))
dev.off()

