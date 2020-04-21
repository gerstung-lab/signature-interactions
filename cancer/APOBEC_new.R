# TCGA APOBEC+REV1/UNG1 analysis

library(VariantAnnotation)
library(MASS)
source('../useful_functions.R')
library(ggplot2)
library(greta)
library(reshape2)
source('../plotting_functions.R')
types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")

# Get the new signature set from here: https://www.synapse.org/#!Synapse:syn11967914
new_cancer_signatures <- read.csv('sigProfiler_exome_SBS_signatures.csv')

# Read in the matrix with mutation counts and with information on samples
PATH_TO_TCGA_MATRIX='path_to_table_with_mutation_counts_for_TCGA'
bigmat <- read.table(PATH_TO_TCGA_MATRIX,sep='\t')
colnames(bigmat) <- c(types.full, indels)
PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)

# restrictions
bigmat <- bigmat[rowSums(bigmat[,1:96])>50 & rowSums(bigmat[,1:96]) < 10000,]
colnames(bigmat) <- c(types.full, indels)
# sample name adjustment (same samples may start with TCGA or with another code)
alternative_rownames_bigmat <- paste0('TCGA',substr(rownames(bigmat),5,nchar(row.names(bigmat))))
# select all samples similar to APOBEC profile
all_apobecs <- alternative_rownames_bigmat[which(sapply(rownames(bigmat), function(x) cosine(as.numeric(bigmat[x,1:96]),
       new_cancer_signatures[,'SBS2'] + new_cancer_signatures[,'SBS13']))>0.8)]

mutations_in_samples <- read.xlsx("Supplementary_Tables/Supplement/Supplementary Table 4.xlsx", sheet = 1, startRow = 2)
genes = unique(sapply(mutations_in_samples$Mutated_genes, function(x) unlist(strsplit(x, split = ','))[1]))
samples <- lapply(genes, function(x) unique(mutations_in_samples$Sample[grep(x, mutations_in_samples$Mutated_genes)]))

# Check the fraction of C>A and C>G mutations versus C>T mutations in mutated/non-mutated samples (heterozygous+homozygous)
REVUNG <- as.numeric(substr(all_apobecs,1,12) %in% c(samples[['REV1_hom']],samples[['UNG_hom']],
                                               samples[['UNG_het']],samples[['REV1_het']]))
data <- rbind(bigmat[match(all_apobecs[REVUNG==0],alternative_rownames_bigmat),],
              bigmat[match(all_apobecs[REVUNG>0],alternative_rownames_bigmat),])
data$`T[C>R]N` <- as.matrix(data[-1]) %*% grepl("^T.C>(A|G)", colnames(data)[-1])
data$`T[C>A]N` <- as.matrix(data[-1]) %*% grepl("^T.C>A", colnames(data)[-1])
data$`T[C>G]N` <- as.matrix(data[-1]) %*% grepl("^T.C>G", colnames(data)[-1])
data$`T[C>T]N` <- as.matrix(data[-1]) %*% grepl("^T.C>T", colnames(data)[-1])
data$REVUNG <- factor(c(rep(0,sum(REVUNG==0)),rep(1,sum(REVUNG==1))))
levels(data$REVUNG) <- c("wt","mut")
fit <- glm(cbind(data$`T[C>G]N`,data$`T[C>T]N`) ~ data$REVUNG,  family=stats::binomial)
summary(fit)

library(beeswarm)
pdf('APOBEC_boxplot.pdf',3,4)
f <- beeswarm(data$`T[C>G]N` / (data$`T[C>G]N` + data$`T[C>T]N`) ~ data$REVUNG, pch = 16, col = 'gray86', 
              bty = 'n', ylab = 'Fraction of T[C>G]N', xaxt = 'n', xlab = 'REV1/UNG mutations', cex=.6,
              corral = 'wrap',
              method = 'hex',
              corralWidth = 0.3, las = 2)
axis(side = 1, at = c(1,2), labels = c('WT (792)','Mut (48)'), col = 'white')
boxplot(data$`T[C>G]N` / (data$`T[C>G]N` + data$`T[C>T]N`) ~ data$REVUNG, frame = F, add = T,
        col=rgb(1,1,1,alpha=0.2), outline = F, xaxt = 'n', yaxt = 'n', staplewex=0, lty = 1,boxwex=.5)
axis(side=3, tcl=0.25, at=1:2, labels=NA)
mtext(side=3, paste0("P = ",signif(summary(fit)$coef[2,4],1)))
dev.off()

# Prepare the data for the effect extractino model
apobec_data <- list()
apobec_data[['y']] <- bigmat[match(all_apobecs,alternative_rownames_bigmat),]
apobec_data[['X']] <- matrix(as.numeric(REVUNG),ncol=1)
apobec_data[['N']] <- nrow(apobec_data[['y']])
apobec_data[['R']] <- ncol(apobec_data[['y']])
apobec_data[['K']] <- 1
apobec_data[['M']] <- max(rowSums(apobec_data[['y']]))
apobec_data[['S']] <- 1

S <- variable(lower = 0, upper = 1, dim = c(apobec_data[['S']],apobec_data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = apobec_data[['R']]))
E <- variable(lower = 0, upper = apobec_data[['M']], dim = c(apobec_data[['N']],apobec_data[['S']]))

sigma = 0.5
beta = normal(mean = 0, sd = sigma, dim = c(apobec_data[['K']],apobec_data[['R']]))

size = 50

mu1 = (E %*% S) * exp(apobec_data[['X']] %*% beta) /
  ((exp(apobec_data[['X']] %*% beta) %*% t(S[apobec_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(apobec_data[['y']]) = negative_binomial(size = matrix(1, nrow = apobec_data[['N']], ncol = apobec_data[['R']]) * size,
                                                  prob = prob)
m <- model(S,E,beta)

# sampling
draws <- mcmc(model = m, warmup=500, n_samples=500, chains=4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(apobec_data[['S']]*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(1:(apobec_data[['K']]*104),9)]])

draws_all <- draws[[1]]
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = apobec_data[['S']], ncol = apobec_data[['R']])
S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = apobec_data[['S']], ncol = apobec_data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = apobec_data[['S']], ncol = apobec_data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = apobec_data[['S']], ncol = apobec_data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:apobec_data[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = apobec_data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = apobec_data[['K']], ncol = apobec_data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = apobec_data[['K']], ncol = apobec_data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = apobec_data[['K']], ncol = apobec_data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = apobec_data[['K']], ncol = apobec_data[['R']])

plot_subindel_wb(cbind(t(S_est),S_est[1,] * exp(beta_est[1,])))
interaction_effect_plot_human(c(beta_est[1,],0), lwd = 2, CI = T,
                              low = beta_low[1,],
                              high = beta_high[1,], at = c(-1,0,1),
                              plot_main = 'APOBEC change w.r.t. REV1/UNG mutations')

mu_apobec <- (E_est %*% S_est) * exp(apobec_data[['X']] %*% beta_est) / 
  ((exp(apobec_data[['X']] %*% beta_est) %*% t(S_est[apobec_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

pdf('~/APOBEC_quality.pdf',7,7)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(as.numeric(as.matrix(apobec_data[['y']])),
     as.vector(mu_apobec), pch = 16,
     main = 'APOBEC samples',
     bty = 'n')
abline(a=0,b=1,col='red',lty=2)
dev.off()

pdf('APOBEC_LFC.pdf',12,4)
interaction_effect_plot_human(c(beta_est[1,],0), lwd = 2, CI = T,
                              low = beta_low[1,],
                              high = beta_high[1,], at = c(-1,0,1),
                              labels = c('<0.1',1,10),
                              plot_main = 'APOBEC change w.r.t. REV1/UNG mutations')
dev.off()

S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = apobec_data[['S']], ncol = apobec_data[['R']])
tmp <- data.frame(S_est[1,],
                  S_est[1,] * exp(beta_est[1,]))
var_tmp <- exp(2*beta_est[1,]) * S_var[1,] + S_est[1,]**2 * exp(2*beta_est[1,]) * beta_var[1,] # variance of S_est[1,] * exp(beta_est[1,])
tmp_low <- data.frame(S_low[1,],
                      S_est[1,] * exp(beta_est[1,]) - 1.96 * sqrt(var_tmp))
tmp_up <- data.frame(S_high[1,],
                     S_est[1,] * exp(beta_est[1,]) + 1.96 * sqrt(var_tmp))
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('APOBEC','APOBEC + muts')
tmp_low$`APOBEC + muts` <- tmp_low$`APOBEC + muts` / sum(tmp$`APOBEC + muts`)
tmp_up$`APOBEC + muts` <- tmp_up$`APOBEC + muts` / sum(tmp$`APOBEC + muts`)
tmp$`APOBEC + muts` <- tmp$`APOBEC + muts` / sum(tmp$`APOBEC + muts`)

pdf('APOBEC_signatures_220619.pdf',12,5)
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())
dev.off()

sd_ap_lfc <- sqrt(1/tmp$`APOBEC + muts`**2 * var_tmp+ 
                    1/tmp$APOBEC**2  * S_var[1,])

pdf('APOBEC_LFC_220619.pdf',12,4)
interaction_effect_plot_human(c(log(tmp$`APOBEC + muts` / tmp$APOBEC),0), lwd = 2, CI = T,
                              low = c(log(tmp$`APOBEC + muts` / tmp$APOBEC) - 1.96*sd_ap_lfc,0),
                              high =c(log(tmp$`APOBEC + muts` / tmp$APOBEC) + 1.96*sd_ap_lfc,0), at = c(-1,log10(0.5),0,log10(2),1),
                              labels = c('<0.1',0.5,1,2,10),
                              plot_main = 'APOBEC change w.r.t. REV1/UNG mutations')
dev.off()

# Same per year
metadata <- read.table('TCGA_caveman_patients_to_files.dat', sep = '\t', header = T)
ids <- read.table('~/Downloads/ICGCtoTCGA.tsv', sep = '\t', header = T)
AGE <- ids$donor_age_at_diagnosis[match(metadata$patient_id[match(rownames(apobec_data[['y']]),
                                                                  metadata$tumour)], ids$submitted_donor_id)]
no_muts_average <- median(rowSums(apobec_data[['y']])[apobec_data[['X']]==0] / AGE[apobec_data[['X']]==0], na.rm = T)
tmp$APOBEC <- tmp$APOBEC * no_muts_average
tmp_low$APOBEC <- tmp_low$APOBEC * no_muts_average
tmp_up$APOBEC <- tmp_up$APOBEC * no_muts_average

muts_average <- median(rowSums(apobec_data[['y']])[apobec_data[['X']]==1]/AGE[apobec_data[['X']]==1], na.rm = T)
tmp_low$`APOBEC + muts` <- tmp_low$`APOBEC + muts` * muts_average / sum(tmp$`APOBEC + muts`)
tmp_up$`APOBEC + muts` <- tmp_up$`APOBEC + muts` * muts_average / sum(tmp$`APOBEC + muts`)
tmp$`APOBEC + muts` <- tmp$`APOBEC + muts` * muts_average / sum(tmp$`APOBEC + muts`)

interaction_effect_plot_human(c(log10(tmp$`APOBEC + muts` / tmp$APOBEC),0),
                              CI = F, log = T, at = c(-1,0,1))

plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = 0.2, norm = F) + theme(panel.grid = element_blank())

sd_ap_lfc <- sqrt(1/tmp$`APOBEC + muts`**2 * var_tmp * muts_average**2 + 
                    1/tmp$APOBEC**2  * S_var[1,] * no_muts_average**2)

APOBEC_LFC <- mean(tmp$`APOBEC + muts` / tmp$APOBEC) # 1.05
APOBEC_LFC_sd <- sqrt(sum(sd_ap_lfc**2) / 104**2) # 0.016
# cosine 0.99

pdf('APOBEC_LFC.pdf',12,4)
interaction_effect_plot_human(c(log(tmp$`APOBEC + muts` / tmp$APOBEC),0), lwd = 2, CI = T,
                              low = c(log(tmp$`APOBEC + muts` / tmp$APOBEC) - 1.96*sd_ap_lfc,0),
                              high =c(log(tmp$`APOBEC + muts` / tmp$APOBEC) + 1.96*sd_ap_lfc,0), at = c(-1,log10(0.5),0,log10(2),1),
                              labels = c('<0.1',0.5,1,2,10),
                              plot_main = 'APOBEC change w.r.t. REV1/UNG mutations')
dev.off()

