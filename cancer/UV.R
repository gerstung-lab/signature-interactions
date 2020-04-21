########################################################
## MCMC sampling model for interaction effects in UV  ##
## N. Volkova, EBI, November 2018                     ##
########################################################
library(VariantAnnotation)
library(MASS)
source('~/Desktop/Git/phd/useful_functions.R')
library(ggplot2)
library(rstan)
library(reshape2)
source('~/Desktop/Git/phd/plot_sigs.R')
library(greta)
library(openxlsx)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")

# Similarity
cosine <- function(x,y) {
  x %*% y / sqrt(sum(x**2)) / sqrt(sum(y**2))
}
# COSMIC signatures
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
types <- as.character(cancer_signatures$Trinucleotide) # trinucleotide classes
types.full <- as.character(cancer_signatures$Somatic.Mutation.Type) # substitution types
row.names(cancer_signatures) <- types.full
cancer_signatures = as.matrix(cancer_signatures[,4:33])
# Data

PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)
skcm <- as.character(metadata$tumour[metadata$project=='SKCM'])
# Mutations
PATH_TO_SKCM='path_to_matrix_with_mutation_counts_for_96_subs_8_indel_types_and_DNV'
uv <- read.table(PATH_TO_SKCM,sep='\t',header=T)
uv <- uv[rowSums(uv[,1:96])>0, ]
# indels =  "D.1.5"    "D.5.50"   "D.50.400" "DI.small" "DI.large" "I.1.5"    "I.5.50"   "I.50.400"
types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")
colnames(uv) <- c(types.full,indels,'DNV')
# Select the ones looking like UV

skcm <- rownames(uv)[sapply(rownames(uv), function(x) cosine(as.numeric(uv[x,1:96]),cancer_signatures[,7]))>0.8]
uv <- uv[skcm,]

#############################
# Gather info on mutations (or get it from Supplementary Table)

mutations_in_samples <- read.xlsx("Supplementary_Tables/Supplement/Supplementary Table 4.xlsx", sheet = 1, startRow = 2)
samples <- lapply(ner.genes, function(x) unique(mutations_in_samples$Sample[grep(x, mutations_in_samples$Mutated_genes)]))

# Germline mutations (if relevant)
SNPs <- list()
PATH_TO_GERMLINE_MUTS_IN_NER='path_to_table_with_germline_mutations_in_NER'
snps_per_sample <- read.table(PATH_TO_GERMLINE_MUTS_IN_NER, sep = '\t', header = T)

thr = 0.2
# Check defects per gene
X <- data.frame(sapply(ner.genes),
            function(x) {
              return(skcm %in% unlist(samples[[x]]) | skcm %in% metadata$tumour[metadata[,x] < thr & !is.na(metadata[,x])])
}))

# Data for sampling
skcm.data <- list(
  N = nrow(uv),
  R = 105,
  S = 1,
  K = ncol(X),
  X = model.matrix( ~ ., data = X)[,-1],
  y = as.matrix(uv),
  M = max(rowSums(uv))
)

# Sampling model
S <- variable(lower = 0, upper = 1, dim = c(skcm.data[['S']],skcm.data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = skcm.data[['R']]))
E <- variable(lower = 0, upper = skcm.data[['M']], dim = c(skcm.data[['N']],skcm.data[['S']]))
sigma = 0.1
beta = normal(mean = 0, sd = sigma, dim = c(skcm.data[['K']],skcm.data[['R']]))
size = 50

mu1 = (E %*% S) * exp(skcm.data[['X']] %*% beta) / 
    ((exp(skcm.data[['X']] %*% beta) %*% t(S)) %*% matrix(1,nrow = 1, ncol=104))
prob = size/(size + mu1)
distribution(skcm.data[['y']]) = negative_binomial(size = size * matrix(1,nrow=skcm.data[['N']],ncol=skcm.data[['R']]),prob = prob)
                                                    
m <- model(S,E,beta)

# Sampling
draws <- greta::mcmc(m, n_samples = 200, warmup = 500, chains = 4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(skcm.data[['S']]*c(1:104),12)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(skcm.data[['K']]*c(1:104),12)]])

draws_all <- do.call('rbind',draws)

S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',c(1:skcm.data[['S']]))
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = skcm.data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = skcm.data[['K']], ncol = skcm.data[['R']])

pdf('~/UV_real_predicted.pdf',7,7)
mu_uv <- E_est %*% S_est * exp(skcm.data[['X']] %*% beta_est) / 
  ((exp(skcm.data[['X']] %*% beta_est) %*% t(S_est)) %*% matrix(1,nrow = 1, ncol=105))
plot(as.vector(skcm.data[['y']]),
     as.vector(mu_uv), pch = 16, xlab = 'Observed', ylab = 'Predicted', bty = 'n')
abline(a=0,b=1,col='red',lty=2)
boxplot(log10(E_est[,1]) ~ skcm.data[['X']][,1],bty = 'n', xaxt = 'n', frame = F, yaxt = 'n')
axis(side = 1, at = c(1,2), labels = c('Normal','NER def.'), tick = T, col = 'white', col.ticks = 'black')
axis(side = 2, at = c(2,3,4), labels = c(100,1000,10000), las = 2, tick = T, col = 'white', col.ticks = 'black')
dev.off()

pdf('UV_LFC.pdf',12,4)
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T,
                              low = beta_low[1,],
                              high = beta_high[1,], at = c(-1,0,1),
                              labels = c('<0.1',1,10),
                              plot_main = 'UV change w.r.t. NER mutations')
dev.off()

S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
tmp <- data.frame(S_est[1,],
                  S_est[1,] * exp(beta_est[1,]))
var_tmp <- exp(2*beta_est[1,]) * S_var[1,] + S_est[1,]**2 * exp(2*beta_est[1,]) * beta_var[1,] # variance of S_est[1,] * exp(beta_est[2,])
tmp_low <- data.frame(S_low[1,],
                      S_est[1,] * exp(beta_est[1,]) - 1.96 * sqrt(var_tmp))
tmp_up <- data.frame(S_high[1,],
                     S_est[1,] * exp(beta_est[1,]) + 1.96 * sqrt(var_tmp))
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('UV','UV + muts')
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())


AGE = ids$donor_age_at_diagnosis[match(substr(rownames(skcm.data[['y']]),1,12), ids$submitted_donor_id)]
no_muts_average <- mean(rowSums(skcm.data[['y']])[skcm.data[['X']][,1]==0] / AGE[skcm.data[['X']][,1]==0], na.rm = T)
tmp$UV <- tmp$UV * no_muts_average
tmp_low$UV <- tmp_low$UV * no_muts_average
tmp_up$UV <- tmp_up$UV * no_muts_average

muts_average <- mean(rowSums(skcm.data[['y']])[skcm.data[['X']][,1]==1] / AGE[skcm.data[['X']][,1]==1], na.rm = T)
tmp_low$`UV + muts` <- tmp_low$`UV + muts` * muts_average / sum(tmp$`UV + muts`)
tmp_up$`UV + muts` <- tmp_up$`UV + muts` * muts_average / sum(tmp$`UV + muts`)
tmp$`UV + muts` <- tmp$`UV + muts` * muts_average / sum(tmp$`UV + muts`)

interaction_effect_plot_human(log10(tmp$`UV + muts` / tmp$UV),
                              CI = F, log = T, at = c(-1,0,1))

pdf('UV_signatures.pdf',12,5)
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())
dev.off()

sd_uv_lfc <- sqrt(1/tmp$`UV + muts`**2 * var_tmp * muts_average**2 + 
                    1/tmp$UV**2  * S_var[1,] * no_muts_average**2)

UV_LFC <- mean(tmp$`UV + muts` / tmp$UV) # 1.814657
UV_LFC_sd <- sqrt(sum(sd_uv_lfc**2) / 104**2) # 0.04520428

pdf('UV_LFC.pdf',12,4)
interaction_effect_plot_human(log(tmp$`UV + muts` / tmp$UV), lwd = 2, CI = T,
                              low = log(tmp$`UV + muts` / tmp$UV) - 1.96*sd_uv_lfc,
                              high = log(tmp$`UV + muts` / tmp$UV) + 1.96*sd_uv_lfc, at = c(-1,0,log10(2),1),
                              labels = c('<0.1',1,2,10),
                              plot_main = 'UV change w.r.t. TC-NER mutations')
dev.off()

# Checking the effect on mutational burden
counts <- rowSums(skcm.data[['y']])
counts <- skcm.data[['y']][,105]
AGE = ids$donor_age_at_diagnosis[match(substr(rownames(skcm.data[['y']]),1,12), ids$submitted_donor_id)]
mdata <- data.frame(y = counts / AGE, model.matrix( ~ ., data = X)[,-1])
model <- glm(y ~ NERTRUE + MMRTRUE + TLSTRUE + BERTRUE, family = gaussian(link = 'log'), data = mdata)
model <- glm(y ~ NERTRUE + MMRTRUE + TLSTRUE + BERTRUE, family = gaussian(), data = mdata)
summary(model)

counts <- skcm.data[['y']][,105] / rowSums(skcm.data[['y']])
mdata <- data.frame(y = counts, model.matrix( ~ ., data = X)[,-1])
model <- glm(y ~ NERTRUE + MMRTRUE + TLSTRUE + BERTRUE, family = gaussian(), data = mdata)
summary(model)

counts <- rowSums(skcm.data[['y']])
mdata <- data.frame(y = counts / AGE, log(metadata[match(skcm,metadata$tumour),c(ner.genes,mmr[-3],tls)]))
model <- glm(y ~ ., family = gaussian(), data = mdata)
summary(model)

library(VGAM)
x <- as.matrix(mdata[,-1])
model <- gam(y ~ s(XPC, bs = 're') + 
               s(ERCC1, bs = 're') +
             s(ERCC2, bs = 're') +
               s(ERCC3, bs = 're') + s(ERCC4, bs = 're') + s(ERCC5, bs = 're'), data = mdata, family = gaussian, smooth = x)
summary(model)

############################################
## Analysis of XP and non-XP cSCC genomes ##
############################################

PATH_TO_XP_MUTATIONS <- 'path_to_table_with_mutation_counts_per_sample'
bigmat <- as.matrix(read.table(PATH_TO_XP_MUTATIONS,sep=' '))
colnames(bigmat) <- c('XP1','XP2','XP3','XP4','XP5',paste0('N',c(1:8)))
metadata <- matrix(0,nrow = ncol(bigmat), ncol = 1)
metadata[1:5,1] <- 1
colnames(metadata) <- 'XPC'

# Data for sampling
skcm.data <- list(
  N = ncol(bigmat),
  R = 96,
  S = 1,
  K = 1,
  X = metadata,
  y = t(bigmat),
  M = max(colSums(bigmat))
)

# Sampling model
S <- variable(lower = 0, upper = 1, dim = c(skcm.data[['S']],skcm.data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = skcm.data[['R']]))
E <- variable(lower = 0, upper = skcm.data[['M']], dim = c(skcm.data[['N']],skcm.data[['S']]))
sigma = 0.1
beta = normal(mean = 0, sd = sigma, dim = c(skcm.data[['K']],skcm.data[['R']]))
size = 50
mu1 = (E %*% S) * exp(skcm.data[['X']] %*% beta) / 
  ((exp(skcm.data[['X']] %*% beta) %*% t(S)) %*% matrix(1,nrow = 1, ncol=skcm.data[['R']]))
prob = size/(size + mu1)
distribution(skcm.data[['y']]) = negative_binomial(size = size * matrix(1,nrow=skcm.data[['N']],ncol=skcm.data[['R']]),prob = prob)
skcm.model <- model(S,E,beta)

# Sampling
draws <- greta::mcmc(skcm.model, n_samples = 500, warmup = 500)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(skcm.data[['S']]*c(1:104),12)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(skcm.data[['K']]*c(1:104),12)]])

draws_all <- do.call('rbind',draws)

S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',c(1:skcm.data[['S']]))
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = skcm.data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = skcm.data[['K']], ncol = skcm.data[['R']])

plot_all_changes(plot_sig_wb, S_est, beta_est)

pdf('UV_XP_boxplot.pdf',4,6)
library(beeswarm)
AGE = c(10,10,10,10,10,76,87,84,61,83,85,58,63)
boxplot(log10(E_est[,1] / AGE) ~ skcm.data[['X']][,1], frame = F, yaxt = 'n', xaxt = 'n', xlab = '', ylab = 'Number of mutations / year', ylim = c(2,5))
beeswarm(log10(E_est[,1] / AGE) ~ metadata[,1], pch = 16, col = 'gray', add = T)
axis(side = 1, lty = 0, at = c(1,2), labels = c('wt','XP'))
axis(side = 2, las = 2, at = c(2,3,4,5), labels = c("10^2","10^3","10^4","10^5"), lty = 0)
axis(side = 2, las = 2, at = c(log10(c(1:9)*100),log10(c(1:9)*1000),log10(c(1:9)*10000), 5), labels = rep("",28))
dev.off()

interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1))

pdf('~/UV_real_predicted.pdf',7,7)
mu_uv <- E_est %*% S_est * exp(skcm.data[['X']] %*% beta_est) / 
  ((exp(skcm.data[['X']] %*% beta_est) %*% t(S_est)) %*% matrix(1,nrow = 1, ncol=skcm.data[['R']]))
plot(as.vector(skcm.data[['y']]),
     as.vector(mu_uv), pch = 16, xlab = 'Observed', ylab = 'Predicted', bty = 'n')
abline(a=0,b=1,col='red',lty=2)
boxplot(log10(E_est[,1]) ~ skcm.data[['X']][,1],bty = 'n', xaxt = 'n', frame = F, yaxt = 'n')
axis(side = 1, at = c(1,2), labels = c('Normal','NER def.'), tick = T, col = 'white', col.ticks = 'black')
axis(side = 2, at = c(2,3,4), labels = c(100,1000,10000), las = 2, tick = T, col = 'white', col.ticks = 'black')
dev.off()

pdf('UV_XP_LFC.pdf',12,4)
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T,
                              low = beta_low[1,],
                              high = beta_high[1,], at = c(-1,0,1),
                              labels = c('<0.1',1,10),
                              plot_main = 'UV change w.r.t. XPC deficiency')
dev.off()

#S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
tmp <- data.frame(S_est[1,],
                  S_est[1,] * exp(beta_est[1,]))
var_tmp <- exp(2*beta_est[1,]) * S_var[1,] + S_est[1,]**2 * exp(2*beta_est[1,]) * beta_var[1,] # variance of S_est[1,] * exp(beta_est[2,])
tmp_low <- data.frame(S_low[1,],
                      S_est[1,] * exp(beta_est[1,]) - 1.96 * sqrt(var_tmp))
tmp_up <- data.frame(S_high[1,],
                     S_est[1,] * exp(beta_est[1,]) + 1.96 * sqrt(var_tmp))
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('UV','UV + muts')



tmz_sim <- 1 - cosine(tmp[,1],tmp[,2]) # 0.567953
prod = sum(tmp[,1] * tmp[,2])
len1 <- sqrt(sum(tmp[,1]**2))
len2 <- sqrt(sum(tmp[,2]**2))
var <- sum(((tmp[,1] * len2**2 - tmp[,2] * prod) / len2**3 / len1)**2 *tmp[,2])
tmz_sim_sd <- var # 0.03880034



pdf('UV_XP_signatures.pdf',12,5)
plot_sig_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())
dev.off()
plot_sig_wb(cbind(tmp,sbs_cancer_signatures[,'SBS7a']+sbs_cancer_signatures[,'SBS7b']),norm = F) + theme(panel.grid = element_blank())

AGE = c(2,8,10,10,10,76,87,84,61,83,85,58,63) # , Durinck et al. 2011
no_muts_average <- mean(rowSums(skcm.data[['y']])[skcm.data[['X']][,1]==0] / AGE[skcm.data[['X']][,1]==0], na.rm = T)
tmp$UV <- tmp$UV * no_muts_average
tmp_low$UV <- tmp_low$UV * no_muts_average
tmp_up$UV <- tmp_up$UV * no_muts_average

muts_average <- mean(rowSums(skcm.data[['y']])[skcm.data[['X']][,1]==1] / AGE[skcm.data[['X']][,1]==1], na.rm = T)
tmp_low$`UV + muts` <- tmp_low$`UV + muts` * muts_average / sum(tmp$`UV + muts`)
tmp_up$`UV + muts` <- tmp_up$`UV + muts` * muts_average / sum(tmp$`UV + muts`)
tmp$`UV + muts` <- tmp$`UV + muts` * muts_average / sum(tmp$`UV + muts`)

interaction_effect_plot_human(log10(tmp$`UV + muts` / tmp$UV),
                              CI = F, log = T, at = c(-1,0,1))

pdf('UV_signatures_adjusted.pdf',12,5)
plot_sig_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F,
            diff_scale = T) + theme(panel.grid = element_blank())
dev.off()

sd_uv_lfc <- sqrt(1/tmp$`UV + muts`**2 * var_tmp * muts_average**2 + 
                    1/tmp$UV**2  * S_var[1,] * no_muts_average**2)

UV_LFC <- mean(tmp$`UV + muts` / tmp$UV) # 1.814657
UV_LFC_sd <- sqrt(sum(sd_uv_lfc**2) / 104**2) # 0.04520428

pdf('UV_LFC_adjusted.pdf',12,3)
interaction_effect_plot_human(log(tmp$`UV + muts` / tmp$UV), lwd = 2, CI = T,
                              low = log(tmp$`UV + muts` / tmp$UV) - 1.96*sd_uv_lfc,
                              high = log(tmp$`UV + muts` / tmp$UV) + 1.96*sd_uv_lfc,
                              at = c(0,1,2),
                              labels = c(1,10,100),
                              plot_main = 'UV change w.r.t. TC-NER mutations')
dev.off()
