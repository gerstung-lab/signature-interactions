# Full NB model in greta


# Preparations and data upload

library(ggplot2)
library(openxlsx)
library(reshape2)
library(greta)
library(bayesplot)
library(MASS)

source('useful_functions.R')
source('plotting_functions.R')


data <- openxlsx::read.xlsx(xlsxFile = "Supplementary_tables/Supplement/Supplementary Table 1.xlsx", 
                            sheet = 2, cols = 1:8, start = 2)
data$Sample <- as.character(data$Sample)
data$Genotype <- as.character(data$Genotype)

CD2Mutant <- sapply(1:nrow(data), function(i) {
  if (data$Type[i] == 'mut.acc.') return(paste0(data$Genotype[i],':',data$Generation[i]))
  return(paste0(data$Genotype[i],':',data$Mutagen[i],':',data$Drug.concentration[i]))
})
names(CD2Mutant) <- data$Sample
CD2Mutant <- CD2Mutant[sort(names(CD2Mutant))]
CD2Mutant[1:5]

new.spectrum <- openxlsx::read.xlsx(xlsxFile = "Supplementary_tables/Supplement/Supplementary Table 1.xlsx", 
                          sheet = 3, rowNames = T, start = 2)
CD2Mutant <- CD2Mutant[rownames(new.spectrum)] # remove reference worms
data <- data[match(names(CD2Mutant),data$Sample),]
head(new.spectrum[,1:10])
Y <- new.spectrum

#########################################################################################################

# Design matrix


m=119;p=54;r=12;s=196;n=nrow(Y)

generation.function <- function(N) {
  if (is.na(N)) return(N)
  if (N==0) return(0)
  if (N==1) return(1)
  alpha = sum(sapply(1:N, function(i) 1/(2^(i-1)))) +
    0.25*sum(sapply(1:(N-1), function(i) sum(sapply(1:i, function(j) 1/(2^(j-1))))))
  return(alpha)
}
g = t(t(sapply(data$Generation,generation.function)))
g[g[,1] == 0,1] <- 1 # set 1 because otherwise 0-generation worms will be excluded from the model, and they may matter

# raw dose
doses <- data$Drug.concentration
doses[is.na(doses)] <- 0 # for MA experiments, replace NA dose with 0

# Mutagens
mutagens <- data$Mutagen
mutagens[is.na(mutagens)] <- 'no'
mutagens[mutagens=='MMS/EMS/DMS'] <- 'no' # N2 zero-exposure samples
# Interactions
interactions <- sapply(which(!is.na(data$Mutagen)), function(j) {
  m <- as.character(data$Mutagen[j])
  return(paste(data$Genotype[j],m,sep=":"))
})
names(interactions) <- data$Sample[which(!is.na(data$Mutagen))]
# Put all together
X <- data.frame(sapply(unique(data$Genotype), function(x) as.numeric(data$Genotype == x)),
                sapply(unique(mutagens)[-1], function(x) as.numeric(mutagens == x)),
                sapply(unique(interactions), function(x) sapply(rownames(Y), function(z) 
                  ifelse(is.na(data$Mutagen[match(z,data$Sample)]), 0, as.numeric(interactions[z] == x))))) 
rownames(X) <- data$Sample


G <-  X[,1:p]
G1 <- X[,1:p]

W2 <- X[,(p+r+1):ncol(X)]
# Genotype per interaction matrix
l = s + r
zero_exposure_samples <- data$Sample[!is.na(data$Mutagen) & data$Drug.concentration==0]
for (j in 1:length(zero_exposure_samples)) {
  if (interactions[zero_exposure_samples[j]]=="N2:MMS/EMS/DMS")
    next
  colname <- interactions[zero_exposure_samples[j]]
  colname <- gsub(x = colname, replacement = '.', pattern = '[ ]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[,]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[:]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[-]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[)]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[(]')
  W2[match(zero_exposure_samples[j],data$Sample),colname] <- 1
}
if ("N2.MMS.EMS.DMS" %in% colnames(W2))
  W2 <- W2[,-match("N2.MMS.EMS.DMS" , colnames(W2))]

X <- X[,-grep('N2.',colnames(X),fixed=T)]
X <- X[,colSums(X)>0]
W <- X[,(p+r+1):ncol(X)]

M1 <- X[,(p+1):(p+r)]
M1[M1>0] <- 1

Mall <- X[,(p+1):(p+r)] * doses
avdose <- NULL
for (j in 1:ncol(Mall)) {
  avdose <- c(avdose, mean(Mall[Mall[,j]>0,j]))
  Mall[,j] = Mall[,j] / mean(Mall[Mall[,j]>0,j])
}
names(avdose) <- colnames(Mall)
# new doses
doses = t(t(rowSums(Mall)))

EMS_batch <- as.numeric(data$Comments)
EMS_batch[is.na(EMS_batch)] <- 0

#########################################################################################################

# Run the model

##################################

# 1. short model

short.Y <- cbind(rowSums(Y[,1:16]),rowSums(Y[,17:32]),rowSums(Y[,33:48]),rowSums(Y[,49:64]),rowSums(Y[,65:80]),rowSums(Y[,81:96]),
                 rowSums(Y[,97:98]),rowSums(Y[,99:104]),rowSums(Y[,105:106]),rowSums(Y[,107:112]),rowSums(Y[,113:119]))
Y <- short.Y

m = ncol(Y)
r = ncol(Mall)
n = nrow(Y)
p = ncol(G1)

sigma_G <- gamma(shape = 1, rate = 1, dim = c(1,p))
beta_GH <- lognormal(meanlog = matrix(0, nrow = m, ncol = p), sdlog = matrix(1, nrow = m, ncol = 1) %*% sigma_G, dim = c(m, p))

sigma_M <- gamma(shape = 1, rate = 1, dim = c(1,r))
beta_M = lognormal(meanlog = matrix(0,nrow = m,ncol=r), sdlog = matrix(1,nrow=m,ncol=1) %*% sigma_M, dim = c(m, r))

# generations
sigma_G = 0.5
alpha_G = normal(mean=0, sd = sigma_G, dim = c(l,1), truncation = c(0,Inf))
r_G = W2 %*% alpha_G

# genotype aplification by mutagens
sigma_G_M = 0.5
alpha_G_M = normal(mean = 0, sd = sigma_G_M, dim=c(s,1), truncation = c(0,10))
r_GM = W %*% alpha_G_M

# EMS
EMS_batch[EMS_batch==1] <- 0
W_EMS <- matrix(0,nrow = n, ncol = 4)
for (j in 1:4)
  W_EMS[which(EMS_batch == j+1),j] <- 1
EMS_dose_adjustment = normal(mean = 0, sd = 0.5, dim = c(4,1))
r_doses = exp(W_EMS %*% EMS_dose_adjustment)

# Interactions
sigma2_I_2 = gamma(shape = 1, rate = 1, dim = c(1,s))
beta_I_all = greta::laplace(matrix(0, nrow = m, ncol = s),
			    sigma = matrix(1,nrow=m,ncol=1) %*% sigma2_I_2)

# Model
mu_GH = as.matrix(G1) %*% t(beta_GH)
mu_M = M1 %*% t(beta_M)
mu1 = exp(W %*% t(beta_I_all)) 
mu =  mu_GH * ((g + r_G + (doses * r_doses)*r_GM) %*% matrix(1,nrow=1,ncol=m)) + (mu1 * mu_M) * ((doses * r_doses) %*% matrix(1,nrow=1,ncol=m))
prob = size / (size + mu)
distribution(Y) = negative_binomial(size = size, prob = prob)

model_full <- model(beta_GH, beta_M, beta_I_all, alpha_G, alpha_G_M, EMS_dose_adjustment, sigma_M)
draws.step1 <- mcmc(model_full, warmup = 2000, n_samples = 2000)
draws_all <- do.call('rbind',draws.step1)

beta_I_short <- matrix(colMeans(draws_all[,grep('beta_I',colnames(draws.step2[[1]]))]), nrow = m, ncol = s)
beta_I_short_var <- matrix(apply(draws_all[,grep('beta_I',colnames(draws.step2[[1]]))],2,var), nrow = m, ncol = s)
beta_I_short_low <- matrix(apply(draws_all[,grep('beta_I',colnames(draws.step2[[1]]))],2,quantile,0.025), nrow = m, ncol = s)
beta_I_short_high <- matrix(apply(draws_all[,grep('beta_I',colnames(draws.step2[[1]]))],2,quantile,0.975), nrow = m, ncol = s)
colnames(beta_I_short) = colnames(beta_I_short_var) = colnames(beta_I_short_low) = colnames(beta_I_short_high) = colnames(W)

short_draws <- draws_all[, c(grep('beta_G',colnames(draws_all)),
grep('beta_M',colnames(draws_all)),
grep('beta_I',colnames(draws_all)),grep('alpha_G_M',colnames(draws_all)))]

changes_subs <- list()
changes_indels <- list()
changes_sv <- list()

k <- 1
unit <- c(1, 10, 10, 100, 0.1, 100, 1, 10, 1, 100, 10, 1)
#  1 microM, 10 microM, 10 microM, 100 Gray, 100 microM, 100 Gray, 1 mM, 10 milliM, 1 mM, 100 Joule, 10 microM, 1 mM
unit_name <- c('muM','muM','muM','Gy', 'muM', 'Gy', 'mM', 'mM', 'mM','J/m2', 'muM','mM')
names(unit_name) = names(unit) = names(avdose) <- colnames(Mall)
cosine <- function(x,y) {
  return(sum(x * y) / sqrt(sum(x**2)) / sqrt(sum(y**2)))
}
for (j in sample(c(1:8000),1000)) {
  beta_GH_greta_full_tmp <- matrix(as.matrix(short_draws[j,grep('beta_G',colnames(draws))]), nrow = m, ncol = p)
  beta_M_greta_full_tmp <- matrix(as.matrix(short_draws[j,grep('beta_M',colnames(draws))]), nrow = m, ncol = r)
  alpha_GM_greta_tmp <- as.matrix(short_draws[j,grep('alpha_G_M',colnames(draws))])
  beta_I_greta_full_tmp <- matrix(as.matrix(short_draws[j,grep('beta_I',colnames(draws))]), nrow = m, ncol = s)
  changes_subs[[j]] <- rep(NA,ncol(W))
  names(changes_subs[[j]]) <- colnames(W)
  changes_indels[[j]] <- rep(NA,ncol(W))
  names(changes_indels[[j]]) <- colnames(W)
  changes_sv[[j]] <- rep(NA,ncol(W))
  names(changes_sv[[j]]) <- colnames(W)
  for (z in colnames(W)) {
    zg = paste(unlist(strsplit(z,split='[.]'))[1:2],collapse='.')
    if (substr(zg,1,3) == 'bub') {
      zg <- substr(z,1,14)
    }
    if (substr(zg,1,7) == 'rad.54B') {
      zg <- substr(z,1,16)
    }
    zm = unlist(strsplit(z,split='[.]'))[length(unlist(strsplit(z,split='[.]')))]
    if (zm == 'B1') zm <- 'Aflatoxin.B1'
    mu_0 <- beta_GH_greta_full_tmp[,match(zg,colnames(beta_GH_greta_full))] + 
      beta_M_greta_full_tmp[,match(zm,colnames(beta_M_greta_full))] / avdose[zm] * unit[zm]
    mu_1 <- beta_GH_greta_full_tmp[,match(zg,colnames(beta_GH_greta_full))] +#*
      #(1 + alpha_GM_greta_tmp[match(z,colnames(W)),] / avdose[zm] * unit[zm]) +
      beta_M_greta_full_tmp[,match(zm,colnames(beta_M_greta_full))] * 
      exp(beta_I_greta_full_tmp[,match(z,colnames(beta_I_greta_full))]) / avdose[zm] * unit[zm]
    changes_subs[[j]][z] <- (sum(mu_1[1:6]) / sum(mu_0[1:6]))
    changes_indels[[j]][z] <- (sum(mu_1[8:10]) / sum(mu_0[8:10]))
    changes_sv[[j]][z] <- (sum(mu_1[11]) / sum(mu_0[11]))
  }
  print(j)
}
save(changes_subs, changes_indels, changes_sv, file = 'Fold_changes_short_model.RData'))

##################################

# 119-long model

ma.samples <- data$Sample[is.na(data$Mutagen) & data$Generation>0]
ma.X <- X[ma.samples,] * sapply(data$Generation[match(ma.samples, data$Sample)],generation.function)
ma.X <- ma.X[,colSums(ma.X)>0]
ma.Y <- Y[ma.samples,]

sigma_GH <- variable(lower = 0)
beta_GH = lognormal(meanlog = 0, sdlog = sigma_GH, dim = c(m, ncol(ma.X))) # lognormal with fitted variance
mu_GH = (ma.X %*% t(beta_GH))
size = 100
prob = size / (size + mu_GH)
distribution(ma.Y) = negative_binomial(prob = prob, size = size * matrix(1,nrow = nrow(ma.Y), ncol = ncol(ma.Y)))
ma.model <- model(beta_GH,sigma_GH)

ma.draws <- mcmc(ma.model, warmup = 2000, n_samples = 5000, thin = 2) # takes a lot of memory and time - better do on HPC

mcmc_trace(ma.draws[,grep('beta_GH',colnames(ma.draws[[1]]))[sample(1:(m*ncol(ma.X)),9)]])
mcmc_trace(ma.draws[,grep('sigma_GH',colnames(ma.draws[[1]])),drop=F])

beta_GH_greta <- matrix(colMeans((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))])), nrow = m, ncol = ncol(ma.X))
beta_GH_greta_low <- matrix(apply((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]),2,quantile,0.025), nrow = m, ncol = ncol(ma.X))
beta_GH_greta_high <- matrix(apply((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]),2,quantile,0.975), nrow = m, ncol = ncol(ma.X))
beta_GH_greta_var <- matrix(apply((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]),2,var), nrow = m, ncol = ncol(ma.X))
colnames(beta_GH_greta) = colnames(beta_GH_greta_low) = colnames(beta_GH_greta_high) =colnames(beta_GH_greta_var) <- colnames(ma.X)

beta.draws <- do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]
shapes <- matrix(sapply(1:ncol(beta.draws), function(i) { print(i); return(fitdistr(beta.draws[,i],'gamma')[[1]][1]) }), nrow = 119, ncol = 51)
rates <- matrix(sapply(1:ncol(beta.draws), function(i) { print(i); return(fitdistr(beta.draws[,i],'gamma')[[1]][2]) }), nrow = 119, ncol = 51)
# now feed these priors into the full genotypes structure
new_shapes <- matrix(1,nrow = 119, ncol = ncol(G1))
new_rates <- matrix(1,nrow = 119, ncol = ncol(G1))
for (j in 1:ncol(G1)) {
  if (colnames(G1)[j] %in% colnames(ma.X)) {
    new_shapes[,j] <- shapes[,match(colnames(G1)[j],colnames(ma.X))]
    new_rates[,j] <- rates[,match(colnames(G1)[j],colnames(ma.X))]
  }
  else {
    new_shapes[,j] <- shapes[,match('N2',colnames(ma.X))]
    new_rates[,j] <- rates[,match('N2',colnames(ma.X))]
  }
}
# Gamma priors fitted to posterior samples from Step 1
beta_GH <- gamma(shape = new_shapes, rate = new_rates, dim = c(m, ncol(G1)))

sigma_M <- gamma(shape = 1, rate = 1, dim = c(1,r))
beta_M = lognormal(meanlog = matrix(0,nrow = m,ncol=r), sdlog = matrix(1,nrow=m,ncol=1) %*% sigma_M, dim = c(m, ncol(Mall)))

sigma_G = 0.5
alpha_G = lognormal(mean = 0, sd=sigma_G, truncation = c(0,Inf), dim = c(l,1))
r_G = W2 %*% alpha_G

W_EMS <- matrix(0,nrow = n, ncol = 4) # 5 batches total, so 4 adjustments
for (j in 1:ncol(W_EMS))
  W_EMS[which(EMS_batch == j+1),j] <- 1
EMS_dose_adjustment = normal(mean = 0, sd = 0.5, dim = c(4,1))
r_doses = exp(W_EMS %*% EMS_dose_adjustment)

sigma_G_M = 0.5
alpha_G_M = lognormal(mean = 0, sd=sigma_G_M, truncation = c(0,Inf), dim=c(s,1))
r_GM = W %*% alpha_G_M

linemat <- function(n1,n2,l) {
  tmptmp <- matrix(0,nrow=n1,ncol=n2)
  tmptmp[l,] <- 1
  return(tmptmp)
}
tmp = cbind(linemat(11,16,1),linemat(11,16,2),linemat(11,16,3),linemat(11,16,4),linemat(11,16,5),linemat(11,16,6),
            linemat(11,2,7),linemat(11,6,8),linemat(11,2,8),linemat(11,6,8),linemat(11,7,9))

sigma2_I_2 = variable(lower = 0, upper = 10, dim= c(1,s))
beta_I_all = greta::laplace(mu = t(tmp) %*% beta_I_short,
			    sigma = matrix(1,nrow=m,ncol=1) %*% sigma2_I_2)

mu_GH = as.matrix(G1) %*% t(beta_GH) # genetic background
M1 <- Mall>0
mu_M = as.matrix(M1) %*% t(beta_M) # mutagen contribution
mu1 = exp(as.matrix(W) %*% t(beta_I_2)) # multiplicative genotype - mutation interaction
mu =  mu_GH * ((g + r_G + (r_doses*doses) * r_GM) %*% matrix(1,nrow=1,ncol=m)) +
  mu1 * mu_M * ((doses*r_doses) %*% matrix(1,nrow=1,ncol=m))
size = 100 # account for a bit of overdispersion
prob = size / (size + mu)
distribution(Y) = negative_binomial(size = size, prob = prob) # fit negative binomial distribution to the counts
model_full <- model(beta_M, beta_GH, sigma_M, beta_I_2, sigma2_I_2, beta_GH, alpha_G, alpha_G_M, EMS_dose_adjustment) 

draws <- mcmc(model_full, warmup = 2000, n_samples = 5000, thin = 2)

#########################################################################################################

# Extract coefficients

# Some convergence diagnostics
mcmc_trace(draws[,grep('beta_M',colnames(draws[[1]]))[sample(1:(r*m),16)]])
mcmc_trace(draws[,grep('beta_I',colnames(draws[[1]]))[sample(1:(s*m),16)]])
mcmc_trace(draws[,grep('alpha_G',colnames(draws[[1]]))[sample(1:l,16)]])
mcmc_trace(draws[,grep('alpha_G_M',colnames(draws[[1]]))[sample(1:s,16)]])
mcmc_trace(draws[,grep('sigma2_I_2',colnames(draws[[1]]))[sample(1:s,16)],drop=F])
mcmc_trace(draws[,grep('sigma_M',colnames(draws[[1]])),drop=F])
mcmc_trace(draws[,grep('EMS',colnames(draws[[1]]))[sample(1:sum(EMS_unknown_batch),16)]])


draws_all <- do.call('rbind',draws)
alpha_G_greta <- colMeans(draws_all[,grep('alpha_G',colnames(draws_all))[1:l]])
alpha_G_greta_low <- apply(draws_all[,grep('alpha_G',colnames(draws_all))[1:l]],2,quantile,0.025)
alpha_G_greta_high <- apply(draws_all[,grep('alpha_G',colnames(draws_all))[1:l]],2,quantile,0.975)
alpha_GM_greta <- colMeans(draws_all[,grep('alpha_G_M',colnames(draws_all))])
alpha_GM_greta_low <- apply(draws_all[,grep('alpha_G_M',colnames(draws_all))],2,quantile,0.025)
alpha_GM_greta_high <- apply(draws_all[,grep('alpha_G_M',colnames(draws_all))],2,quantile,0.975)
alpha_GM_greta_var <- apply(draws_all[,grep('alpha_G_M',colnames(draws_all))],2,var)

beta_GH_greta_full <- matrix(colMeans(exp(draws_all[,grep('beta_GH',colnames(draws_all))])), nrow = m, ncol = p)
beta_GH_greta_full_var <- matrix(apply(exp(draws_all[,grep('beta_GH',colnames(draws_all))]),2,var), nrow = m, ncol = p)
beta_GH_greta_full_low <- matrix(apply(exp(draws_all[,grep('beta_GH',colnames(draws_all))]),2,quantile,0.025), nrow = m, ncol = p)
beta_GH_greta_full_high <- matrix(apply(exp(draws_all[,grep('beta_GH',colnames(draws_all))]),2,quantile,0.975), nrow = m, ncol = p)
colnames(beta_GH_greta_full) <- colnames(G1)
colnames(beta_GH_greta_full_low) <- colnames(G1)
colnames(beta_GH_greta_full_high) <- colnames(G1)

beta_M_greta_full <- matrix(colMeans(draws_all[,grep('beta_M',colnames(draws_all))]), nrow = m, ncol = r)
beta_M_greta_full_var <- matrix(apply(draws_all[,grep('beta_M',colnames(draws_all))],2,var), nrow = m, ncol = r)
beta_M_greta_full_low <- matrix(apply(draws_all[,grep('beta_M',colnames(draws_all))],2,quantile,0.025), nrow = m, ncol = r)
beta_M_greta_full_high <- matrix(apply(draws_all[,grep('beta_M',colnames(draws_all))],2,quantile,0.975), nrow = m, ncol = r)

colnames(beta_M_greta_full) = colnames(beta_M_greta_full_var) = colnames(beta_M_greta_full_high) = colnames(beta_M_greta_full_low) = colnames(Mall)

beta_I_greta_full <- matrix(colMeans(draws_all[,grep('beta_I',colnames(draws_all))]), nrow = m, ncol = s)
beta_I_greta_full_var <- matrix(apply(draws_all[,grep('beta_I',colnames(draws_all))],2,var), nrow = m, ncol = s)
beta_I_greta_full_low <- matrix(apply(draws_all[,grep('beta_I',colnames(draws_all))],2,quantile,0.025), nrow = m, ncol = s)
beta_I_greta_full_high <- matrix(apply(draws_all[,grep('beta_I',colnames(draws_all))],2,quantile,0.975), nrow = m, ncol = s)

colnames(beta_I_greta_full) = colnames(beta_I_greta_full_var)=colnames(beta_I_greta_full_low) = colnames(beta_I_greta_full_high) = colnames(W)

EMS_dose_adjustment_greta = colMeans(draws_all[,grep('EMS',colnames(draws.step2[[1]]))])
EMS_dose_adjustment_greta_low = apply(draws_all[,grep('EMS',colnames(draws.step2[[1]]))],2,quantile,0.025)
EMS_dose_adjustment_greta_high = apply(draws_all[,grep('EMS',colnames(draws.step2[[1]]))],2,quantile,0.975)


mu <- (as.matrix(G1) %*% t((beta_GH_greta_full))) * 
  ((g + as.matrix(W2) %*% t(t(alpha_G_greta)) + (doses * exp(W_EMS %*% t(t(EMS_dose_adjustment_greta)))) * (as.matrix(W) %*% t(t(alpha_GM_greta)))) 
   %*% matrix(1,nrow = 1,ncol=119)) +
  (as.matrix(M1) %*% t((beta_M_greta_full))) * exp(as.matrix(W) %*% t(beta_I_greta_full)) *
  ((doses * exp(W_EMS %*% t(t(EMS_dose_adjustment_greta)))) %*% matrix(1,nrow = 1,ncol=119))
plot(rowSums(Y), rowSums(mu),pch = 16)
abline(a=0,b=1,col='red')

draws_full <- draws_all[, c(grep('beta_G',colnames(draws_all)), grep('beta_M',colnames(draws_all)), grep('beta_I',colnames(draws_all)),grep('alpha_G_M',colnames(draws_all)))]

similarities <- list()
for (j in 1:nrow(draws)) {
  beta_GH_greta_full_tmp <- matrix(as.matrix(draws_full[j,grep('beta_G',colnames(draws))]), nrow = m, ncol = p)
  beta_M_greta_full_tmp <- matrix(as.matrix(draws_full[j,grep('beta_M',colnames(draws))]), nrow = m, ncol = r)
  alpha_GM_greta_tmp <- as.matrix(draws_full[j,grep('alpha_G_M',colnames(draws))])
  beta_I_greta_full_tmp <- matrix(as.matrix(draws_full[j,grep('beta_I',colnames(draws))]), nrow = m, ncol = s)
  names(changes_sv[[j]]) <- colnames(W)
  similarities[[j]] <- rep(NA,ncol(W))
  names(similarities[[j]]) <- colnames(W)
  for (z in colnames(W)) {
    zg = paste(unlist(strsplit(z,split='[.]'))[1:2],collapse='.')
    if (substr(zg,1,3) == 'bub') {
      zg <- substr(z,1,14)
    }
    if (substr(zg,1,7) == 'rad.54B') {
      zg <- substr(z,1,16)
    }
    zm = unlist(strsplit(z,split='[.]'))[length(unlist(strsplit(z,split='[.]')))]
    if (zm == 'B1') zm <- 'Aflatoxin.B1'
    mu_0 <- beta_GH_greta_full_tmp[,match(zg,colnames(beta_GH_greta_full))] + 
      beta_M_greta_full_tmp[,match(zm,colnames(beta_M_greta_full))] / avdose[zm] * unit[zm]
    mu_1 <- beta_GH_greta_full_tmp[,match(zg,colnames(beta_GH_greta_full))] *
      (1 + alpha_GM_greta_tmp[match(z,colnames(W)),] / avdose[zm] * unit[zm]) +
      beta_M_greta_full_tmp[,match(zm,colnames(beta_M_greta_full))] * 
      exp(beta_I_greta_full_tmp[,match(z,colnames(beta_I_greta_full))]) / avdose[zm] * unit[zm]
    similarities[[j]][z] <- cosine(mu_0, mu_1)
  }
  print(j)
}
# similarities <- similarities[, sample(1:nrow(similarities),1000)] # to make it smaller
save(similarities, file = 'Changes_in_full_model.RData'))


#########################################################################################################

# Save the matrices with coefficients

rownames(beta_GH_greta_full) <- colnames(Y)
rownames(beta_GH_greta_full_low) <- colnames(Y)
rownames(beta_GH_greta_full_high) <- colnames(Y)
rownames(beta_M_greta_full) <- colnames(Y)
rownames(beta_M_greta_full_low) <- colnames(Y)
rownames(beta_M_greta_full_high) <- colnames(Y)
names(alpha_G_greta) <- colnames(W2)
names(alpha_G_greta_low) <- colnames(W2)
names(alpha_G_greta_high) <- colnames(W2)
names(alpha_G_greta) <- colnames(W)
names(alpha_G_greta_low) <- colnames(W)
names(alpha_G_greta_high) <- colnames(W2)
names(EMS_dose_adjustment_greta) = names(EMS_dose_adjustment_greta_low) = names(EMS_dose_adjustment_greta_high) <- paste0('batch',2:5)
l <- list(beta_GH_greta_full, beta_GH_greta_full_low, beta_GH_greta_full_high,
          beta_M_greta_full, beta_M_greta_full_low, beta_M_greta_full_high,
          data.frame(mean=alpha_G_greta,low=alpha_G_greta_low,high=alpha_G_greta_high),
          data.frame(mean=EMS_dose_adjustment_greta,low=EMS_dose_adjustment_greta_low,high=EMS_dose_adjustment_greta_high))
openxlsx::write.xlsx(l, file="~/Supplementary Table 2. Signatures of DNA repair deficiencies and mutagen exposures.xlsx",
                     colNames = T, rowNames = T,
                     sheetName = c('Genotypes-mean', "Genotypes-low","Genotypes-high",
                                   'Mutagens-mean', "Mutagens-low","Mutagens-high",
                                   'Generation_adjustment', 'EMS_dose_adjustment'))

rownames(beta_I_greta_full) <- colnames(Y)
rownames(beta_I_greta_full_low) <- colnames(Y)
rownames(beta_I_greta_full_high) <- colnames(Y)
names(alpha_GM_greta) <- colnames(W)
names(alpha_GM_greta_low) <- colnames(W)
names(alpha_GM_greta_high) <- colnames(W)
l <- list(beta_I_greta_full, beta_I_greta_full_low, beta_I_greta_full_high,
          data.frame(mean=alpha_GM_greta,low=alpha_GM_greta_low,high=alpha_GM_greta_high))
openxlsx::write.xlsx(l, file="~/Supplementary Table 3. Interaction effects beta_I on the mutagen signatures and alpha_GM on genotype signatures.xlsx",
                     colNames = T, rowNames = T,
                     sheetName = c('Log-Beta-means', "Log-Beta-low","Log-Beta-high",'Dose-dependent b_I'))
