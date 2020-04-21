################################################################################
##### Supplementary functions                                              #####
##### January 2016 - ... , EMBL-EBI, N.Volkova (nvolkova@ebi.ac.uk)        #####
################################################################################

# A set of useful functions

# round all numeric columns of a dataframe
round_dataframe <- function(df, n = 2) {
  inds <- which(sapply(df,is.numeric))
  for (j in inds) df[,j] <- round(df[,j],n)
  return(df)
}

# return the top n elements of a vector
top <- function(x,n) {
  return(sort(x,decreasing=T)[1:n])
}
# return indices of top n elements of a vector
which.top <- function(x,n) {
  if (is.null(names(x))) return(order(x,decreasing=T)[1:n])
  return(names(x)[order(x,decreasing=T)][1:n])
}

# return the bottom n elements of a vector
bottom <- function(x,n) {
  return(sort(x,decreasing=F)[1:n])
}
# return indices of the bottom n elements of a vector
which.bottom <- function(x,n) {
  if (is.null(names(x))) return(order(x,decreasing=F)[1:n])
  return(names(x)[order(x,decreasing=F)][1:n])
}

# cosine similarity
cosine <- function(x,y) {
    return(sum(x * y) / sqrt(sum(x**2)) / sqrt(sum(y**2)))
}

# KL-divergence for real values
divergence <- function (a,b) {
  return (sum(a * log ( (a+.Machine$double.eps)/(b + .Machine$double.eps)) - a + b))
}

# KL divergence for probabilities
KL <- function(p,q) {
  return(sum(p * log2((p+0.001)/(q+0.001))))
}

# Jensen-Shaennon distance
JSdistance <- function(p,q) {
  M = 0.5*(p+q)
  return (sqrt(0.5*KL(p,M) + 0.5*KL(q,M)))
}

# read in a set of VCF files
read_ce_vcf <- function(file, genome = "WBcel235") {
  out <- tryCatch(
    {
      readVcf(file, genome=genome) 
    },
    error=function(cond) {
      message(paste("\nBad file in: ", file))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("\nFile caused a warning: ", file))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}

# reverse complement
RevCom <- function(x) {
  return(as.character(reverseComplement(DNAString(x))))
}

# plot all a signature and all its chances
plot_all_changes <- function(plot_fc,S_est,beta_est,n = F) {
  to.show <- data.frame(t(S_est))
  for (j in 1:nrow(beta_est)) {
    to.show <- cbind(to.show, S_est[nrow(S_est),] * exp(beta_est[j,]))
  }
  colnames(to.show) <- c(paste0('Sig',1:nrow(S_est)), paste0('Sig*beta',1:nrow(beta_est)))
  plot_fc(to.show, norm = n)
}

# collect trinucleotide substitutions from a C. elegans vcf
getTrinucleotideSubs <- function(vcf, ref_genome="BSgenome.Celegans.UCSC.ce11") {
  seqlevels(vcf) <- seqnames(get(ref_genome))
  if (ref_genome==worm_ref_genome)
    seqlevels(vcf) <- c("chrM","chrIV","chrIII","chrX","chrI","chrV","chrII")
  vcf <- vcf[as.character(seqnames(vcf)) != "chrM"]
  ranges = resize(vcf, 3, fix = "center")
  tnc = as.character(getSeq(get(ref_genome), seqnames(vcf), 
                            start(vcf) - 1, end(vcf) + 1))
  s <- paste0(substr(tnc,1,1),"[",as.character(ref(vcf)), ">",
              as.character(unlist(alt(vcf))),"]", substr(tnc,3,3))
  n <- c("A","C","G","T")
  f <- paste0(rep(n, each=4), "[", rep(n, each=96/2), ">", c(rep(c("C","G","T"), each=48/3),rep(c("A","G","T"), each=48/3),rep(c("A","C","T"), each=48/3), rep(c("A","C","G"), each=48/3)), "]", n)
  s <- factor(s, levels=f)
  return(s)
}

# collect trinucleotide substitutions from human vcf file
getTrinucleotideSubs_human_table <- function(subs, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19', chr = 0) {
  g = get(ref_genome)
  if (is.na(subs)) return(NA)
  if (chr) subs <- subs[subs$Chrom != "chrM",] else subs <- subs[subs$Chrom != "MT",]
  if (nrow(subs)==0) return(NA)
  if (chr)
    tnc = as.character(getSeq(g, subs$Chrom, subs$Pos - 1, subs$Pos + 1))
  else
    tnc = as.character(getSeq(g, paste('chr',subs$Chrom,sep=''), 
                            subs$Pos - 1, subs$Pos + 1))
  s <- paste0(substr(tnc,1,1),"[",as.character(subs$Ref), ">",
              as.character(subs$Alt),"]", substr(tnc,3,3))
  n <- c("A","C","G","T")
  f <- paste0(rep(n, each=4), "[", 
              rep(n, each=96/2), ">", 
              c(rep(c("C","G","T"), each=16),
                rep(c("A","G","T"), each=16),
                rep(c("A","C","T"), each=16), 
                rep(c("A","C","G"), each=16)),
              "]", n)
  s <- factor(s, levels=f)
  return(s)
}

# adjust trinucleotide variants to pyrimidine reference
tncToPyrimidine <- function(nucl) {
  ind <- sort(c(grep('[G',nucl,fixed = T),grep('[A',nucl,fixed = T)))
  newnucl <- as.character(reverseComplement(DNAStringSet(paste(substr(nucl[ind],1,1),
                                                               substr(nucl[ind],5,5),
                                                               substr(nucl[ind],3,3),
                                                               substr(nucl[ind],7,7),sep=''))))
  nucl[ind] <- paste(substr(newnucl,1,1),'[',
                     substr(newnucl,2,2),'>',
                     substr(newnucl,3,3),']',
                     substr(newnucl,4,4),sep='')
  return(nucl)
}

# reverse complement with mutation type
RevComMutType <- function(nucl) {
  tmp <- unlist(strsplit(nucl, split = ''))
  return(paste0(RevCom(tmp[7]),'[',RevCom(tmp[3]),'>',RevCom(tmp[5]),']',RevCom(tmp[1])))
}

# collect dinucleotide variant types from VCF
get_DNV_type <- function(vcf, ref_genome="BSgenome.Celegans.UCSC.ce11", dnv.types) {
  s <- paste0(paste(as.character(ref(vcf)),collapse=''), ">", paste(as.character(unlist(alt(vcf))), collapse=''))
  if (!(s %in% dnv.types)) {
    s <- paste0(as.character(reverseComplement(DNAString(paste(as.character(ref(vcf)),collapse='')))), ">",
                as.character(reverseComplement(DNAString(paste(as.character(unlist(alt(vcf))), collapse='')))))
  }
  return(s)
}

# check if a set of variants is MNV or not
isMNV <- function(vcf) {
  d <- diff(start(vcf)) == 1 & abs(diff(geno(vcf)$PM[,"TUMOUR"] )) <= 0.05
  w <- c(FALSE, d) | c(d, FALSE)
  return(w)
}

# check if variants are in a repeat or not
is.in.repeat <- function(vcf, ref_genome="BSgenome.Celegans.UCSC.ce11") {
  seqlevels(vcf) <- c("chrM","chrIV","chrIII","chrX","chrI","chrV","chrII")
  vcf <- vcf[as.character(seqnames(vcf)) != "chrM"]
  #ranges = resize(vcf, 3, fix = "center")
  in.repeat <- sapply(1:length(vcf), function(j) {
    tmp <- vcf[j]
    if (info(tmp)$PC=="DI") return(FALSE)
    tnc = as.character(getSeq(get(ref_genome), seqnames(tmp), 
                              start(tmp) - info(vcf)$LEN[j]*5,
                              end(tmp) + info(vcf)$LEN[j]*5))
    if (info(tmp)$PC=="I")
      return(grepl(paste0('(',substr(as.character(unlist(alt(tmp))), 2, nchar(as.character(unlist(alt(tmp))))),'){3,}'), tnc))
    if (info(tmp)$PC=="D")
      return(grepl(paste0('(',substr(as.character(ref(tmp)), 2, nchar(as.character(ref(tmp)))),'){3,}'), tnc))
    return(FALSE)
  })
  return(in.repeat)
}

# derive minimal range from a set of ranges
min_ranges <- function(gr) {
  if (length(gr) == 1) return(gr)
  to.merge <- 1
  while (sum(to.merge)!=0) {
    gr_start <- start(gr)
    gr_end <- end(gr)
    m = 1
    to.merge <- vector('numeric',length(gr))
    for (j in 2:length(gr)) {
      if (gr_start[j] < gr_end[j-1]) {
        if (gr_end[j] <= gr_end[j-1]) to.merge[j] <- -1
        else { 
          to.merge[j] <- m
          to.merge[j-1] <- m
          m <- m+1
        }
      }
    }
    if (sum(to.merge==-1)>0) {
      gr <- gr[-which(to.merge == -1)]
      to.merge <- to.merge[-which(to.merge == -1)]
    }
    if (sum(to.merge>0)>0) {
      for (j in 1:max(to.merge)) {
        start(gr[to.merge == j]) <- min(gr_start[to.merge == j])
        end(gr[to.merge == j]) <- min(gr_end[to.merge == j])
      }
    }
    gr <- unique(gr)
  }
  
  return(gr)
}

# GRanges overlap over 50%
`%over.50%` <- function(query, subject) overlapsAny(query, subject, minoverlap = 0.5 * min(width(query),width(subject)))

# filtering breakpoints
filter.breakpoints <- function(SV, reference) {
  
  reference$CHR1 <- as.character(reference$CHR1)
  reference$CHR2 <- as.character(reference$CHR2)
  reference$TYPE <- as.character(reference$TYPE)
  
  if (is.na(SV) || nrow(SV)==0) return(SV)
  
  bad <- NULL
  
  SV$CHR1 <- as.character(SV$CHR1)
  SV$CHR1 <- as.character(SV$CHR1)
  SV$TYPE <- as.character(SV$TYPE)
  
  `%over.50%` <- function(query, subject) overlapsAny(query, subject, minoverlap = 0.5 * min(width(query),width(subject)))
  
  for (j in 1:nrow(SV)) {
    
    tmp <- reference[reference$CHR1 == SV$CHR1[j] & 
                       reference$CHR2 == SV$CHR2[j] & 
                       reference$TYPE == SV$TYPE[j],,drop = F]
    if (nrow(tmp) > 0) {
      big.ranges.1 <- IRanges(start=as.numeric(tmp$POS1)-200,end=as.numeric(tmp$POS1)+200)
      big.ranges.2 <- IRanges(start=as.numeric(tmp$POS2)-200,end=as.numeric(tmp$POS2)+200)
      inds <- intersect(which(big.ranges.1 %over.50% IRanges(start=as.numeric(SV$POS1[j])-200,end=as.numeric(SV$POS1[j])+200)),
                        which(big.ranges.2 %over.50% IRanges(start=as.numeric(SV$POS2[j])-200,end=as.numeric(SV$POS2[j])+200)))
      if (length(inds) > 0)
        bad <- c(bad,j)
    }
  }
  if (length(bad)>0)
    return(SV[-bad,,drop=F])
  return(SV)
}

# Poisson regression for signature fitting
nmSolve <- function(D, P, maxIter = 10000, tol=1e-5, div.err=1e-7) {
  n <- nrow(D)
  mask <- !is.na(D)
  m <- ncol(D)
  s <- ncol(P)
  rP <- rep(colSums(P), m)
  tP <- t(P)
  D <- as.matrix(D)
  P <- as.matrix(P)
  E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
  err <- 2*tol
  D[is.na(D)] <- 0
  
  iter <- 1
  divergence.old <- mean(D*log(D/(P %*% (E2 + .Machine$double.eps))) - D + P%*%E2, na.rm=T)
  div.change <- 2 * div.err
  
  while (iter < maxIter & err > tol & abs(div.change) > div.err) {
    E1 <- E2
    E2 <- E1 * (tP %*% ((mask*D)/(mask*(P %*% (E1 + .Machine$double.eps)) + .Machine$double.eps)))/rP
    iter <- iter + 1
    err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
    divergence <- mean(D*log(D/(P %*% (E2 + .Machine$double.eps))) - D + P%*%E2, na.rm=T) # KL distance from D to P%*%E2
    div.change <- divergence.old - divergence
    divergence.old = divergence
    if(iter %% 100 == 0) cat(round(-log10(err)))
  }
  cat("\n")
  if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
  E2
}
