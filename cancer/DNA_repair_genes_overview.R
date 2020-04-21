library(VariantAnnotation)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(openxlsx)
source('../useful_functions.R')

genes <- openxlsx::read.xlsx('../Supplementary_tables/Supplement/Supplementary Table 4. TCGA samples with mutations in DNA pathways.xlsx', sheet = 1)
# Get the new signature set from here: https://www.synapse.org/#!Synapse:syn11967914
new_cancer_signatures <- read.csv('sigProfiler_exome_SBS_signatures.csv')

library(BSgenome.Hsapiens.UCSC.hg19)
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library("org.Hs.eg.db")
library(deconstructSigs)
genome <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
hg19.genes <- genes(txdb, columns="gene_id")
gene_symbol <- AnnotationDbi::select(org.Hs.eg.db, keys=hg19.genes$gene_id, 
                                     columns="SYMBOL", keytype="ENTREZID")
hg19.genes$gene_id <- gene_symbol$SYMBOL
txs <- exonsBy(txdb, "gene")
trxs <- transcriptsBy(txdb, "gene")
ALL_alt <- as.character(genes$Gene)
ALL_alt[ALL_alt == "SHFM1"] <- "SEM1"
ALL_alt[ALL_alt == "TCEB1"] <- "ELOB"
ALL_alt[ALL_alt == "TCEB2"] <- "ELOC"
ALL_alt[ALL_alt == "TCEB3"] <- "ELOA"
ALL_alt[ALL_alt == "APITD1"] <- "CENPS"
ALL_alt[ALL_alt == "BRE"] <- "BABAM2"
ALL_alt[ALL_alt == "STRA13"] <- "CENPX"
ALL_alt[ALL_alt == "STRA13"] <- "CENPX"
ALL_alt[ALL_alt == "MRE11A"] <- "MRE11"
ALL_alt[ALL_alt == "FAM175A"] <- "ABRAXAS1"

ALL_odd <- names(hg19.genes[match(ALL_alt, hg19.genes$gene_id)])
names(ALL_odd) <- ALL_alt
transcr <- txs[ALL_odd]


PATH_TO_TCGA_MATRIX='path_to_table_with_mutation_counts_for_TCGA'
bigmat <- read.table(PATH_TO_TCGA_MATRIX,sep='\t')

PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)

# Need the copy number segments for all samples from TCGA
PATH_TO_COPYNUMBER_FOR_TCGA='path_to_table_with_copynumbers_in_samples'
copynumber <- read.delim(PATH_TO_COPYNUMBER_FOR_TCGA, sep = '\t', header = T)
colnames(copynumber) <- c('Tissue','Sample','Chromosome','Start','End','MajorCN','MinorCN','Cancer_type')
copynumber <- GRanges(copynumber)
seqlevels(copynumber) <- paste0('chr',seqlevels(copynumber))

# Methylation list created in separate script methylation_check.R
load('methylation.list.RData')
methylated <- methylated[sapply(methylated,nrow)>0]

PATHTOMUTATION='/path/to/mutation/tables/per/gene/' # to make life easier, I generated mutation tables per gene across all the cancers in TCGA

samples <- list()
for (gene in c(genes$Gene[genes$CORE],
               'POLE','POLD1')) {
  
  methylation <- NA
  
  if (gene == 'SLX1A') gene <- 'SLX1'
  
  
  polk_subs <- read.delim(paste0(PATHTOMUTATION,gene,'_mutations.txt'), sep = '\t', header = T)
  polk_dels <- read.table(paste0(PATHTOMUTATION,gene,'_indels.txt'), sep = '\t', header = T, comment.char = "")
  polk_subs <- polk_subs[polk_subs$Type == 'Sub',c('Sample', 'Chrom', 'Pos','Ref','Alt','Gene','Type','Effect','PM.Norm','PM.Tum')]
  polk_dels$PM.Tum <- (polk_dels$PU.Tum + polk_dels$NU.Tum) / (polk_dels$PR.Tum + polk_dels$NR.Tum)
  polk_dels$PM.Norm <- (polk_dels$PU.Norm + polk_dels$NU.Norm) / (polk_dels$PR.Norm + polk_dels$NR.Norm)
  polk_dels <- polk_dels[,c('Sample', 'Chrom', 'Pos','Ref','Alt','Gene','Type','Effect','PM.Norm','PM.Tum')]
  polk_mutations <- rbind(polk_subs, polk_dels)
  polk_mutations$Sample <- paste0('TCGA',substr(as.character(polk_mutations$Sample),5,
                                              nchar(as.character(polk_mutations$Sample))))
  if (gene == 'SLX1') gene <- 'SLX1A'
  polk_mutations <- polk_mutations[polk_mutations$PM.Tum>0.4 &
                                     polk_mutations$Gene == gene &
                                     polk_mutations$Effect %in% c('missense','nonsense','stop_lost',
                                                                  'start_lost','frameshift','ess_splice'),]
  
  if (gene == 'POLE') {
    polk_subs <- read.delim(paste0(PATHTOMUTATION,gene,'_mutations.txt'), sep = '\t', header = T)
    polk_dels <- read.table(paste0(PATHTOMUTATION,gene,'_indels.txt'), sep = '\t', header = T, comment.char = "")
    polk_subs <- polk_subs[polk_subs$Type == 'Sub',c('Sample', 'Chrom', 'Pos','Ref','Alt','Gene','Type','Effect','PM.Norm','PM.Tum','Protein')]
    polk_dels$PM.Tum <- (polk_dels$PU.Tum + polk_dels$NU.Tum) / (polk_dels$PR.Tum + polk_dels$NR.Tum)
    polk_dels$PM.Norm <- (polk_dels$PU.Norm + polk_dels$NU.Norm) / (polk_dels$PR.Norm + polk_dels$NR.Norm)
    polk_dels <- polk_dels[,c('Sample', 'Chrom', 'Pos','Ref','Alt','Gene','Type','Effect','PM.Norm','PM.Tum','Protein')]
    polk_mutations <- rbind(polk_subs, polk_dels)
    polk_mutations$Sample <- paste0('TCGA',substr(as.character(polk_mutations$Sample),5,
                                                  nchar(as.character(polk_mutations$Sample))))
    polk_mutations <- polk_mutations[polk_mutations$PM.Tum>0.2 &
                                       polk_mutations$Gene == gene &
                                       polk_mutations$Effect %in% c('missense','nonsense','stop_lost',
                                                                    'start_lost','frameshift','ess_splice'),]
    polk_mutations <- polk_mutations[which(as.numeric(substr(as.character(polk_mutations$Protein),4,6)) < 472 & 
                                            as.numeric(substr(as.character(polk_mutations$Protein),4,6)) > 267),]
  }
  if (gene == 'POLD1') {
    polk_subs <- read.delim(paste0(PATHTOMUTATION,gene,'_mutations.txt'), sep = '\t', header = T)
    polk_dels <- read.table(paste0(PATHTOMUTATION,gene,'_indels.txt'), sep = '\t', header = T, comment.char = "")
    polk_subs <- polk_subs[polk_subs$Type == 'Sub',c('Sample', 'Chrom', 'Pos','Ref','Alt','Gene','Type','Effect','PM.Norm','PM.Tum','Protein')]
    polk_dels$PM.Tum <- (polk_dels$PU.Tum + polk_dels$NU.Tum) / (polk_dels$PR.Tum + polk_dels$NR.Tum)
    polk_dels$PM.Norm <- (polk_dels$PU.Norm + polk_dels$NU.Norm) / (polk_dels$PR.Norm + polk_dels$NR.Norm)
    polk_dels <- polk_dels[,c('Sample', 'Chrom', 'Pos','Ref','Alt','Gene','Type','Effect','PM.Norm','PM.Tum','Protein')]
    polk_mutations <- rbind(polk_subs, polk_dels)
    polk_mutations$Sample <- paste0('TCGA',substr(as.character(polk_mutations$Sample),5,
                                                  nchar(as.character(polk_mutations$Sample))))
    polk_mutations <- polk_mutations[polk_mutations$PM.Tum>0.4 &
                                       polk_mutations$Gene == gene &
                                       polk_mutations$Effect %in% c('missense','nonsense','stop_lost',
                                                                    'start_lost','frameshift','ess_splice'),]
    polk_mutations <- polk_mutations[which(as.numeric(substr(as.character(polk_mutations$Protein),4,6)) < 518 & 
                                             as.numeric(substr(as.character(polk_mutations$Protein),4,6)) > 303),]
  }
  
  mutated <- substr(unique(as.character(polk_mutations$Sample)),1,12)
  very_mutated <- substr(unique(as.character(polk_mutations$Sample[polk_mutations$PM.Tum > 0.8])),1,12)
  
  if (gene %in% names(methylated)) {
    
    methylation <- methylated[[gene]]
    methylation$V2 <- as.numeric(as.character(methylation$V2))
    methylation$V1 <- as.character(methylation$V1)
    
  }
  
  w <- gene
  
  if (gene == 'SHFM1') w <- 'SEM1'
  
  cn <- copynumber[queryHits(findOverlaps(copynumber, transcr[[ALL_odd[w]]]))]
  
  single_copy <- unique(as.character(cn$Tissue[cn$MinorCN==0 & cn$MajorCN==1]))
  
  both_copies <- unique(as.character(cn$Tissue[cn$MajorCN==0]))
  
  single_copy_and_mutation <- intersect(mutated,single_copy)
  
  if (!is.na(methylation)) {
    
    mutated <- unique(c(mutated, substr(methylation$V1[methylation$V2 < 0.75],1,12)))
    
    single_copy_and_mutation <- unique(c(single_copy_and_mutation,
                                         intersect(single_copy,substr(methylation$V1[methylation$V2 > 0.75],1,12))))
    
    very_mutated <- unique(c(very_mutated, substr(methylation$V1[methylation$V2 > 0.75],1,12)))
    
  }

  samples[[paste0(gene,'_het')]] <- setdiff(c(mutated,single_copy), c(single_copy_and_mutation,very_mutated,both_copies))
  
  samples[[paste0(gene,'_hom')]] <- c(both_copies,single_copy_and_mutation,very_mutated)
  
  if (gene %in% c('POLE','POLD1')) {
    
    samples[[paste0(gene,'_het')]] <- setdiff(mutated, c(single_copy_and_mutation,very_mutated))
    
    samples[[paste0(gene,'_hom')]] <- c(single_copy_and_mutation,very_mutated)
    
  }
  
}
alternative_rownames_bigmat <- paste0('TCGA',substr(rownames(bigmat),5,nchar(row.names(bigmat))))
samples <- lapply(samples, function(l) l[l %in% substr(alternative_rownames_bigmat,1,12)])
one_copy_list <- lapply(one_copy_list, function(l) l[l %in% substr(alternative_rownames_bigmat,1,12)])
sapply(samples,length)


#########################################3

PATH_TO_TCGA_MATRIX='path_to_table_with_mutation_counts_for_TCGA'
bigmat <- read.table(PATH_TO_TCGA_MATRIX,sep='\t')

PATH_TO_METADATA='path_to_table_with_sample_names_and_projects_and_median_normalised_expression'
metadata <- read.table(PATH_TO_METADATA, sep='\t', header=T)

types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
bigmat <- bigmat[rowSums(bigmat)<20000 & rowSums(bigmat) > 50,]
alternative_rownames_bigmat <- paste0('TCGA', substr(rownames(bigmat),5,nchar(rownames(bigmat))))
colnames(bigmat)[1:96] <- types.full

list_of_burden_changes <- list()
ctype <- sapply(alternative_rownames_bigmat, function(x) metadata$project[match(x,metadata$tumour)])

ids <- read.table('ICGCtoTCGA.tsv', sep = '\t', header = T) # clinical info for all the samples - available in ICGC

for (pw in c('MMR','NER','NHEJ','HR','FA','DS','DR','TLS','BER','POLE','DOUBLE')) {
  
  for (can in unique(metadata$project)) {
    
    if (pw == 'MMR') genes <- MMR.core
    if (pw == 'NER') genes <- NER.core
    if (pw == 'BER') genes <- BER.core
    if (pw == 'TLS') genes <- TLS.core
    if (pw == 'HR') genes <- HR.core
    if (pw == 'NHEJ') genes <- NHEJ.core
    if (pw == 'FA') genes <- FA.core
    if (pw == 'DS') genes <- DS.core
    if (pw == 'DR') genes <- DR.core
    if (pw == 'POLE') genes <- 'POLE'
    
    homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% unlist(samples[paste0(genes,'_hom')]) & ctype == can]
    hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% unlist(samples[paste0(genes,'_het')]) & ctype == can]
    nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(genes,'_hom')],samples[paste0(genes,'_het')],
                                                                                          samples['POLE_het']))) & ctype == can]
    
    if (pw == 'MMR') {
      homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples[paste0(MMR.core,'_hom')]),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
      hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples[paste0(MMR.core,'_het')]),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
      nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(MMR.core,'POLE'),'_hom')],
                                                                                            samples[paste0(c(MMR.core,'POLE'),'_het')],
                                                                                            samples['POLE_het']))) & ctype == can]
      
    }
    
    
    if (pw == 'POLE') {
      homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples['POLE_hom']),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
      hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples['POLE_het']),unlist(samples[c(paste0(MMR.core,'_het'),paste0(MMR.core,'_hom'))])) & ctype == can]
      nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(MMR.core,'POLE'),'_hom')],samples[paste0(c(MMR.core,'POLE'),'_het')]))) & ctype == can]
    }
    
    if (pw == 'DOUBLE') {
      homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% intersect(unlist(samples[c('POLD1_hom','POLE_hom')]),unlist(samples[paste0(MMR.core,'_hom')])) & ctype == can]
      hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% intersect(unlist(samples[c('POLD1_het','POLE_het',
                                                                                                          'POLD1_hom','POLE_hom')]),
                                                                                           unlist(samples[c(paste0(MMR.core,'_het'),paste0(MMR.core,'_hom'))])) & ctype == can]
      nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(MMR.core,'POLE'),'_hom')],samples[paste0(c(MMR.core,'POLE'),'_het')]))) & ctype == can]
    }
    
    if (length(hetnames) > 4 && length(nonnames) > 4) {

      if (length(homnames) > 4) {
        if (sum(!is.na(ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(homnames,5,12)), ids$submitted_donor_id)])) > 4)
          list_of_burden_changes[[paste(pw,can,sep='_')]] <- list(
            hom = rowSums(bigmat[homnames,]) / ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(homnames,5,12)), ids$submitted_donor_id)],
            het = rowSums(bigmat[hetnames,]) / ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(hetnames,5,12)), ids$submitted_donor_id)],
            none = rowSums(bigmat[nonnames,]) / ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(nonnames,5,12)), ids$submitted_donor_id)]
          )

      } else {

        if (sum(!is.na(ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(hetnames,5,12)), ids$submitted_donor_id)])) > 4)
          list_of_burden_changes[[paste(pw,can,sep='_')]] <- list(
            het = rowSums(bigmat[hetnames,]) / ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(hetnames,5,12)), ids$submitted_donor_id)],
            none = rowSums(bigmat[nonnames,]) / ids$donor_age_at_diagnosis[match(paste0('TCGA',substr(nonnames,5,12)), ids$submitted_donor_id)]
          )

      }
      
    }
            
  }
  
}


################################################

col_vector <- c(ACC = '#00CD66', BLCA = '#EEAD0E', BRCA = '#CD6090', CESC = '#79CDCD',
                COAD = '#191970', ESCA = '#1E90FF', GBM = '#3D3D3D', HNSC = '#8B2323',
                KICH = '#B32F0B', KIRC = '#FF4500', KIRP = '#FFC100', LAML = '#CD6600', LGG = '#B0B0B0',
                LIHC = '#006400', LUAD = '#FFFFFF', LUSC = '#FDF5E6', MESO = '#698B22',
                OV = '#008B8B', PAAD = '#7A378B', PCPG = '#E066FF', PRAD = '#87CEFA',
                READ = 'blue', SARC = '#FFD700', SKCM = '#000000', STAD = '#BFEFFF',
                TGCT = '#995493', THCA = '#9370DB',
                THYM = '#FFEC8B', UCEC =  '#FF8C69', UCS = '#DDCDCD')


pv <- p.adjust(unlist(lapply(list_of_burden_changes, function(x) { 
  if (length(x) == 3) {
    pvals <- c(wilcox.test(x$hom,x$none,alternative = 'greater')$p.value,wilcox.test(x$het,x$none,alternative = 'greater')$p.value)  
    return(pvals)
  }
  return(wilcox.test(x$het,x$none,alternative = 'greater')$p.value)
})),method = 'BH')

names(pv)[pv < 0.05]

# for interesting pathways only
pdf('Interesting_burden_per_year_per_pathway.pdf', 9,5)
par(mar = c(8,4,2,2))

to.show <- list_of_burden_changes[unique(sapply(names(pv)[pv < 0.05], function(x) ifelse(grepl('1',x,fixed=T) || grepl('2',x,fixed=T), substr(x,1,nchar(x)-1),x)))]
df_burden <- data.frame(unlist(to.show))
colnames(df_burden)[1] <- 'value'
df_burden <- df_burden[!is.na(df_burden$value),,drop = F]
df_burden$pathway <- sapply(rownames(df_burden), function(x) unlist(strsplit(x,split='[_]'))[1])
df_burden$cancer <- sapply(rownames(df_burden), function(x) {
  y <- unlist(strsplit(x,split = '[.]'))[1]
  return(unlist(strsplit(y,split = '[_]'))[2])
})
df_burden$mode <- sapply(rownames(df_burden), function(x) unlist(strsplit(x,split = '[.]'))[2])
df_burden$value <- log10(df_burden$value)
df_burden$name <- sapply(rownames(df_burden), function(x) unlist(strsplit(x,split='[.]'))[1])
df_median <- data.frame(value = log10(unlist(sapply(to.show, function(x) sapply(x,median,na.rm = T)))))
df_median$name <- sapply(rownames(df_median), function(x) unlist(strsplit(x,split = '[.]'))[1])
df_median$mode <- sapply(rownames(df_median), function(x) unlist(strsplit(x,split = '[.]'))[2])
df_median$cancer <- sapply(df_median$name, function(x) unlist(strsplit(x,split = '[_]'))[2])
df_median$pathway <- sapply(df_median$name, function(x) unlist(strsplit(x,split = '[_]'))[1])

plot(NA,NA,xlim = c(0.6,length(unique(df_burden$name))), ylim = c(-2,1.5), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = 'Mutations per year per Mbase')
axis(side = 2, at = c(-2,-1,0,1),
     labels = c(0.01,0.1,1,10), las = 2)
axis(side = 1, col= 'white',at = c(1:length(unique(df_burden$name))), labels='')
     #labels = unique(df_burden$name), srt = 60)
text(seq(1, length(unique(df_burden$name)), by=1), par("usr")[3] - 0.2, labels = unique(df_burden$name), srt = 60, pos = 1, xpd = TRUE)
curx = 1
med_colors <- brewer.pal(n = 3, name = 'Set2')
for (xpoint in unique(df_burden$name)) {
  
  points(rnorm(n = sum(df_burden$name==xpoint & df_burden$mode == 'none'), mean = 0,sd = 0.05) + curx, df_burden$value[df_burden$name == xpoint & df_burden$mode == 'none'], col = 'gray', cex = 0.2)
  points(rnorm(n = sum(df_burden$name==xpoint & df_burden$mode == 'het'), mean = 0,sd = 0.05) + curx + 0.3, df_burden$value[df_burden$name == xpoint & df_burden$mode == 'het'], col = 'gray', cex = 0.2)
  points(rnorm(n = sum(df_burden$name==xpoint & df_burden$mode == 'hom'), mean = 0,sd = 0.05) + curx + 0.6, df_burden$value[df_burden$name == xpoint & df_burden$mode == 'hom'], col = 'gray', cex = 0.2)
  
  lines(x = c(curx,curx+0.3),y = c(df_median$value[df_median$name == xpoint & df_median$mode == 'none'],df_median$value[df_median$name == xpoint & df_median$mode == 'het']))
  lines(x = c(curx-0.05,curx+0.05), y = rep(df_median$value[df_median$name == xpoint & df_median$mode == 'none'],2), col = med_colors[1], lwd = 4)
  lines(x = c(curx+0.25,curx+0.35), y = rep(df_median$value[df_median$name == xpoint & df_median$mode == 'het'],2), col = med_colors[2], lwd = 4)
  if (sum(df_burden$name==xpoint & df_burden$mode == 'hom') > 0) {
    lines(x = c(curx+0.3,curx+0.6),y = c(df_median$value[df_median$name == xpoint & df_median$mode == 'het'],df_median$value[df_median$name == xpoint & df_median$mode == 'hom']))
    lines(x = c(curx+0.55,curx+0.65), y = rep(df_median$value[df_median$name == xpoint & df_median$mode == 'hom'],2), col = med_colors[3], lwd = 4)
    lines(x = c(curx+0.25,curx+0.35), y = rep(df_median$value[df_median$name == xpoint & df_median$mode == 'het'],2), col = med_colors[2], lwd = 4)
  }
  
  curx = curx + 1
}
legend('topleft',legend = c('No mutations','Heterozygous','Homozygous'), fill = med_colors, bty = 'n', cex = 1.2)


dev.off()
###########################################

types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
bigmat <- read.table('TCGA.caveman.matrix.dat')
bigmat <- bigmat[rowSums(bigmat)<20000 & rowSums(bigmat) > 50,]
alternative_rownames_bigmat <- paste0('TCGA', substr(rownames(bigmat),5,nchar(rownames(bigmat))))
colnames(bigmat)[1:96] <- types.full
ctype <- sapply(alternative_rownames_bigmat, function(x) metadata$project[match(x,metadata$tumour)])


profdf <- data.frame()
for (pw in c('MMR','NER','NHEJ','HR','FA','DS','DR','TLS','BER','POLE','DOUBLE')) {
  
  if (pw == 'MMR') genes <- MMR.core
  if (pw == 'NER') genes <- NER.core
  if (pw == 'BER') genes <- BER.core
  if (pw == 'TLS') genes <- TLS.core
  if (pw == 'HR') genes <- HR.core
  if (pw == 'NHEJ') genes <- NHEJ.core
  if (pw == 'FA') genes <- FA.core
  if (pw == 'DS') genes <- DS.core
  if (pw == 'DR') genes <- DR.core
  if (pw == 'POLE') genes <- 'POLE'
  if (pw == 'DOUBLE') genes <- c('POLE','POLD1',MMR.core)
  
  for (can in unique(metadata$project)) {
    
    homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples[paste0(genes,'_hom')]),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
    hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples[paste0(genes,'_het')]),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
    nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(genes,'POLE'),'_hom')],samples[paste0(c(genes,'POLE'),'_het')]))) & ctype == can]
    
    
    if (pw == 'MMR') {
      homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples[paste0(MMR.core,'_hom')]),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
      hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples[paste0(MMR.core,'_het')]),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
      nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(MMR.core,'POLE'),'_hom')],samples[paste0(c(MMR.core,'POLE'),'_het')]))) & ctype == can]
      
    }
    
    
    if (pw == 'POLE') {
      homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples['POLE_hom']),unlist(samples[c('POLE_hom','POLE_het')])) & ctype == can]
      hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% setdiff(unlist(samples['POLE_het']),unlist(samples[c(paste0(MMR.core,'_het'),paste0(MMR.core,'_hom'))])) & ctype == can]
      nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(MMR.core,'POLE'),'_hom')],samples[paste0(c(MMR.core,'POLE'),'_het')]))) & ctype == can]
    }
    
    if (pw == 'DOUBLE') {
      homnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% intersect(unlist(samples[c('POLD1_hom','POLE_hom')]),unlist(samples[paste0(MMR.core,'_hom')])) & ctype == can]
      hetnames <- rownames(bigmat)[substr(alternative_rownames_bigmat,1,12) %in% intersect(unlist(samples[c('POLD1_het','POLE_het',
                                                                                                            'POLD1_hom','POLE_hom')]),
                                                                                           unlist(samples[c(paste0(MMR.core,'_het'),paste0(MMR.core,'_hom'))])) & ctype == can]
      nonnames <- rownames(bigmat)[!(substr(alternative_rownames_bigmat,1,12) %in% unlist(c(samples[paste0(c(MMR.core,'POLE'),'_hom')],samples[paste0(c(MMR.core,'POLE'),'_het')]))) & ctype == can]
    }
    
    
    # regress out the apobec and signatures 1 and 5
    
    newmat <- bigmat[c(homnames,hetnames,nonnames),]
    
    decomp <- nmFit(D = newmat[,1:96], E = t(as.matrix(new_cancer_signatures[,-c(1:2)])))
    newmat[,1:96] <- newmat[,1:96] - decomp[,c('SBS1','SBS5','SBS2','SBS13','SBS30')] %*% t(new_cancer_signatures[,c('SBS1','SBS5','SBS2','SBS13','SBS30')])
    
    newmat[newmat<0] <- 0
    
    homnames <- homnames[rowSums(newmat[homnames,])>100 & rowSums(newmat[homnames,])<20000]
    hetnames <- hetnames[rowSums(newmat[hetnames,])>100 & rowSums(newmat[hetnames,])<20000]
    nonnames <- nonnames[rowSums(newmat[nonnames,])>100 & rowSums(newmat[nonnames,])<20000]
    
    if (length(hetnames) > 3 && length(nonnames) > 3) {
      none_profile <- apply(newmat[nonnames,] / rowSums(newmat[nonnames,]),2,median)
      het_profile <- apply(newmat[hetnames,] / rowSums(newmat[hetnames,]),2,median)
      if (length(homnames) > 3) {
        hom_profile <- apply(newmat[homnames,] / rowSums(newmat[homnames,]),2,median)
        profdf <- rbind(profdf, data.frame(
          cancer = can,
          pathway = pw,
          mode = c('hom','het'),
          distance = c(1 - cosine(hom_profile/sum(hom_profile),none_profile/sum(none_profile)),
                       1 - cosine(het_profile/sum(het_profile),none_profile/sum(none_profile)))
        ))
      } else {
        profdf <- rbind(profdf, data.frame(
          cancer = can,
          pathway = pw,
          mode = 'het',
          distance = 1 - cosine(het_profile/sum(het_profile),
                                none_profile/sum(none_profile))
        ))
      }
      
    }
    
  }
  
}

df_median <- data.frame(value = log10(unlist(sapply(list_of_burden_changes, function(x) sapply(x,median,na.rm = T)))))
df_median$name <- sapply(rownames(df_median), function(x) unlist(strsplit(x,split = '[.]'))[1])
df_median$mode <- sapply(rownames(df_median), function(x) unlist(strsplit(x,split = '[.]'))[2])
df_median$cancer <- sapply(df_median$name, function(x) unlist(strsplit(x,split = '[_]'))[2])
df_median$pathway <- sapply(df_median$name, function(x) unlist(strsplit(x,split = '[_]'))[1])

plot(x = df_median$value, y = profdf$distance, bty = 'n', xlab = 'Median burden folchange', ylab = 'Average distance')

profdf$pathway <- as.character(profdf$pathway)
profdf$pathway[profdf$pathway == 'POLE+MMR'] <- 'POLE+\nMMR'
profdf$pathway <- factor(profdf$pathway, levels = c('BER','DR','DS','FA','HR','MMR','POLE','POLE+\nMMR','NER','NHEJ','TLS'))
write.table(profdf,file='~/Distances_for_profiles_between_altered_and_wildtype_DNA_repair.dat')
pdf('Profile_change_per_pathway.pdf', 9,5.5)
library(beeswarm)
beeswarm(distance ~ pathway, method = 'center', corralWidth = 0.8, corral = 'wrap',
         data = profdf,
         pwpch = ifelse(profdf$mode == 'het',21,23),
         pwbg = col_vector[as.character(profdf$cancer)],
         lty = 0,
         xlab = '', ylab = 'Cosine distance', main = 'Change of the average profile')
axis(side = 2, labels = NA)
legend('topleft', legend = c('Heterozygous','Homozygous'), pch = c(21,23),pt.bg = 'black', bty = 'n', cex = 1.3)
dev.off()




