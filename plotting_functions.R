# Function for plotting 96-long signatures / profiles with or without CIs
plot_sig_wb <- function (mut_matrix, colors=NA, norm=T, ymax = NA, ymin = NA, flip=F, ind=F, font=11, size=12,
                         CI = F, low = NULL, high = NULL, DNV = FALSE, diff_scale = F,
                         rownames = T, diff_limits = NULL, thr = 0.05, size.all = 10, err.bar.size = 0.5)
{
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  IND_TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  ylab = 'Number of mutations'
  if (norm) {
    mult <- colSums(mut_matrix)
    mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ylab = 'Relative contribution'
  }
  if (is.na(ymax)) ymax = max(mut_matrix)
  if (is.na(ymin)) ymin = 0
  
  
  if (diff_scale)
    scales <- 'free'
  else scales <- 'free_x'
  
  TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  TRIPLETS_112 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3),IND_TRIPLETS)
  if (norm) norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x)) else norm_mut_matrix <- mut_matrix
  if (missing(colors)) {
    if (!ind)  colors = c("#2EBAED", "#000000", "#DE1C14",
                          "#D4D2D2", "#ADCC54", "#F0D0CE")
    else colors = c("#2EBAED", "#000000", "#DE1C14", "yellow",
                    "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  context = TRIPLETS_96
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  if (ind) {
    context = TRIPLETS_112
    substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G","INDEL"), each = 16)
  }
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  if (flip==T) norm_mut_matrix = -norm_mut_matrix
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (ymax>0) breaks = c(0,0.1) else breaks = c(-0.1,0)
  if (ymax>0) lbls = c('0','0.1') else lbls = c('0.1','0')

  
  if (!is.null(diff_limits)) {
    df3 <- rbind(df3, data.frame(substitution = 'Zfake', context = '0', variable = unique(df3$variable), value = diff_limits))
  }
  
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution)) +
    geom_bar(stat = "identity", colour = "black", size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = scales) +
    ylab(ylab) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = size.all,vjust = 1), 
          axis.text.y = element_text(size = size.all), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = size.all/2, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size.all), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) plot <- plot + coord_cartesian(ylim = c(0,ymax), clip='off') 
  if (!rownames) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  

  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df3$value_min = df_lower$value
    df3$value_max = df_upper$value
    
    if (!is.null(diff_limits)) {
      df3$value_min = c(df_lower$value, rep(0,ncol(mut_matrix)))
      df3$value_max = c(df_upper$value, rep(0,ncol(mut_matrix)))
    } else {
      df3$value_min = df_lower$value
      df3$value_max = df_upper$value
    }
    
    plot <- plot + geom_linerange(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', size = err.bar.size)
  }
  
  return(plot)
  
  
}

# plot_id_only <- function(mut_matrix, colors = NA, norm = T, ymax = NA, ytitle = 'Number of mutations') {
#   
#   ID <- c('DEL_C_1','DEL_T_1','INS_C_1','INS_T_1','DEL_repeats_2','DEL_repeats_3','DEL_repeats_4','DEL_repeats_5+',
#           'INS_repeats_2','INS_repeats_3','INS_repeats_4','INS_repeats_5+','DEL_MH_2','DEL_MH_3','DEL_MH_4','DEL_MH_5+')
#   
#   id.types.full <- c('DEL_C_1_0','DEL_C_1_1','DEL_C_1_2','DEL_C_1_3','DEL_C_1_4','DEL_C_1_5+','DEL_T_1_0','DEL_T_1_1','DEL_T_1_2',
#                      'DEL_T_1_3','DEL_T_1_4','DEL_T_1_5+','INS_C_1_0','INS_C_1_1','INS_C_1_2','INS_C_1_3','INS_C_1_4','INS_C_1_5+',
#                      'INS_T_1_0','INS_T_1_1','INS_T_1_2','INS_T_1_3','INS_T_1_4','INS_T_1_5+','DEL_repeats_2_0','DEL_repeats_2_1',
#                      'DEL_repeats_2_2','DEL_repeats_2_3','DEL_repeats_2_4','DEL_repeats_2_5+','DEL_repeats_3_0','DEL_repeats_3_1',
#                      'DEL_repeats_3_2','DEL_repeats_3_3','DEL_repeats_3_4','DEL_repeats_3_5+','DEL_repeats_4_0','DEL_repeats_4_1',
#                      'DEL_repeats_4_2','DEL_repeats_4_3','DEL_repeats_4_4','DEL_repeats_4_5+','DEL_repeats_5+_0','DEL_repeats_5+_1',
#                      'DEL_repeats_5+_2','DEL_repeats_5+_3','DEL_repeats_5+_4','DEL_repeats_5+_5+','INS_repeats_2_0','INS_repeats_2_1',
#                      'INS_repeats_2_2','INS_repeats_2_3','INS_repeats_2_4','INS_repeats_2_5+','INS_repeats_3_0','INS_repeats_3_1',
#                      'INS_repeats_3_2','INS_repeats_3_3','INS_repeats_3_4','INS_repeats_3_5+','INS_repeats_4_0','INS_repeats_4_1',
#                      'INS_repeats_4_2','INS_repeats_4_3','INS_repeats_4_4','INS_repeats_4_5+','INS_repeats_5+_0','INS_repeats_5+_1',
#                      'INS_repeats_5+_2','INS_repeats_5+_3','INS_repeats_5+_4','INS_repeats_5+_5+','DEL_MH_2_1','DEL_MH_3_1','DEL_MH_3_2',
#                      'DEL_MH_4_1','DEL_MH_4_2','DEL_MH_4_3','DEL_MH_5+_1','DEL_MH_5+_2','DEL_MH_5+_3','DEL_MH_5+_4','DEL_MH_5+_5+')
#   
#   if (is.na(colors)) {
#     library(RColorBrewer)
#     colors <- c( brewer.pal(10,'Set3'), 'orange')
#   }
#   
#   if (norm) {
#     norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
#     ytitle <- 'Relative contribution'
#   } else norm_mut_matrix <- mut_matrix
#   
#   if (is.na(ymax)) ymax <- max(norm_mut_matrix)
#   
#   context = dnv.types.full
#   substitution = rep(DNV, c(9,6,9,6,9,6,6,9,9,9))
#   df = data.frame(substitution = substitution, context = context)
#   df2 = cbind(df, as.data.frame(norm_mut_matrix))
#   df3 = melt(df2, id.vars = c("substitution", "context"))
#   df3$substitution = factor(df3$substitution, levels=DNV)
#   plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
#     geom_bar(stat = "identity", colour = "black",size = 0.2) +
#     scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = 'free_x') +
#     ylab(ytitle) + coord_cartesian(ylim = c(0,ymax)) + 
#     guides(fill = FALSE) + theme_bw() +
#     theme(text = element_text(family='ArialMT'),
#           axis.title.y = element_text(size = 12,vjust = 1), 
#           axis.text.y = element_text(size = 8), 
#           axis.title.x = element_blank(), 
#           axis.text.x = element_text(size = 6, angle = 90, vjust = 0.4), 
#           strip.text.x = element_text(size = 12), 
#           strip.text.y = element_text(size = 12), 
#           panel.grid.major.x = element_blank(), 
#           strip.background = element_blank(), 
#           panel.border = element_rect(colour="white"),
#           panel.spacing = unit(0.1,'lines'))
#   return(plot)
#   
#   
#   
# }


# Function for saving a list of ggplots as pages of a PDF file
GG_save_pdf = function(list, filename, ...) {
  #start pdf
  pdf(filename, ...)
  #loop
  count <- 1
  for (p in list) {
    print(p)
    print(count)
    count = count +1
  }
  #end pdf
  dev.off()
  invisible(NULL)
}




dnv.types.full <- c(paste0(rep('AC',9),'>',unlist(lapply(c('C','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                    paste0(rep('AT',6),'>',c('CA','CC','CG','GA','GC','TA')),
                    paste0(rep('CC',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                    paste0(rep('CG',6),'>',c('AT','GC','GT','TA','TC','TT')),
                    paste0(rep('CT',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','C','G'))))),
                    paste0(rep('GC',6),'>',c('AA','AG','AT','CA','CG','TA')),
                    paste0(rep('TA',6),'>',c('AT','CG','CT','GC','GG','GT')),
                    paste0(rep('TC',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                    paste0(rep('TG',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','T'))))),
                    paste0(rep('TT',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','G'))))))

# Plotting dinucleotide signatures/profiles with or without CI
plot_dnv_only <- function(mut_matrix, colors=NA, norm=T, ymax=NA, ytitle = 'Number of mutations',
                          CI = F, low = NA, high = NA, rownames = T, diff_scale = F,diff_limits = NULL, 
                          size.all = 10, err.bar.size = 0.5) {
  
  DNV <- c('AC','AT','CC','CG','CT','GC','TA','TC','TG','TT')
  
  dnv.types.full <- c(paste0(rep('AC',9),'>',unlist(lapply(c('C','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                      paste0(rep('AT',6),'>',c('CA','CC','CG','GA','GC','TA')),
                      paste0(rep('CC',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                      paste0(rep('CG',6),'>',c('AT','GC','GT','TA','TC','TT')),
                      paste0(rep('CT',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','C','G'))))),
                      paste0(rep('GC',6),'>',c('AA','AG','AT','CA','CG','TA')),
                      paste0(rep('TA',6),'>',c('AT','CG','CT','GC','GG','GT')),
                      paste0(rep('TC',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                      paste0(rep('TG',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','T'))))),
                      paste0(rep('TT',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','G'))))))
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c( brewer.pal(10,'Set3'), 'orange')
  }
  
  if (norm) {
    mult <- colSums(mut_matrix)
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ytitle <- 'Relative contribution'
  } else norm_mut_matrix <- mut_matrix
  
  if (is.na(ymax)) ymax <- max(norm_mut_matrix)
  
  scales <- 'free_x'
  if (diff_scale)
    scales <- 'free'
  
  context = dnv.types.full
  substitution = rep(DNV, c(9,6,9,6,9,6,6,9,9,9))
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$substitution = factor(df3$substitution, levels=DNV)
  
  
  if (!is.null(diff_limits)) {
    df3 <- rbind(df3, data.frame(substitution = 'Zfake', context = '0', variable = unique(df3$variable), value = diff_limits))
  }
  
  
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = scales) +
    ylab(ytitle) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = size.all,vjust = 1), 
          axis.text.y = element_text(size = size.all), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = size.all/2, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size.all), 
          strip.text.y = element_text(size = size.all), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) plot <- plot + coord_cartesian(ylim = c(0,ymax), clip='off')
  if (!rownames) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    if (!is.null(diff_limits)) {
      df3$value_min = c(df_lower$value, rep(0,ncol(mut_matrix)))
      df3$value_max = c(df_upper$value, rep(0,ncol(mut_matrix)))
    } else {
      df3$value_min = df_lower$value
      df3$value_max = df_upper$value
    }
    plot <- plot + geom_linerange(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', size = err.bar.size)
  }
  
  return(plot)
  
}

# Plotting base substitution and indel profiles/signatures with/without CIs
plot_subindel_wb <- function(mut_matrix, colors=NA, norm=T, ymax = NA, ymin = NA, 
                             CI = F, low = NULL, high = NULL, DNV = FALSE, diff_scale = F,
                             rownames = T, diff_limits = NULL, thr = 0.05, size.all = 10, err.bar.size = 0.5)
{
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE","orange",'darkmagenta',"purple",'brown')
  }
  
  TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  INDELS = c(
    '1-4 bp', "5-49 bp","50-400 bp",
    '1-49 bp','50-400 bp',
    '1-4 bp', '5-49 bp','50-400 bp')
  
  ylab = 'Number of mutations'
  if (norm) {
    mult <- colSums(mut_matrix)
    mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ylab = 'Relative contribution'
  }
  
  if (diff_scale)
    scales <- 'free'
  else scales <- 'free_x'
  
  
  if (is.na(ymax)) ymax = max(mut_matrix)
  if (is.na(ymin)) ymin = 0
  
  TRIPLETS_96 = rep(TRIPLETS, 6)
  context = c(TRIPLETS_96,INDELS)
  if (DNV) context <- c(context,'DNV')
  substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),
                   rep('Deletions',3),rep('Complex\n indels',2),rep('Insertions',3))
  if (DNV) substitution <- c(substitution,'DNV')
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  if (!DNV) df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",
                                                       "Deletions","Complex\n indels","Insertions"))
  if (DNV) df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",
                                                               "Deletions","Complex\n indels","Insertions",'DNV'))
  df3$width <- 0.6
  df3$width[df3$substitution=='Deletions' | df3$substitution=='Insertions'] <- 0.3;
  if (DNV) df3$width[df3$substitution == 'DNV'] <- 0.1
  df3$width[df3$substitution=='Complex\n indels'] <- 0.2
  
  if (!is.null(diff_limits)) {
    df3 <- rbind(df3, data.frame(substitution = 'Zfake', context = '0', variable = unique(df3$variable), value = diff_limits, width = 0.2))
  }
  
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "black", size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = scales) +
    ylab(ylab) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = size.all,vjust = 1), 
          axis.text.y = element_text(size = size.all), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = size.all/2, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size.all), 
          strip.text.y = element_text(size = size.all), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) plot <- plot + coord_cartesian(ylim = c(0,ymax), clip='off') 
  if (!rownames) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df3$value_min = df_lower$value
    df3$value_max = df_upper$value
    
    if (!is.null(diff_limits)) {
      df3$value_min = c(df_lower$value, rep(0,ncol(mut_matrix)))
      df3$value_max = c(df_upper$value, rep(0,ncol(mut_matrix)))
    } else {
      df3$value_min = df_lower$value
      df3$value_max = df_upper$value
    }
    
    plot <- plot + geom_linerange(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', size = err.bar.size)
  }
  
  return(plot)
  
}


############################################################################################################

# visualisation for the huge profile! (worms)

TRIPLETS = c(
  "A*A", "A*C", "A*G", "A*T",
  "C*A", "C*C", "C*G", "C*T",
  "G*A", "G*C", "G*G", "G*T",
  "T*A", "T*C", "T*G", "T*T")
SVS = c("Tandem dup.","Deletions","Inversions","Complex events","Translocations","Interchr. events","Foldbacks")
INDELS = c(
  '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
  "6-50 bp","over 50 bp",
  '1-49 bp','50-400 bp',
  '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
  '6-50 bp','over 50 bp')
library(grid)
library(ggplot2)
library(reshape2)

# Plotting full 119-long signature/profile
plot_fullsig_wb <- function(mut_matrix, colors=NA, norm=T, ymax=NA, ytitle = 'Number of mutations', 
                               CI = F, low = NA, high = NA, diff_scale = F, rownames = T, diff_limits = NULL,
                               significance = NULL, thr = 0.05, size.all = 10, err.bar.size = 0.5, yspace = NULL)
{
  types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[',
                      rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']',
                      rep(c('A','C','G','T'), 24), sep='')
  
  ID <- c('Deletions','Complex\n indels','Insertions')
  
  if (is.na(colors[1])) {
    library(RColorBrewer)
    colors <- c(brewer.pal(6, "Set1"), 'brown', brewer.pal(3,"Set2"), 'darkmagenta')
  }
  
  TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  INDELS = c(
    '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
    "6-50 bp","over 50 bp",
    '1-49 bp','50-400 bp',
    '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
    '6-50 bp','over 50 bp')
  
  SVS = c("Tandem dup.","Deletions","Inversions","Complex events","Translocations","Interchr. events","Foldbacks")
  
  if (norm) {
    mult <- colSums(mut_matrix)
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ytitle <- 'Relative contribution'
  } else norm_mut_matrix <- mut_matrix
  
  
  if (is.na(ymax)) ymax <- max(norm_mut_matrix)
  if(is.null(yspace)) yspace = 'fixed'
  
  TRIPLETS_96 = rep(TRIPLETS, 6)
  context = c(TRIPLETS_96,'2 bp','over 2 bp',INDELS,SVS)
  substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),
                   'MNV',
                   'MNV',
                   rep(ID,c(6,2,6)),
                   rep('Structural\n Variants',7))
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$context <- as.character(df3$context)
  df3$context <- factor(df3$context, levels = unique(context))
  if (!is.null(significance)) {
    significance <- significance[,colnames(mut_matrix)]
    dfsig <- melt(significance)
    df3$significance <- dfsig$value
    df3$significance[df3$significance >= thr] <- 0.4
    df3$significance[df3$significance < thr] <- 1
  } else
    df3$significance <- 1
  df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",'MNV',ID,"Structural\n Variants"))
  df3$width = rep(0.7,nrow(df3));
  df3$width[df3$substitution == 'MNV' | df3$substitution == 'Complex\n indels'] <- 0.2
  if (diff_scale)
    scales <- 'free'
  else scales <- 'free_x'
  
  if (!is.null(diff_limits)) {
    df3 <- rbind(df3, data.frame(substitution = 'Zfake', context = '0', variable = unique(df3$variable), value = diff_limits, significance = 0.4, width = 0.2))
  }
  
  p = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "black", size = NA, alpha = df3$significance[order(df3$variable)]) +
    scale_fill_manual(values = c(colors,'white')) + facet_grid(variable ~ substitution, scales = scales, space = yspace) +
    ylab(ytitle) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = size.all,vjust = 1), 
          axis.text.y = element_text(size = size.all), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = size.all/2, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size.all), 
          strip.text.y = element_text(size = size.all),
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) p <- p + coord_cartesian(ylim = c(0,ymax), clip='off') 
  if (!rownames) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    if (!is.null(diff_limits)) {
      df3$value_min = c(df_lower$value, rep(0,ncol(mut_matrix)))
      df3$value_max = c(df_upper$value, rep(0,ncol(mut_matrix)))
    } else {
      df3$value_min = df_lower$value
      df3$value_max = df_upper$value
    }
    p <- p + geom_linerange(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', size = err.bar.size)
#    plot <- plot + geom_pointrange(data=df3, aes(ymin=value, ymax=value_max, colour=substitution),
#                                  fatten = 0.001, size=0.5, show.legend = F) +
#      scale_color_manual(values=c(colors,"white")) +
#      geom_pointrange(data=df3, aes(ymin=value_min,ymax=value),
#                      size=0.5, fatten = 0.001,show.legend = F, colour="white")
  }

  return(p)
  
}

intersperse_119 <- function(y) {
  return(c(y[1:16],NA,
           y[17:32],NA,
           y[33:48],NA,
           y[49:64],NA,
           y[65:80],NA,
           y[81:96],NA,
           y[97],NA,
           y[98],NA,
           y[99:104],NA,
           y[105],NA,
           y[106],NA,
           y[107:112],NA,
           y[113:119]))
}

intersperse_104 <- function(y) {
  return(c(y[1:16],NA,
           y[17:32],NA,
           y[33:48],NA,
           y[49:64],NA,
           y[65:80],NA,
           y[81:96],NA,
           y[97],NA,
           y[98],NA,
           y[99],NA,
           y[100],NA,
           y[101],NA,
           y[102],NA,
           y[103],NA,
           y[104]))
}


##################

# Plotting interaction effects for worms
interaction_effect_plot <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high, cex = 3,  lwd = 4, labels=NA, lwd.means = 4, cex.main = 1,
                           col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta"),
                           means =NULL, means_low = NULL, means_high = NULL) {
  if (length(y) != 119) return('Working progress!')
  fff <- data.frame(bar = intersperse_119(y))
  if (log) fff$bar <- fff$bar / log(10)
  line_X_axis <- cumsum(c(rep(1,101),4,2,6,2,4,rep(2.5,6),5,2,5,2,5,rep(2.5,6),2,rep(2,7)))
  plot(x = line_X_axis, y = rep(NA,length(line_X_axis)),
       xlim = c(0,max(line_X_axis)),
       ylim = c(min(c(at,min(fff$bar,na.rm=T)),na.rm=T),max(c(at,max(fff$bar,na.rm=T)),na.rm=T)),
       yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty='n')
  
  fff[fff < -1] <- -1
  
  if (CI) {
    fff$low <- intersperse_119(low)
    fff$high <- intersperse_119(high)
    if (log) {
      fff$low <- fff$low / log(10)
      fff$high <- fff$high / log(10)
    }
    
    fff[fff < -1] <- -1
    
    col1 <- sapply(col, function(x) {
      tmpf <- colorRampPalette(c(x, "white"))
      return(tmpf(5)[4])
    })
    
    pos = 1
    for (j in 1:6) {
      for (k in pos:(pos+15))
        polygon(x = c(line_X_axis[k] - 1, rep(line_X_axis[k], 2),
                      rep(line_X_axis[k] - 1,2)),
                y = c(rep(fff$high[k],2),
                      rep(fff$low[k],2), fff$high[k]),
                col = col1[j], border = NA)
      pos = pos + 17
    }
    
    for (k in c(103,105))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[7], border = NA)
      
    for (k in c(114,116))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[9], border = NA)
    
    for (k in c(107:112))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[8], border = NA)
    
    for (k in c(118:123))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[10], border = NA)

    for (k in c(125:131))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[11], border = NA)

  }
  # All means
  col1 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "white"))
    return(tmpf(5)[2])
  })
  pos = 1
  for (j in 1:6) {
    for (k in pos:(pos+15))
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col1[j], lwd = lwd.means/2)
    pos = pos + 17
  }
  for (k in c(103,114))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[ifelse(k==103,7,9)], lwd = lwd.means/2)
  for (k in c(105,116))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[ifelse(k==105,7,9)], lwd = lwd.means/2)
  for (k in c(107:112,118:123))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[ifelse(k<117,8,10)], lwd = lwd.means/2)
  for (k in c(125:131))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[11], lwd = lwd.means/2)
  
  # Significant means
  #col2 <- sapply(col, function(x) {
  #  tmpf <- colorRampPalette(c(x, "black"))
  #  return(tmpf(5)[1])
  #})
  col2 <- col
  pos = 1
  for (j in 1:6) {
    pos1 <- c(pos:(pos+15))[which(fff$high[pos:(pos+15)] < 0 | fff$low[pos:(pos+15)] > 0)]
    for (k in pos1)
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col2[j], lwd = lwd.means)
    pos = pos + 17
  }
  pos1 <- which(fff$high < 0 | fff$low > 0)
  pos1 <- pos1[pos1>101]
  
  for (k in pos1) {
    if (k %in% c(103,105))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[7], lwd = lwd.means)
    if (k %in% c(107:112))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[8], lwd = lwd.means)
    else if (k %in% c(114,116))
      lines(x = c(line_X_axis[k] -2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[9], lwd = lwd.means)
    else if (k %in% c(118:123))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[10], lwd = lwd.means)
    else 
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[11], lwd = lwd.means)
  }
  
  
  # Class means
  if(is.null(means)) {
    pos = 1
    for (j in 1:6) {
     lines(x = c(line_X_axis[pos]-1,line_X_axis[pos+15]), y = rep(mean(fff$bar[pos:(15+pos)]),2))
      pos = pos + 17
    }
    lines(x = c(line_X_axis[103]-2,line_X_axis[105]), y = rep(mean(fff$bar[c(103,105)]),2))
    lines(x = c(line_X_axis[107]-2,line_X_axis[112]), y = rep(mean(fff$bar[107:112]),2))
    lines(x = c(line_X_axis[114]-2,line_X_axis[116]), y = rep(mean(fff$bar[c(114,116)]),2))
    lines(x = c(line_X_axis[118]-2,line_X_axis[123]), y = rep(mean(fff$bar[118:123]),2))
    lines(x = c(line_X_axis[125]-2,line_X_axis[131]), y = rep(mean(fff$bar[125:131]),2))
    
  } else {
    if (log) {
      means <- means / log(10)
      means_low <- means_low / log(10)
      means_high <- means_high / log(10)
    }
    pos = 1
    for (j in 1:6) {
      polygon(x = c(line_X_axis[pos] - 1, rep(line_X_axis[pos+15], 2),
                    rep(line_X_axis[pos] - 1,2)),
              y = c(rep(means_high[j],2),
                    rep(means_low[j],2), means_high[j]),
              col = adjustcolor('gray86', alpha.f = 0.4), border = NA)
      lines(x = c(line_X_axis[pos]-1,line_X_axis[pos+15]), y = rep(means[j],2))
      pos = pos + 17
    }
    
    polygon(x = c(line_X_axis[103]-2, rep(line_X_axis[105], 2),
                  rep(line_X_axis[103]-2,2)),
            y = c(rep(means_high[7],2),
                  rep(means_low[7],2), means_high[7]),
            col = adjustcolor('gray86', alpha.f = 0.4), border = NA)
    lines(x = c(line_X_axis[103]-2,line_X_axis[105]), y = rep(means[7],2))
    
    polygon(x = c(line_X_axis[107]-2, rep(line_X_axis[112], 2),
                  rep(line_X_axis[107]-2,2)),
            y = c(rep(means_high[8],2),
                  rep(means_low[8],2), means_high[8]),
            col = adjustcolor('gray86', alpha.f = 0.4), border = NA)
    lines(x = c(line_X_axis[107]-2,line_X_axis[112]), y = rep(means[8],2))
    
    polygon(x = c(line_X_axis[114]-2, rep(line_X_axis[116], 2),
                  rep(line_X_axis[114]-2,2)),
            y = c(rep(means_high[9],2),
                  rep(means_low[9],2), means_high[9]),
            col = adjustcolor('gray86', alpha.f = 0.4), border = NA)
    lines(x = c(line_X_axis[114]-2,line_X_axis[116]), y = rep(means[9],2))
    
    polygon(x = c(line_X_axis[118]-2, rep(line_X_axis[123], 2),
                  rep(line_X_axis[118]-2,2)),
            y = c(rep(means_high[10],2),
                  rep(means_low[10],2), means_high[10]),
            col = adjustcolor('gray86', alpha.f = 0.4), border = NA)
    lines(x = c(line_X_axis[118]-2,line_X_axis[123]), y = rep(means[10],2))
    
    polygon(x = c(line_X_axis[125]-2, rep(line_X_axis[131], 2),
                  rep(line_X_axis[125] - 2,2)),
            y = c(rep(means_high[11],2),
                  rep(means_low[11],2), means_high[11]),
            col = adjustcolor('gray86', alpha.f = 0.4), border = NA)
    lines(x = c(line_X_axis[125]-2,line_X_axis[131]), y = rep(means[11],2))
    
  }
  
  
  title(main=plot_main, cex.main = cex.main)
  if (is.na(labels)) {
    if (log) labels = 10**at else labels = at
    axis(side = 2, at = c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) 
      c(2:9) * 10**i)))), labels = c(rep('',8*(max(at) - min(at)))), las = 2, pos = -0.5, tck = -0.02)
    axis(side = 2, at = c(at), labels = c(as.character(labels)), las = 2, pos = -0.5, tck = -0.03)
    #axis(side = 2, col = 'white', tcl = 0, labels = NA, lwd = 3, pos = -0.5)
    if (log) lines(x = c(-0.5,185.5), y = c(0,0),lty = 2, col = 'red')
    else  lines(x = c(-0.5,185.5), y = c(1,1),lty = 2, col = 'red')
  }  else {
    axis(side = 2, at = c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) 
      c(2:9) * 10**i)))), labels = c(rep('',8*(max(at) - min(at)))), las = 2, pos = -0.5, tck = -0.02)
    axis(side = 2, at = c(at), labels = c(as.character(labels)), las = 2, pos = -0.5, tck = -0.03)
    #axis(side = 2, col = 'white', lwd = 3, tcl = 0, labels = NA, pos = -0.5)
    if (log) lines(x = c(-0.5,185.5), y = c(0,0),lty = 2, col = 'red')
    else  lines(x = c(-0.5,185.5), y = c(1,1),lty = 2, col = 'red')
  }
  
}

# Plotting interaction effects in humans
interaction_effect_plot_human <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high, cex = 3, lwd = 4, labels=NA, lwd.means = 4,
                                    col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","orange",'darkmagenta',"purple"),
                                    cex.main = 1) {

  fff <- data.frame(bar = intersperse_104(y))
  if (log) fff$bar <- fff$bar / log(10)
  line_X_axis <- cumsum(c(rep(1,101),3,2,3,2,3,2,6,2,5,2,6,2,3,2,3,2)) - 1
  plot(x = line_X_axis, y = rep(NA,length(line_X_axis)),
       xlim = c(0,max(line_X_axis)),
       ylim = c(min(at),max(at)),
       yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty='n')
  fff[fff < -1] <- -1
  if (CI) {
    fff$low <- intersperse_104(low)
    fff$high <- intersperse_104(high)
    if (log) {
      fff$low <- fff$low / log(10)
      fff$high <- fff$high / log(10)
    }
    fff[fff < -1] <- -1
    col1 <- sapply(col, function(x) {
      tmpf <- colorRampPalette(c(x, "white"))
      return(tmpf(7)[6])
    })
    pos = 1
    for (j in 1:6) {
      for (k in pos:(pos+15))
        polygon(x = c(line_X_axis[k] - 1, rep(line_X_axis[k], 2),
                      rep(line_X_axis[k] - 1,2)),
                y = c(rep(fff$high[k],2),
                      rep(fff$low[k],2), fff$high[k]),
                col = col1[j], border = NA)
      pos = pos + 17
    }
    for (k in c(103,105,107))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[7], border = NA)
    
    for (k in c(109,111))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[8], border = NA)
    
    for (k in c(113,115,117))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[9], border = NA)
    
  }
  
  col1 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "white"))
    return(tmpf(7)[4])
  })
  pos = 1
  for (j in 1:6) {
    for (k in pos:(pos+15))
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col1[j], lwd = lwd.means/2)
    pos = pos + 17
  }
  for (k in c(103,105,107))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[7], lwd = lwd.means/2)
  for (k in c(109,111))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[8], lwd = lwd.means/2)
  for (k in c(113,115,117))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[9], lwd = lwd.means/2)
  
  col2 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "black"))
    return(tmpf(5)[2])
  })
  pos = 1
  for (j in 1:6) {
    pos1 <- c(pos:(pos+15))[which(fff$high[pos:(pos+15)] < 0 | fff$low[pos:(pos+15)] > 0)]
    for (k in pos1)
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col2[j], lwd = lwd.means)
    pos = pos + 17
  }
  pos1 <- which(fff$high < 0 | fff$low > 0)
  pos1 <- pos1[pos1>101]
  for (k in pos1) {
    if (k %in% c(103,105,107))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[7], lwd = lwd.means)
    else if (k %in% c(109,111))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[8], lwd = lwd.means)
    else if (k %in% c(113,115,117))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[9], lwd = lwd.means)
  }
  
  pos = 1
  for (j in 1:6) {
    lines(x = c(line_X_axis[pos]-1,line_X_axis[pos+15]), y = rep(median(fff$bar[pos:(15+pos)]),2))
    pos = pos + 17
  }
  lines(x = c(line_X_axis[103]-2,line_X_axis[107]), y = rep(median(fff$bar[c(103,105,107)]),2))
  lines(x = c(line_X_axis[109]-2,line_X_axis[111]), y = rep(median(fff$bar[c(109,111)]),2))
  lines(x = c(line_X_axis[113]-2,line_X_axis[117]), y = rep(median(fff$bar[c(113,115,117)]),2))

  title(main=plot_main, cex.main = cex.main)
  if (is.na(labels)) {
    if (log) labels = 10**at else labels = at
    axis(side = 2, at = c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) c(2:9) * 10**i)))), 
         labels = c(rep('',length(c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) c(2:9) * 10**i))))))), 
         las = 2, pos = -0.5, tck = -0.02)
    axis(side = 2, at = c(at), labels = c(as.character(labels)), las = 2, pos = -0.5, tck = -0.03)
    #axis(side = 2, col = 'white', tcl = 0, labels = NA, lwd = 3, pos = -0.5)
    if (log) lines(x = c(-0.5,148), y = c(0,0),lty = 2, col = 'red')
    else  lines(x = c(-0.5,148), y = c(1,1),lty = 2, col = 'red')
  }  else {
    axis(side = 2, at = c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) c(2:9) * 10**i)))), 
         labels = c(rep('',length(c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) c(2:9) * 10**i))))))), 
         las = 2, pos = -0.5, tck = -0.02)
    axis(side = 2, at = c(at), labels = c(as.character(labels)), las = 2, pos = -0.5, tck = -0.03)
    #axis(side = 2, col = 'white', lwd = 3, tcl = 0, labels = NA, pos = -0.5)
    if (log) lines(x = c(-0.5,148), y = c(0,0),lty = 2, col = 'red')
    else  lines(x = c(-0.5,148), y = c(1,1),lty = 2, col = 'red')
  }
  
  
}

