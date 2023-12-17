library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(circlize)
library(ComplexHeatmap)
library(AnnotationDbi)
library(karyoploteR)
library(chicane)
library(GenomicInteractions)
library(Gviz)
library(plyr)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(ggpubr)
library(grid)
library(dplyr)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

gencode.v41.annotation.gene = data.frame(fread("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf"))
gencode.v41.annotation.gene$gene_id = sapply(gencode.v41.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][1], " ")[[1]][2])})
gencode.v41.annotation.gene$gene_type = sapply(gencode.v41.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][2], " ")[[1]][2])})
gencode.v41.annotation.gene$gene_name = sapply(gencode.v41.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][3], " ")[[1]][2])})
gencode.v41.annotation.gene = gencode.v41.annotation.gene[,-9]

Gr_gencode.v41.annotation.gene <- GRanges(
  seqnames = Rle(as.character(gencode.v41.annotation.gene[,1])),
  ranges = IRanges(gencode.v41.annotation.gene[,4], end = gencode.v41.annotation.gene[,5], names = c(1:nrow(gencode.v41.annotation.gene))),
  strand = Rle(strand('*')),
  gene_id = gencode.v41.annotation.gene[,9],
  gene_type = gencode.v41.annotation.gene[,10],
  gene_name = gencode.v41.annotation.gene[,11])


############### Figure 1a (compartment_analysis.r)
generate_karyoploteR_data <- function(tags, genome_gr, window_size, amplifier, threshold, names) {
  genome_window <- tileGenome(seqinfo(genome_gr), tilewidth = window_size, cut.last.tile.in.chrom = T)
  profile_list <- list()
  for(i in 1:length(tags)){
    profile_tmp  <- countOverlaps(genome_window, tags[[i]]) * amplifier[i]
    if (!is.na(threshold)){
      profile_tmp[profile_tmp > threshold[i]] <- threshold[i]
    }
    profile_list <- c(profile_list, list(profile_tmp))
  }
  genome_window@elementMetadata <- setNames(DataFrame(profile_list), nm = names)
  return(data.frame(genome_window))
}

plot_compartment_repetitive_RAL_tracks <- function(input_bedpe_file="",
                                                   eigen_file,
                                                   repetitive_RAL_file,
                                                   sample_name,
                                                   chromosomes,
                                                   resolution=500000,
                                                   tick_dist=50000000,
                                                   minor_tick_dist=10000000,
                                                   normalized = F,
                                                   normalized_by_RAL = F,
                                                   minoverlap = "30",
                                                   add_on_label = "", # additional custom label to be added to the figure filenames
                                                   only_interchromosomal_margi = F
)
{
  if (add_on_label == "_with_global_RAL"){
    
    ### Load iMARGI data
    input_bedpe = data.frame(fread(input_bedpe_file))
    input_data = input_bedpe[which(input_bedpe[,4] %in% chromosomes),1:6]
    rm(input_bedpe)
    
    input_data[,2] = input_data[,2] + 1 # 1-based 
    input_data[,5] = input_data[,5] + 1 # 1-based 
    input_data[,1] = as.character(input_data[,1])
    input_data[,4] = as.character(input_data[,4])
    colnames(input_data) = paste0("V",seq(1,6))
    
    if (only_interchromosomal_margi == T){
      input_data = input_data[which(input_data$V1 != input_data$V4),]
    }
    
    Gr_RNA_end = GRanges(
      seqnames = Rle(as.character(input_data[,1])),
      ranges = IRanges(input_data[,2], end = input_data[,3], names = c(1:nrow(input_data))),
      strand = Rle(strand('*')))
    
    Gr_DNA_end = GRanges(
      seqnames = Rle(as.character(input_data[,4])),
      ranges = IRanges(input_data[,5], end = input_data[,6], names = c(1:nrow(input_data))),
      strand = Rle(strand('*')))
    
    rm(input_data)
    
    cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
                                         genome_gr = Gr_hg38, 
                                         window_size = resolution,
                                         amplifier = c(1,1), 
                                         threshold = NA, 
                                         names = c("RNA_end","DNA_end"))
    
    cov_data = cov_data[which(cov_data$seqnames != "chrY"),]
  }
  
  ### Load compartment data
  options(scipen=999)
  compartment = read.table(eigen_file, header = T)
  compartment[which(compartment$eigen == 0), "eigen"] = NA
  
  compartment$eigen_pos = compartment$eigen
  compartment$eigen_neg = compartment$eigen
  
  compartment[,4][which(compartment[,4]<0)] = NA
  compartment[,5][which(compartment[,5]>0)] = NA
  
  colnames(compartment) = c("chr", "coord", "eigen", "eigen_pos", "eigen_neg")
  compartment = compartment[which(compartment$chr != "chrY"),]
  
  
  ##### Repetitive elements
  hg38_repeatMasker = data.frame(fread("/mnt/extraids/SDSC_NFS/wenxingzhao/database/Repeat_masker_human/hg38.fa.out.bed"))
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V11 == "SINE/Alu"),]
  repeat_data = repeat_data[which(repeat_data$V1 %in% hg38_chromosomes),]
  Alu_norm_factor = sum(repeat_data$V3-repeat_data$V2) / sum(hg38_lengths)
  
  Gr_Alu = GRanges(
    seqnames = Rle(repeat_data[,1]),
    ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
    strand = Rle(strand(repeat_data[,6])))
  
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V11 == "LINE/L1"),]
  repeat_data = repeat_data[which(repeat_data$V1 %in% hg38_chromosomes),]
  L1_norm_factor = sum(repeat_data$V3-repeat_data$V2) / sum(hg38_lengths)
  
  Gr_L1 = GRanges(
    seqnames = Rle(repeat_data[,1]),
    ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
    strand = Rle(strand(repeat_data[,6])))
  
  cov_data_repetitive = generate_karyoploteR_data(tags = list(Gr_Alu,Gr_L1),
                                                  genome_gr = Gr_hg38, 
                                                  window_size = resolution,
                                                  amplifier = c(1,1), 
                                                  threshold = NA, 
                                                  names = c("Alu","L1"))
  cov_data_repetitive = cov_data_repetitive[which(cov_data_repetitive$seqnames!="chrY"),]
  
  cov_data_repetitive$log2_Alu_L1 = log2(cov_data_repetitive$Alu / cov_data_repetitive$L1)
  
  cov_data_repetitive$log2_Alu_L1_pos = cov_data_repetitive$log2_Alu_L1
  cov_data_repetitive[which(cov_data_repetitive$log2_Alu_L1_pos < 0),"log2_Alu_L1_pos"] = NA
  
  cov_data_repetitive$log2_Alu_L1_neg = cov_data_repetitive$log2_Alu_L1
  cov_data_repetitive[which(cov_data_repetitive$log2_Alu_L1_neg > 0),"log2_Alu_L1_neg"] = NA
  
  cor_index = which(!is.na(compartment$eigen) & !is.na(cov_data_repetitive$log2_Alu_L1) & !is.infinite(cov_data_repetitive$log2_Alu_L1), arr.ind = T)
  out_cor = c()
  out_cor = c(out_cor, cor(compartment$eigen[cor_index], cov_data_repetitive$log2_Alu_L1[cor_index], method = "spearman"))
  
  
  ######## Repetitive elements RAL
  cov_data_RAL = read.table(paste0(repetitive_RAL_file), header = T)
  cov_data_RAL = cov_data_RAL[which(cov_data_RAL$seqnames!="chrY"),]
  
  ##### Add coverage data without RAL from Alu and L1 elements
  cov_data$RAL_without_Alu_L1 = cov_data$DNA_end - (cov_data_RAL$RAL_Alu + cov_data_RAL$RAL_L1)
  
  ####
  if (normalized_by_RAL == T){
    cov_data_RAL$RAL_Alu = cov_data_RAL$RAL_Alu / cov_data$DNA_end * 100
    cov_data_RAL$RAL_L1 = cov_data_RAL$RAL_L1 / cov_data$DNA_end * 100
  }
  
  if (normalized == T){
    cov_data_RAL$RAL_Alu = cov_data_RAL$RAL_Alu / (Alu_norm_factor * 100)
    cov_data_RAL$RAL_L1 = cov_data_RAL$RAL_L1 / (L1_norm_factor * 100)
    cov_data_RAL$log2_Alu_L1 = log2(cov_data_RAL$RAL_Alu / cov_data_RAL$RAL_L1)
  }
  
  cov_data_RAL$log2_Alu_L1_pos = cov_data_RAL$log2_Alu_L1
  cov_data_RAL[which(cov_data_RAL$log2_Alu_L1_pos < 0),"log2_Alu_L1_pos"] = NA
  
  cov_data_RAL$log2_Alu_L1_neg = cov_data_RAL$log2_Alu_L1
  cov_data_RAL[which(cov_data_RAL$log2_Alu_L1_neg > 0),"log2_Alu_L1_neg"] = NA
  
  cor_index = which(!is.na(compartment$eigen) & !is.na(cov_data_RAL$log2_Alu_L1) & !is.infinite(cov_data_RAL$log2_Alu_L1), arr.ind = T)
  out_cor = c(out_cor, cor(compartment$eigen[cor_index], cov_data_RAL$log2_Alu_L1[cor_index], method = "spearman"))
  
  names(out_cor) = c("cor_rep_density","cor_rep_RAL")
  
  
  ##### SPIN states
  ma_H1 <- fread("/mnt/extraids/SDSC_NFS/wenxingzhao/database/4DN/SPIN/H1_new.SPIN.JAWG.25kb.9_state.bed")
  STATES <- c("Speckle",  "Interior_Act1","Interior_Act2","Interior_Act3","Interior_Repr1", "Interior_Repr2","Near_Lm1","Near_Lm2","Lamina")
  color_ <-  c("#8b254a","#c14e4c","#ec7b57","#f2b579","#dbd291","#a8d29f", "#5fbba2","#7d9a98","#54508b")
  H1_spin <- ma_H1 %>% dplyr::filter(!grepl("NAN",V4)) %>% dplyr::mutate(V2=V2+1) %>% 
    mutate(spinidx = seq(1,nrow(.))) %>% 
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T,
                                            seqnames.field = "V1",
                                            start.field = "V2",
                                            end.field = "V3")
  H1_spin$color <- color_[match(H1_spin$V4,STATES)]
  
  for (my_chr in chromosomes){
    options(scipen=999)
    print(paste0(sample_name, " - ", my_chr))
    start_coord = 0
    end_coord = hg38_lengths[my_chr]
    detail.region <- toGRanges(data.frame(my_chr, start_coord, end_coord))
    
    temp_cov_data_RAL = cov_data_RAL[which(cov_data_RAL$seqnames==my_chr),]
    
    y_text = 0.05 # track label height
    
    line_width = 1 # width of the line track
    label_cex = 0.4
    
    if (add_on_label == "_with_global_RAL"){
      # RAL area
      r1_DNA_end = 0.99
      r0_DNA_end = r1_DNA_end - 0.08
      
      # RAL area without ALu and L1
      # r1_RAL_without_Alu_L1 = r0_DNA_end - 0.06
      # r0_RAL_without_Alu_L1 = r1_RAL_without_Alu_L1 - 0.08
      
      # Alu area
      r1_Alu = r0_DNA_end - 0.06
      r0_Alu = r1_Alu - 0.08
      
      # L1 area
      r1_L1 = r0_Alu - 0.06
      r0_L1 = r1_L1 - 0.08
      
      # log2 RAL repetitive element area
      r1_repetitive_RAL = r0_L1 - 0.06
      r0_repetitive_RAL = r1_repetitive_RAL - 0.08
      
      # Compartment area
      r1_comp = r0_repetitive_RAL - 0.06
      r0_comp = r1_comp - 0.08
      
      # SPIN states area
      # r1_spin = r0_comp - 0.06
      # r0_spin = r1_spin - 0.03
      
    } else if (add_on_label == ""){
      # Alu area
      r1_Alu = 0.99
      r0_Alu = r1_Alu - 0.12
      
      # L1 area
      r1_L1 = r0_Alu - 0.06
      r0_L1 = r1_L1 - 0.12
      
      # log2 RAL repetitive element area
      r1_repetitive_RAL = r0_L1 - 0.06
      r0_repetitive_RAL = r1_repetitive_RAL - 0.12
      
      # Compartment area
      r1_comp = r0_repetitive_RAL - 0.06
      r0_comp = r1_comp - 0.12
      
      # SPIN states area
      r1_spin = r0_comp - 0.06
      r0_spin = r1_spin - 0.04
    }
    
    
    # RAL area
    if (is.na(y_max_DNA)) {
      y_max_DNA_end = round(max(cov_data[which(cov_data$seqnames==my_chr),"DNA_end"]))
      y_max_RAL_without_Alu_L1 = round(max(cov_data[which(cov_data$seqnames==my_chr),"RAL_without_Alu_L1"]))
    } else {
      y_max_DNA_end = y_max_DNA
    }
    
    # log2 repetitive element area
    # r1_repetitive = 0.99
    # r0_repetitive = r1_repetitive - 0.15
    # temp = cov_data_repetitive[which(cov_data_repetitive$seqnames==my_chr),]
    # ymin_repetitive = round(min(temp[which(!is.na(temp$log2_Alu_L1) & !is.infinite(temp$log2_Alu_L1)),"log2_Alu_L1"]),1)
    # ymax_repetitive = round(max(temp[which(!is.na(temp$log2_Alu_L1) & !is.infinite(temp$log2_Alu_L1)),"log2_Alu_L1"]),1)
    # y_repetitive = max(ymax_repetitive,-ymin_repetitive)
    
    # log2 RAL repetitive element area
    ymin_repetitive_RAL = round(min(temp_cov_data_RAL[which(!is.na(temp_cov_data_RAL$log2_Alu_L1) & !is.infinite(temp_cov_data_RAL$log2_Alu_L1)),"log2_Alu_L1"]),1)
    ymax_repetitive_RAL = round(max(temp_cov_data_RAL[which(!is.na(temp_cov_data_RAL$log2_Alu_L1) & !is.infinite(temp_cov_data_RAL$log2_Alu_L1)),"log2_Alu_L1"]),1)
    y_repetitive_RAL = max(ymax_repetitive_RAL,-ymin_repetitive_RAL)
    
    # Compartment area
    ymin_comp = round(min(compartment[which(compartment$chr==my_chr),"eigen_neg"][which(!is.na(compartment[which(compartment$chr==my_chr),"eigen_neg"]))]),1)
    ymax_comp = round(max(compartment[which(compartment$chr==my_chr),"eigen_pos"][which(!is.na(compartment[which(compartment$chr==my_chr),"eigen_pos"]))]),1)
    y_comp = max(ymax_comp,-ymin_comp)
    
    options(scipen=999)
    if (normalized == T){
      pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks_norm/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_norm",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    } else if (normalized_by_RAL == T){
      pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_norm_by_RAL",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    } else if (only_interchromosomal_margi == T){
      pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_inter_",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    } else {
      pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution),"/",sample_name,"_",my_chr,"_tracks_RAL",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    }
    
    plot_params <- getDefaultPlotParams(plot.type=1)
    plot_params$data1inmargin <- -40
    plot_params$ideogramheight <- 10
    plot_params$bottommargin <- 30
    plot_params$topmargin <- -5
    plot_params$leftmargin <- 0.12
    plot_params$rightmargin <- 0.02
    kp <- plotKaryotype(genome=Gr_hg38, plot.type=1, plot.params = plot_params, zoom = detail.region, cex=0.5)
    kpAddBaseNumbers(kp, tick.dist = tick_dist, tick.len = 6, tick.col="#4d4d4d", cex=0.5,
                     minor.tick.dist = minor_tick_dist, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
    #kpAddMainTitle(kp, main="test", col="red")
    
    if (add_on_label == "_with_global_RAL"){
      kpDataBackground(kp, data.panel = 1, r0=r0_DNA_end, r1=r1_DNA_end, color = "#ffffff")
      kpAxis(kp, ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, col="gray50", cex=label_cex, numticks = 3)
      # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
      #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
      kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$DNA_end,
              col="purple", ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, lwd=line_width)
      
      # kpDataBackground(kp, data.panel = 1, r0=r0_RAL_without_Alu_L1, r1=r1_RAL_without_Alu_L1, color = "#ffffff")
      # kpAxis(kp, ymin=0, ymax=y_max_RAL_without_Alu_L1, r0=r0_RAL_without_Alu_L1, r1=r1_RAL_without_Alu_L1, col="gray50", cex=label_cex, numticks = 3)
      # # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
      # #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
      # kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$RAL_without_Alu_L1,
      #         col="purple", ymin=0, ymax=y_max_RAL_without_Alu_L1, r0=r0_RAL_without_Alu_L1, r1=r1_RAL_without_Alu_L1, lwd=line_width)
    }
    
    # Alu barplot
    kpDataBackground(kp, data.panel = 1, r0=r0_Alu, r1=r1_Alu, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=round(max(temp_cov_data_RAL[which(!is.na(temp_cov_data_RAL$RAL_Alu)),"RAL_Alu"]),1), r0=r0_Alu, r1=r1_Alu, col="gray50", cex=label_cex, numticks = 3)
    kpBars(kp, chr=temp_cov_data_RAL$seqnames, x0=temp_cov_data_RAL$start, x1=temp_cov_data_RAL$end, y1=temp_cov_data_RAL$RAL_Alu,
           col="orange", ymin=0, ymax=max(temp_cov_data_RAL[which(!is.na(temp_cov_data_RAL$RAL_Alu)),"RAL_Alu"]), r0=r0_Alu, r1=r1_Alu, border = "orange")
    # kpText(kp, chr=cov_data_RAL$seqnames, x=mean(c(start(detail.region),end(detail.region))),
    #        y=y_text, col="#000000", r0=r1_Alu + 0.01, r1=r1_Alu + 0.02, labels="Alu_RAL", cex=label_cex)
    
    # L1 barplot
    kpDataBackground(kp, data.panel = 1, r0=r0_L1, r1=r1_L1, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=round(max(temp_cov_data_RAL[which(!is.na(temp_cov_data_RAL$RAL_L1)),"RAL_L1"]),1), r0=r0_L1, r1=r1_L1, col="gray50", cex=label_cex, numticks = 3)
    kpBars(kp, chr=temp_cov_data_RAL$seqnames, x0=temp_cov_data_RAL$start, x1=temp_cov_data_RAL$end, y1=temp_cov_data_RAL$RAL_L1,
           col="blue", ymin=0, ymax=max(temp_cov_data_RAL[which(!is.na(temp_cov_data_RAL$RAL_L1)),"RAL_L1"]), r0=r0_L1, r1=r1_L1, border = "blue")
    # kpText(kp, chr=cov_data_RAL$seqnames, x=mean(c(start(detail.region),end(detail.region))),
    #        y=y_text, col="#000000", r0=r1_L1 + 0.01, r1=r1_L1 + 0.02, labels="L1_RAL", cex=label_cex)
    
    # Repetitive log2ratio RAL barplot
    kpDataBackground(kp, data.panel = 1, r0=r0_repetitive_RAL, r1=r1_repetitive_RAL, color = "#ffffff")
    kpAxis(kp, ymin=-y_repetitive_RAL, ymax=y_repetitive_RAL, r0=r0_repetitive_RAL, r1=r1_repetitive_RAL, col="gray50", cex=label_cex, numticks = 3)
    kpBars(kp, chr=temp_cov_data_RAL$seqnames, x0=temp_cov_data_RAL$start, x1=temp_cov_data_RAL$end, y1=temp_cov_data_RAL$log2_Alu_L1_pos,
           col="orange", ymin=-y_repetitive_RAL, ymax=y_repetitive_RAL, r0=r0_repetitive_RAL, r1=r1_repetitive_RAL, border = "orange")
    kpBars(kp, chr=temp_cov_data_RAL$seqnames, x0=temp_cov_data_RAL$start, x1=temp_cov_data_RAL$end, y1=temp_cov_data_RAL$log2_Alu_L1_neg,
           col="blue", ymin=-y_repetitive_RAL, ymax=y_repetitive_RAL, r0=r0_repetitive_RAL, r1=r1_repetitive_RAL, border = "blue")
    # kpText(kp, chr=cov_data_RAL$seqnames, x=mean(c(start(detail.region),end(detail.region))),
    #        y=y_text, col="#000000", r0=r1_repetitive_RAL + 0.01, r1=r1_repetitive_RAL + 0.02, labels="log2(Alu_RAL/L1_RAL)", cex=label_cex)
    
    # Repetitive log2ratio barplot
    # kpDataBackground(kp, data.panel = 1, r0=r0_repetitive, r1=r1_repetitive, color = "#ffffff")
    # kpAxis(kp, ymin=-y_repetitive, ymax=y_repetitive, r0=r0_repetitive, r1=r1_repetitive, col="gray50", cex=label_cex, numticks = 3)
    # kpBars(kp, chr=cov_data_repetitive$seqnames, x0=cov_data_repetitive$start, x1=cov_data_repetitive$end, y1=cov_data_repetitive$log2_Alu_L1_pos,
    #        col="black", ymin=-y_repetitive, ymax=y_repetitive, r0=r0_repetitive, r1=r1_repetitive, border = "black")
    # kpBars(kp, chr=cov_data_repetitive$seqnames, x0=cov_data_repetitive$start, x1=cov_data_repetitive$end, y1=cov_data_repetitive$log2_Alu_L1_neg,
    #        col="gray40", ymin=-y_repetitive, ymax=y_repetitive, r0=r0_repetitive, r1=r1_repetitive, border = "gray40")
    # # kpText(kp, chr=cov_data_RAL$seqnames, x=mean(c(start(detail.region),end(detail.region))),
    # #        y=y_text, col="#000000", r0=r1_repetitive + 0.01, r1=r1_repetitive + 0.02, labels="log2(Alu/L1)", cex=label_cex)
    
    # Compartment barplot
    kpDataBackground(kp, data.panel = 1, r0=r0_comp, r1=r1_comp, color = "#ffffff")
    kpAxis(kp, ymin=-y_comp, ymax=y_comp, r0=r0_comp, r1=r1_comp, col="gray50", cex=label_cex, numticks = 3)
    kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_pos,
           col="#E41A1C", ymin=-y_comp, ymax=y_comp, r0=r0_comp, r1=r1_comp, border = "#E41A1C")
    kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_neg,
           col="#0cad01", ymin=-y_comp, ymax=y_comp, r0=r0_comp, r1=r1_comp, border = "#0cad01")
    # kpText(kp, chr=cov_data_RAL$seqnames, x=mean(c(start(detail.region),end(detail.region))),
    #        y=y_text, col="#000000", r0=r1_comp + 0.01, r1=r1_comp + 0.02, labels="A/B compartments", cex=label_cex)
    
    # SPIN states
    # kpPlotRegions(kp, data=H1_spin, col=H1_spin$color, r0=r0_spin , r1=r1_spin, data.panel=1)
    
    
    dev.off()
    
  }
  
  return(out_cor)
}

plot_compartment_repetitive_RAL_tracks(input_bedpe_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                       eigen_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_control/H1_control_500000.txt",
                                       repetitive_RAL_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/500000/H1_control_200k_Alu_L1_RAL_minoverlap_30.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr11"),
                                       resolution=500000,
                                       tick_dist=20000000,
                                       minor_tick_dist=10000000,
                                       add_on_label = "_with_global_RAL")

plot_compartment_repetitive_RAL_tracks(input_bedpe_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                       eigen_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_control/H1_control_500000.txt",
                                       repetitive_RAL_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/500000/H1_control_200k_Alu_L1_RAL_minoverlap_30_inter.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr11"),
                                       resolution=500000,
                                       tick_dist=20000000,
                                       minor_tick_dist=10000000,
                                       add_on_label = "_with_global_RAL",
                                       only_interchromosomal_margi = T)



########### RAL of specific genes
plot_compartment_gene_RAL_tracks <- function(input_bedpe_file="",
                                             eigen_file,
                                             sample_name,
                                             genes,
                                             chromosomes,
                                             resolution=500000,
                                             tick_dist=50000000,
                                             minor_tick_dist=10000000,
                                             normalized = F,
                                             normalized_by_RAL = F,
                                             minoverlap = "30",
                                             only_interchromosomal_margi = F
)
{
  
  
  ### Load compartment data
  options(scipen=999)
  compartment = read.table(eigen_file, header = T)
  compartment[which(compartment$eigen == 0), "eigen"] = NA
  
  compartment$eigen_pos = compartment$eigen
  compartment$eigen_neg = compartment$eigen
  
  compartment[,4][which(compartment[,4]<0)] = NA
  compartment[,5][which(compartment[,5]>0)] = NA
  
  colnames(compartment) = c("chr", "coord", "eigen", "eigen_pos", "eigen_neg")
  compartment = compartment[which(compartment$chr != "chrY"),]
  
  
  ### Load iMARGI data
  input_bedpe = data.frame(fread(input_bedpe_file))
  
  # temp_chr = as.character(seqnames(Gr_annotations_edb_hg38[mcols(Gr_annotations_edb_hg38)["SYMBOL"][,1] %in% genes]))
  # input_data = input_bedpe[which(input_bedpe[,1] %in% temp_chr &
  #                                  input_bedpe[,4] %in% temp_chr),1:6] # we are interested only in intrachromosomal RAL of the selected genes
  
  input_data = input_bedpe[which(input_bedpe[,1] == input_bedpe[,4]),1:6] # only intrachromosomal
  rm(input_bedpe)
  
  input_data[,2] = input_data[,2] + 1 # 1-based 
  input_data[,5] = input_data[,5] + 1 # 1-based 
  input_data[,1] = as.character(input_data[,1])
  input_data[,4] = as.character(input_data[,4])
  colnames(input_data) = paste0("V",seq(1,6))
  
  if (only_interchromosomal_margi == T){
    input_data = input_data[which(input_data$V1 != input_data$V4),]
  }
  
  Gr_RNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,1])),
    ranges = IRanges(input_data[,2], end = input_data[,3], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  Gr_DNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,4])),
    ranges = IRanges(input_data[,5], end = input_data[,6], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  rm(input_data)
  
  cov_data_full = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
                                           genome_gr = Gr_hg38,
                                           window_size = resolution,
                                           amplifier = c(1,1),
                                           threshold = NA,
                                           names = c("RNA_end","DNA_end"))
  
  # overlaps = countOverlaps(Gr_annotations_edb_hg38, Gr_RNA_end, ignore.strand = T)
  # overlaps = sort(overlaps, decreasing = T)
  # genes = data.frame(Gr_annotations_edb_hg38[names(overlaps[1:10])])[,"SYMBOL"]
  
  for (my_gene in genes){
    
    Gr_gene = Gr_annotations_edb_hg38[mcols(Gr_annotations_edb_hg38)["SYMBOL"][,1] == my_gene]
    overlaps = countOverlaps(Gr_RNA_end, Gr_gene, ignore.strand = T)
    
    Gr_RNA_end_temp = Gr_RNA_end[overlaps>0]
    Gr_DNA_end_temp = Gr_DNA_end[overlaps>0]
    
    cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end_temp,Gr_DNA_end_temp),
                                         genome_gr = Gr_hg38,
                                         window_size = resolution,
                                         amplifier = c(1,1),
                                         threshold = NA,
                                         names = c("RNA_end","DNA_end"))
    
    cov_data = cov_data[which(cov_data$seqnames != "chrY"),]
    
    
    options(scipen=999)
    my_chr = as.character(seqnames(Gr_gene))
    print(paste0(sample_name, " - ", my_gene))
    start_coord = 0
    end_coord = hg38_lengths[my_chr]
    detail.region <- toGRanges(data.frame(my_chr, start_coord, end_coord))
    
    ### Save to file source data
    options(scipen=999)
    resolution_char = as.character(resolution)
    
    write.table(cov_data[which(cov_data$seqnames == my_chr),c(1,2,3,7)] %>% `colnames<-`(c("chr", "start", "end","RAL")),
                paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution_char,"/",sample_name,"_",my_gene,"_tracks_RAL_minoverlap_", minoverlap,".txt"),
                row.names = F, col.names = T, quote = F, sep = "\t")
    write.table(cov_data[which(cov_data_full$seqnames == my_chr),c(1,2,3,7)] %>% `colnames<-`(c("chr", "start", "end","cRAL")),
                paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution_char,"/",sample_name,"_",my_gene,"_tracks_cRAL_minoverlap_", minoverlap,".txt"),
                row.names = F, col.names = T, quote = F, sep = "\t")
    write.table(compartment[which(compartment$chr == my_chr),1:3],
                paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution_char,"/",sample_name,"_",my_gene,"_tracks_PC1_minoverlap_", minoverlap,".txt"),
                row.names = F, col.names = T, quote = F, sep = "\t")
    
    # temp_cov_data_RAL = cov_data_RAL[which(cov_data_RAL$seqnames==my_chr),]
    
    y_text = 0.05 # track label height
    line_width = 1 # width of the line track
    label_cex = 0.4
    
    # RAL area
    r1_DNA_end = 0.8
    r0_DNA_end = r1_DNA_end - 0.08
    
    r1_DNA_end_zoom = r0_DNA_end - 0.08
    r0_DNA_end_zoom = r1_DNA_end_zoom - 0.08
    
    r1_DNA_end_full = r0_DNA_end_zoom - 0.08
    r0_DNA_end_full = r1_DNA_end_full - 0.08
    
    # Compartment area
    r1_comp = r0_DNA_end_full - 0.09
    r0_comp = r1_comp - 0.08
    
    
    # RAL area
    if (is.na(y_max_DNA)) {
      y_max_DNA_end = round(max(cov_data[which(cov_data$seqnames==my_chr),"DNA_end"]))
    } else {
      y_max_DNA_end = y_max_DNA
    }
    
    if (is.na(y_max_DNA)) {
      y_max_DNA_end_full = round(max(cov_data_full[which(cov_data_full$seqnames==my_chr),"DNA_end"]))
    } else {
      y_max_DNA_end_full = y_max_DNA
    }
    
    # Compartment area
    ymin_comp = round(min(compartment[which(compartment$chr==my_chr),"eigen_neg"][which(!is.na(compartment[which(compartment$chr==my_chr),"eigen_neg"]))]),1)
    ymax_comp = round(max(compartment[which(compartment$chr==my_chr),"eigen_pos"][which(!is.na(compartment[which(compartment$chr==my_chr),"eigen_pos"]))]),1)
    y_comp = max(ymax_comp,-ymin_comp)
    
    # options(scipen=999)
    # if (normalized == T){
    #   pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks_norm/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_norm",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    # } else if (normalized_by_RAL == T){
    #   pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_norm_by_RAL",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    # } else if (only_interchromosomal_margi == T){
    #   pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_inter_",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    # } else {
    #   pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution),"/",sample_name,"_",my_chr,"_tracks_RAL",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    # }
    
    dir.create(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution_char), recursive = T, showWarnings = F)
    
    pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution_char,"/",sample_name,"_",my_gene,"_tracks_RAL_minoverlap_", minoverlap,".pdf"), width = 3.3, height = 2.5)
    
    plot_params <- getDefaultPlotParams(plot.type=1)
    plot_params$data1inmargin <- -40
    plot_params$ideogramheight <- 7
    plot_params$bottommargin <- 30
    plot_params$topmargin <- -5
    plot_params$leftmargin <- 0.12
    plot_params$rightmargin <- 0.02
    kp <- plotKaryotype(genome=Gr_hg38, plot.type=1, plot.params = plot_params, zoom = detail.region, cex=label_cex)
    kpAddBaseNumbers(kp, tick.dist = tick_dist, tick.len = 6, tick.col="#4d4d4d", cex=label_cex,
                     minor.tick.dist = minor_tick_dist, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
    #kpAddMainTitle(kp, main="test", col="red")
    
    # RAL plot
    temp_y_max = y_max_DNA_end
    kpDataBackground(kp, data.panel = 1, r0=r0_DNA_end, r1=r1_DNA_end, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=temp_y_max, r0=r0_DNA_end, r1=r1_DNA_end, col="gray50", cex=label_cex, numticks = 3)
    # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
    #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
    kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$DNA_end,
            col="purple", ymin=0, ymax=temp_y_max, r0=r0_DNA_end, r1=r1_DNA_end, lwd=line_width)
    
    
    temp_y_max = round(y_max_DNA_end / 5)
    cov_data_temp = cov_data
    cov_data_temp[which(cov_data_temp$DNA_end >= temp_y_max), "DNA_end"] = temp_y_max
    
    kpDataBackground(kp, data.panel = 1, r0=r0_DNA_end_zoom, r1=r1_DNA_end_zoom, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=temp_y_max, r0=r0_DNA_end_zoom, r1=r1_DNA_end_zoom, col="gray50", cex=label_cex, numticks = 3)
    # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
    #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
    kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data_temp$DNA_end,
            col="purple", ymin=0, ymax=temp_y_max, r0=r0_DNA_end_zoom, r1=r1_DNA_end_zoom, lwd=line_width)
    
    
    temp_y_max = y_max_DNA_end_full
    kpDataBackground(kp, data.panel = 1, r0=r0_DNA_end_full, r1=r1_DNA_end_full, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=temp_y_max, r0=r0_DNA_end_full, r1=r1_DNA_end_full, col="gray50", cex=label_cex, numticks = 3)
    # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
    #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
    kpLines(kp, chr=cov_data_full$seqnames, x=rowMeans(cov_data_full[,2:3]), y=cov_data_full$DNA_end,
            col="purple", ymin=0, ymax=temp_y_max, r0=r0_DNA_end_full, r1=r1_DNA_end_full, lwd=line_width)
    
  
    # Compartment barplot
    kpDataBackground(kp, data.panel = 1, r0=r0_comp, r1=r1_comp, color = "#ffffff")
    kpAxis(kp, ymin=-y_comp, ymax=y_comp, r0=r0_comp, r1=r1_comp, col="gray50", cex=label_cex, numticks = 3)
    kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_pos,
           col="#E41A1C", ymin=-y_comp, ymax=y_comp, r0=r0_comp, r1=r1_comp, border = "#E41A1C")
    kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_neg,
           col="#0cad01", ymin=-y_comp, ymax=y_comp, r0=r0_comp, r1=r1_comp, border = "#0cad01")
    # kpText(kp, chr=cov_data_RAL$seqnames, x=mean(c(start(detail.region),end(detail.region))),
    #        y=y_text, col="#000000", r0=r1_comp + 0.01, r1=r1_comp + 0.02, labels="A/B compartments", cex=label_cex)
    
    dev.off()
    
  }
  
}

plot_compartment_gene_RAL_tracks(input_bedpe_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                       eigen_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_control/H1_control_500000.txt",
                                       sample_name="H1_control_200k",
                                       genes=c("PVT1","JARID2"),
                                       resolution=500000,
                                       tick_dist=20000000,
                                       minor_tick_dist=10000000)

plot_compartment_gene_RAL_tracks(input_bedpe_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                 eigen_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_control/H1_control_500000.txt",
                                 sample_name="H1_control_200k",
                                 resolution=500000,
                                 tick_dist=20000000,
                                 minor_tick_dist=10000000)

# plot_compartment_gene_RAL_tracks(input_bedpe_file = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe.gz",
#                                  eigen_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_control/H1_control_500000.txt",
#                                  sample_name="H1_control_1k",
#                                  genes=c("NCAM2", "NRIP1", "ATXN1", "ZMYND8"),
#                                  resolution=500000,
#                                  tick_dist=20000000,
#                                  minor_tick_dist=10000000)


make_gene_RAL_heatmap <- function(input_bedpe_file,
                                             sample_name,
                                             genes,
                                             resolution=500000,
                                             minoverlap = "30",
                                             only_interchromosomal_margi = F
)
{
  
  
  ### Load iMARGI data
  input_bedpe = data.frame(fread(input_bedpe_file))
  
  # temp_chr = as.character(seqnames(Gr_annotations_edb_hg38[mcols(Gr_annotations_edb_hg38)["SYMBOL"][,1] %in% genes]))
  # input_data = input_bedpe[which(input_bedpe[,1] %in% temp_chr &
  #                                  input_bedpe[,4] %in% temp_chr),1:6] # we are interested only in intrachromosomal RAL of the selected genes
  
  input_data = input_bedpe[which(input_bedpe[,1] == input_bedpe[,4]),1:6] # only intrachromosomal
  rm(input_bedpe)
  
  input_data[,2] = input_data[,2] + 1 # 1-based 
  input_data[,5] = input_data[,5] + 1 # 1-based 
  input_data[,1] = as.character(input_data[,1])
  input_data[,4] = as.character(input_data[,4])
  colnames(input_data) = paste0("V",seq(1,6))
  
  if (only_interchromosomal_margi == T){
    input_data = input_data[which(input_data$V1 != input_data$V4),]
  }
  
  Gr_RNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,1])),
    ranges = IRanges(input_data[,2], end = input_data[,3], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  Gr_DNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,4])),
    ranges = IRanges(input_data[,5], end = input_data[,6], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  rm(input_data)
  
  # Top10 expressed genes
  # overlaps = countOverlaps(Gr_annotations_edb_hg38, Gr_RNA_end, ignore.strand = T)
  # overlaps = sort(overlaps, decreasing = T)
  # genes = data.frame(Gr_annotations_edb_hg38[names(overlaps[1:10])])[,"SYMBOL"]
  
  # Bottom10 expressed genes
  # overlaps = countOverlaps(Gr_gencode.v41.annotation.gene, Gr_RNA_end, ignore.strand = T)
  # overlaps = sort(overlaps)
  # overlaps = overlaps[which(overlaps > 1000)]
  # 
  # temp = data.frame(Gr_gencode.v41.annotation.gene)
  # temp = temp[which(temp$width > 200000),]
  #   
  # genes = temp[names(overlaps) %in% rownames(temp),][1:3,"gene_name"]
  
  
  for (my_gene in genes){
    
    Gr_gene = Gr_gencode.v41.annotation.gene[mcols(Gr_gencode.v41.annotation.gene)["gene_name"][,1] == my_gene]
    
    for (i in c(10000)){
      Gr_gene_tile = tile(Gr_gene, width = i)[[1]]
      names(Gr_gene_tile) = seq(1:length(Gr_gene_tile))
      
      temp_fun <- function(x){
        overlaps = countOverlaps(Gr_RNA_end, Gr_gene_tile[x], ignore.strand = T)
        
        Gr_RNA_end_temp = Gr_RNA_end[overlaps>0]
        Gr_DNA_end_temp = Gr_DNA_end[overlaps>0]
        
        cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end_temp,Gr_DNA_end_temp),
                                             genome_gr = Gr_hg38,
                                             window_size = resolution,
                                             amplifier = c(1,1),
                                             threshold = NA,
                                             names = c("RNA_end","DNA_end"))
        cov_data = cov_data[which(cov_data$seqnames == as.character(seqnames(Gr_gene))),]
        return(cov_data$DNA_end)
      }
      
      out_list = pbmcapply::pbmclapply(1:length(Gr_gene_tile), temp_fun, mc.cores = 32)
      df_heatmap = do.call("rbind", out_list)
      
      options(scipen=999)
      dir.create(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",as.character(resolution),"/heatmap/data/"), recursive = T, showWarnings = F)
      filename = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",as.character(resolution),"/heatmap/data/",sample_name,"_",my_gene,"_", i, "_heatmap_RAL_minoverlap_", minoverlap,".txt")
      write.table(df_heatmap, filename, row.names = F, col.names = F, sep = "\t", quote = F)
    
    }

  }
  
}


make_gene_RAL_heatmap(input_bedpe_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                  sample_name="H1_control_200k",
                                  genes=c("GRID2",  "PTPRG",  "ROR1",   "PVT1",   "UNC5D",  "JARID2", "FOXN3",  "ADCY2",  "PBX1",   "SHANK2", "NCAM2", "NRIP1", "ATXN1", "ZMYND8"))



plot_gene_RAL_heatmap <- function(sample_name,
                                  resolution=500000,
                                  minoverlap=30,
                                  genes,
                                  gene_window_sizes){
  for (my_gene in genes){
    for (i in gene_window_sizes){
      
      options(scipen=999)
      df_heatmap = as.matrix(read.table(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution,"/heatmap/data/",sample_name,"_",my_gene,"_", i, "_heatmap_RAL_minoverlap_", minoverlap,".txt")))
      
      col_fun = colorRamp2(c(0, 100, 200), c("blue","white","red"))
      # col_fun = colorRamp2(c(0, 20, 40), c("blue","white","red")) # to highlight the stripes farther from the genes
      p <- Heatmap(df_heatmap, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, border = T, show_column_names = F, show_row_names = F, heatmap_legend_param = list(title = ""))

      # p <- Heatmap(df_heatmap, cluster_rows = FALSE, cluster_columns = FALSE, border = T, show_column_names = F, show_row_names = F, heatmap_legend_param = list(title = ""))
      dir.create(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution,"/heatmap/figures/"), showWarnings = F, recursive = T)
      png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/genes/",resolution,"/heatmap/figures/",sample_name,"_",my_gene,"_", i, "_heatmap_RAL_minoverlap_", minoverlap,".png"), width = 5, height = 1.4, res = 300, units = "in")
      print(p)
      dev.off()
    }
  }
}


plot_gene_RAL_heatmap(sample_name="H1_control_200k",
                                  genes=c("ZMYND8","PVT1","JARID2","SHANK2"),
                                  gene_window_sizes=c(10000))


c("GRID2",  "PTPRG",  "ROR1",   "PVT1",   "UNC5D",  "JARID2", "FOXN3",  "ADCY2",  "PBX1",   "SHANK2", "NCAM2", "NRIP1", "ATXN1", "ZMYND8")


plot_gene_RAL_heatmap(sample_name="control",
                      genes=genes,
                      gene_window_sizes=c(10000))


plot_gene_RAL_heatmap(sample_name="control",
                      genes=c("UPP2", "TSPAN9", "TFPI", "SLC25A13", "KLHL13", "HS3ST3B1"),
                      gene_window_sizes=c(10000))


############### Figure 1c (MARGI_heatmap.r)
my_sample = "iMARGI_HFF_control"
loop_file_path_temp = "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_HFF/"
HFF_stripes = read.table("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/stripe_cellLine_narrow/HFF.raw.2dBed")

p <- make_complex_heatmap_margi(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr6_chr6.txt"),
  sparse=T,
  chr_row="chr6",
  chr_col="chr6",
  bin_size=10000,
  select_coord=rep(c(15500000,17500000),2),
  contacts_cutoffs=c(0,0.99),
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = NA,
  loops_file_path = NA,
  EP_loops_file = NA,
  feature_color = "blue",
  RNA_domain_coord = c(16285686,16768731,16246982,17097795))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_RNAdomain_A2.pdf", width = 8, height = 8)
print(p)
dev.off()

# Full chromosome for presentation (no paper plot)
p <- make_complex_heatmap_margi(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/500000_filter1k/chr6_chr6.txt"),
  sparse=T,
  chr_row="chr6",
  chr_col="chr6",
  bin_size=500000,
  select_coord=NA,
  contacts_cutoffs=c(0,0.99),
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = NA,
  loops_file_path = NA,
  EP_loops_file = NA,
  feature_color = "blue",
  RNA_domain_coord = rep(c(15500000,17500000),2))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/chr6_presentation.pdf", width = 8, height = 8)
print(p)
dev.off()



TxDb.Hsapiens.gencode.24 = makeTxDbFromGFF('/home/frankyan/research/refGenome/standard4DN/GRCh38/gencode.v24.primary_assembly.annotation.gtf')

make_annotation_tracks_2 <- function(chr,
                                   start,
                                   end,
                                   sample_name,
                                   panel="A"){
  
  chr.region <- toGRanges(paste0(chr,":",start,"-",end))
  
  kp <- plotKaryotype(zoom = chr.region)
  genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.gencode.24,
                                      karyoplot=kp,
                                      plot.transcripts = T, 
                                      plot.transcripts.structure = T)
  
  mcols(genes.data$genes)[,"gene_id"] = gsub("\\..*", "", mcols(genes.data$genes)[,"gene_id"])
  protein_coding_genes = annotations_edb_hg38[which(annotations_edb_hg38$GENEID %in% mcols(genes.data$genes)[,"gene_id"] &
                                                      annotations_edb_hg38$GENEBIOTYPE == "protein_coding"), "SYMBOL"]
  
  mcols(genes.data$genes)[,"gene_id"] = annotations_edb_hg38[which(annotations_edb_hg38$GENEID %in% mcols(genes.data$genes)[,"gene_id"]), "SYMBOL"]
  
  
  a = which(mcols(genes.data[["genes"]])["gene_id"][,1] == "ATXN1")
  genes.data[["genes"]] = genes.data[["genes"]][a]
  
  genes.data <- mergeTranscripts(genes.data)
  
  gene_track_border = 0.75
  step = 0.8 / length(chipseq_tracks_list) # height of each track
  track_borders = seq(0.8,0, by=-step)
  
  #png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_name, "_", chr, "_", start, "_", end, "_tracks.png"), width = 8.5, height = 6, units = "in", res = 500)
  

  pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/",sample_name, "_", chr,"_", start, "_", end, "_tracks.pdf"), width = 3, height = 1)
  
  plot_params <- getDefaultPlotParams(plot.type=1)
  plot_params$data1inmargin <- 5
  plot_params$ideogramheight <- 10
  plot_params$bottommargin <- 30
  plot_params$topmargin <- 0
  plot_params$leftmargin <- 0.15
  plot_params$rightmargin <- 0.02
  
  kp <- plotKaryotype(zoom = chr.region, plot.params = plot_params)
  kpPlotGenes(kp, data=genes.data, r1=1, r0=0.6, gene.name.position = "left", gene.name.cex = 1, mark.height = 0.4)
  dev.off()
  
}

make_annotation_tracks_2(chr="chr6", start=15500000, end=17500000, sample_name="HFF", panel="A")


############### Figure 1d
temp = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/B_stripe_cellLine_upset.csv")

library(UpSetR)
png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_RNAdomain_B.png", width = 5, height = 3.5, units = "in", res = 400)
upset(temp, 
      sets = c("K562","HFF","H1"),
      sets.x.label = "Number of caRNA Domains",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      #text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
      text.scale = c(1.8, 1.8, 0, 0, 1.6, 0),
      mainbar.y.label = "Shared caRNA domains",
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 3, line.size = 0.5)
dev.off()


############### Figure 1e
temp = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/D_stripe_cellLine_heightWidth_boxPlot_wideVersion.csv", header = F)

temp = t(temp)[-1,]
temp = data.frame(temp)

temp_H1 = temp[,1:2]
colnames(temp_H1) = c("Height", "Width")
temp_H1$sample = "H1"

temp_HFF = temp[,3:4]
colnames(temp_HFF) = c("Height", "Width")
temp_HFF$sample = "HFF"

temp_K562 = temp[,5:6]
colnames(temp_K562) = c("Height", "Width")
temp_K562$sample = "K562"


df = rbind(temp_H1, temp_HFF, temp_K562)
df$Height = as.numeric(as.character(df$Height))
df$Width = as.numeric(as.character(df$Width))

write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_1e.txt", row.names = F, col.names = T, sep = "\t", quote = F)

library(scales)

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_RNAdomain_D.pdf", height = 1.6, width = 2)
ggplot(melt(df), aes(x=sample, y=value, fill=variable)) + 
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.2, position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  labs(x="", y="Rectangular block (kb)") +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none")
dev.off()

temp = df[which(df$sample == "H1"),"Height"]
IQR = as.numeric(quantile(temp, 0.75) - quantile(temp, 0.25))


####### Figure 2g
OR = c(1.268448943, 0.748604623)
err = c(0.094166895, 0.069751158)

df = data.frame(x = OR,
                y = c(2,1),
                err = err)

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HERV_H_paper_data/margi/result/Figure_2g.pdf", width = 1.5, height = 0.8)
ggplot(df, aes(x=x, y=y)) +
  geom_errorbar(aes(xmin=exp(log(x)-err), xmax=exp(log(x)+err)), width=.3, position=position_dodge(.9)) +
  geom_point(size=1.5, color="blue") +
  scale_x_continuous(trans='log', breaks = c(0.5,1,1:5), labels = as.character(c(0.5,1,1:5)), limits = c(0.5,1.5)) +
  #ylim(c(0.5,4.5))+
  # xlab("Odds ratio (log scale)") +
  # ylab("") +
  guides(y = "none") +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_bw() +
  theme(text = element_text(size=6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed")
dev.off()

temp_mat = matrix(data=c(207,401,912,2241), nrow=2, ncol=2, byrow = T)
chisq.test(as.table(temp_mat))


####### Figure 2e
# df = readRDS("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/RNase_Ctrl_boxplot.rds")
# colnames(df) = c("Control", "RNase", "mark")
# df = melt(df)
# 
# pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_2e.pdf", height = 1.7, width = 1.6)
# ggplot(df, aes(x=variable, y=value)) +
#   stat_boxplot(geom = "errorbar", width = 0.3, size = 0.2) +
#   geom_boxplot(outlier.shape = NA, lwd = 0.3) +
#   labs(x="", y="") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=6),
#         axis.text.x = element_text(size = 6),
#         axis.text.y = element_text(size = 6),
#         legend.position = "none")
# dev.off()
  

####### Figure 2h (GSE98671_CTCF_degron.r)
# df = rbind(data.frame(sample="AID disrupted", RAL = CTCF_union_mm10[which(CTCF_union_mm10$untreated == TRUE & CTCF_union_mm10$auxin2days == FALSE), "RAL"]),
#            data.frame(sample="Retained", RAL = CTCF_union_mm10[which(CTCF_union_mm10$untreated == TRUE & CTCF_union_mm10$auxin2days == TRUE), "RAL"]))
# 
# pdf("/mnt/extraids/OceanStor-1/rcalandrelli/phase_separation/paper_plots_revision/Figure_2h.pdf", width = 1.9, height = 1.9)
# ggplot(df, aes(x=sample, y=log(RAL+1))) +
#   stat_boxplot(geom = "errorbar", width = 0.3, size = 0.2) +
#   geom_boxplot(outlier.shape = NA, lwd = 0.3) +
#   labs(x="", y="Log(RAL+1)") +
#   ylim(c(0,4.5)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=8),
#         axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8),
#         legend.position = "none")
# dev.off()


####### Figure 2g
# df = data.frame(Gr_union_tad_boundaries_FULL)
# 
# label_size = 7
# 
# df1 = df[which(df$control == 1),]
# n = nrow(df1)
# df1 = df1[which(!is.na(rowSums(df1[,c("insulation_score_control","insulation_score_RNase")]))),]
# df1 = df1[which(!is.infinite(rowSums(df1[,c("insulation_score_control","insulation_score_RNase")]))),]
# df1$color = factor(ifelse(df1$control == 1 & df1$RNase == 0, "Only control", "Other"))
# df1$delta_insulation = 2^df1$insulation_score_RNase - 2^df1$insulation_score_control
# df1$delta_RAL = df1$RAL_RNase - df1$RAL_control
# df1$log2_delta_RAL = -log2(-df1$delta_RAL)
# 
# intervals = c(min(df1$delta_RAL), rev(-seq(0, 100000, 10000)))
# 
# breaks = seq(1:(length(intervals)-1))
# labels = c()
# for (j in 1:(length(intervals)-1)){
#   labels = c(labels, paste0("[",round(intervals[j]/1000),":",round(intervals[j+1]/1000),")"))
# }
# 
# out_list = list()
# for (j in 1:(length(intervals)-1)){
#   df_temp = df1[which(df1$delta_RAL >= intervals[j] & df1$delta_RAL < intervals[j+1]),]
#   out = data.frame(x = c(j,j),
#                    type = c("Positive delta insulation", "Negative delta insulation"),
#                    delta_insulation_proportion = c(sum(df_temp$delta_insulation > 0)/nrow(df_temp), 1-sum(df_temp$delta_insulation > 0)/nrow(df_temp)))
#   
#   out$x_lab = paste0("[",round(intervals[j]/1000),":",round(intervals[j+1]/1000),")")
#   out_list[[j]] = out
# }
# 
# df = do.call("rbind", out_list)
# df$x_lab = factor(df$x_lab, levels = labels)
# df$type = factor(df$type)
# 
# pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_2g.pdf", width = 2, height = 2)
# ggplot(df, aes(x=x,y=delta_insulation_proportion,fill=type))+
#   geom_bar(stat = "identity", color="white") +
#   labs(x = "", y = "Proportion of boundaries") +
#   scale_x_continuous(name="Delta RAL", breaks = breaks, labels = labels) +
#   scale_fill_manual(values = c("#00b6bc","#f76c62")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=6),
#         axis.text.x = element_text(size=5, angle = 45, hjust = 1),
#         axis.text.y = element_text(size=6),
#         legend.title = element_blank(),
#         legend.position = "none")
# dev.off()
# 
# 
# out_list = list()
# for (j in 1:(length(intervals)-1)){
#   df_temp = df1[which(df1$delta_RAL >= intervals[j] & df1$delta_RAL < intervals[j+1]),]
#   out = data.frame(delta_insulation = c(sum(df_temp$delta_insulation > 0), sum(df_temp$delta_insulation <= 0)))
#   #out = data.frame(delta_insulation = c(sum(df_temp$delta_insulation > 0)/nrow(df_temp), 1-sum(df_temp$delta_insulation > 0)/nrow(df_temp)))
#   out_list[[j]] = out
# }
# 
# temp = data.frame(do.call("cbind", out_list))
# chisq.test(as.table(as.matrix(temp)))


##### Supplementary Figure S-Delta b
df = data.frame(Gr_union_tad_boundaries_FULL)

label_size = 7

df1 = df[which(df$control == 1),]
n = nrow(df1)
df1 = df1[which(!is.na(rowSums(df1[,c("insulation_score_control","insulation_score_RNase")]))),]
df1 = df1[which(!is.infinite(rowSums(df1[,c("insulation_score_control","insulation_score_RNase")]))),]
df1$color = factor(ifelse(df1$control == 1 & df1$RNase == 0, "Only control", "Other"))
df1$delta_insulation = 2^df1$insulation_score_RNase - 2^df1$insulation_score_control
df1$delta_RAL = df1$RAL_RNase - df1$RAL_control
df1$log2_delta_RAL = -log2(-df1$delta_RAL)

intervals = c(min(df1$delta_RAL), rev(-seq(0, 100000, 10000)))
#intervals = round(-log2(-c(min(df1$delta_RAL), rev(-seq(1, 100000, 10000)))),2)
#intervals = rev(-seq(6, 18))

out_list = list()
for (j in 1:(length(intervals)-1)){
  df_temp = df1[which(df1$delta_RAL >= intervals[j] & df1$delta_RAL < intervals[j+1]),]
  out = data.frame(x = j,
                   mean = mean(df_temp$delta_insulation),
                   sd = sd(df_temp$delta_insulation),
                   n = nrow(df_temp))
  out$se = out$sd / sqrt(out$n)
  out$x_lab = paste0("[",intervals[j],":",intervals[j+1],")")
  out_list[[j]] = out
}

df = do.call("rbind", out_list)

breaks = seq(1:(length(intervals)-1))
labels = c()
for (j in 1:(length(intervals)-1)){
  labels = c(labels, paste0("[",round(intervals[j]/1000),":",round(intervals[j+1]/1000),")"))
}

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_S-delta_b.pdf", width = 3.5, height = 3.5)
ggplot(df, aes(y=mean, x=x)) +
  geom_point(size = 2) +
  labs(x = "Delta RAL", y = "Average delta insulation") +
  scale_x_continuous(name="Delta RAL", breaks = breaks, labels = labels) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size=7, angle = 45, hjust = 1),
        axis.text.y = element_text(size=7),
        legend.title = element_blank(),
        legend.position = "top")
dev.off()


for (j in 1:(length(intervals)-1)){
  df_temp = df1[which(df1$delta_RAL >= intervals[j] & df1$delta_RAL < intervals[j+1]),]
  out = data.frame(x = j,
                   delta_RAL = df_temp$delta_RAL)
  out_list[[j]] = out
}

temp = do.call("rbind", out_list)
temp$x = factor(temp$x)
var.test(delta_RAL ~ x, data = temp)
summary(aov(delta_RAL ~ x, data = temp))

  
####### Figure 2h (GSE98671_CTCF_degron.r)
# i = "control"
# df = data.frame(Gr_union_tad_boundaries_CTCF_degron_FULL_norm)
# df = df[which(df$width <= 200000),]
# df = df[which(!is.infinite(rowSums(df[,paste0("insulation_score_",c("control","treated"))])) &
#                 !is.na(rowSums(df[,paste0("insulation_score_",c("control","treated"))]))),]
# df$persistent = 0
# df[which(df$control+df$treated == 2),"persistent"] = 1
# df$persistent = factor(df$persistent)
# df$alpha = 0.5
# 
# df = df[which(df[,i] == 1),]
#   
# pdf("/mnt/extraids/OceanStor-1/rcalandrelli/phase_separation/paper_plots_revision/Figure_2h.pdf", width = 2, height = 1.9)
# ggplot(data=df, aes(x=RAL_mESC)) +
#   geom_histogram(aes(fill=persistent), color ="black", bins = 30, alpha = 0.4, position = "identity", size = 0.2) +
#   geom_density(aes(fill=persistent), alpha=0.6, size = 0.2) +
#   xlab("RAL") +
#   ylab("Boundary numbers") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=8),
#         axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8),
#         legend.position = "none")
# dev.off()




########## Figure 1a old alternative
# plot_compartment_imargi_tracks <- function(eigen_file,
#                                            input_bedpe_file,
#                                            sample_name,
#                                            chromosomes,
#                                            resolution,
#                                            start_coord=0,
#                                            end_coord=0,
#                                            tick_dist=50000000,
#                                            minor_tick_dist=10000000,
#                                            log_margi=F,
#                                            MALAT1_cov=F,
#                                            y_max_RNA=NA,
#                                            y_max_DNA=NA,
#                                            figure) # EIGEN or IF
# {
#   
#   ### Load gene density data
#   gene_density = read.table(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/gene_density_annotation/gene_density_",resolution,".txt"), stringsAsFactors = F)
#   
#   ### Load compartment data
#   compartment = read.table(eigen_file, header = T)
#   compartment = cbind(compartment,compartment[,3])
#   compartment[,3][which(compartment[,3]<0)] = NA
#   compartment[,4][which(compartment[,4]>=0)] = NA
#   colnames(compartment) = c("chr", "coord", "eigen_pos", "eigen_neg")
#   #compartment = compartment[which(compartment$chr != "chrY"),]
#   
#   ### Load iMARGI data
#   # input_bedpe = data.frame(fread(input_bedpe_file))
#   # input_data = input_bedpe[which(input_bedpe[,1] == "chr11" | 
#   #                                  input_bedpe[,4] == "chr11"),1:6]
#   # input_data = input_data[which(input_data[,1] %in% chromosomes & 
#   #                                 input_data[,4] %in% chromosomes),]
#   # rm(input_bedpe)
#   # input_data[,2] = input_data[,2] + 1 # 1-based 
#   # input_data[,5] = input_data[,5] + 1 # 1-based 
#   # input_data[,1] = as.character(input_data[,1])
#   # input_data[,4] = as.character(input_data[,4])
#   # colnames(input_data) = paste0("V",seq(1,6))
#   # 
#   # Gr_RNA_end = GRanges(
#   #   seqnames = Rle(as.character(input_data[,1])),
#   #   ranges = IRanges(input_data[,2], end = input_data[,3], names = c(1:nrow(input_data))),
#   #   strand = Rle(strand('*')))
#   # 
#   # Gr_DNA_end = GRanges(
#   #   seqnames = Rle(as.character(input_data[,4])),
#   #   ranges = IRanges(input_data[,5], end = input_data[,6], names = c(1:nrow(input_data))),
#   #   strand = Rle(strand('*')))
#   # 
#   # rm(input_data)
#   # 
#   # if (MALAT1_cov == T){
#   #   malat1 = annotation_hg38[which(annotation_hg38$gene_name=="MALAT1"),][1,]
#   #   Gr_malat1 = GRanges(
#   #     seqnames = Rle(as.character(malat1[1])),
#   #     ranges = IRanges(as.numeric(malat1[2]), end = as.numeric(malat1[3]), names = 1),
#   #     strand = Rle(strand('*')))
#   #   
#   #   overlap = countOverlaps(Gr_RNA_end, Gr_malat1, ignore.strand = T)
#   #   mcols(Gr_RNA_end)["overlap_malat1"] = overlap
#   #   Gr_RNA_end = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap_malat1"] >= 1]
#   #   
#   #   Gr_DNA_end = Gr_DNA_end[names(Gr_RNA_end)]
#   # }
#   
#   # cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
#   #                                      genome_gr = Gr_hg38, 
#   #                                      window_size = resolution,
#   #                                      amplifier = c(1,1), 
#   #                                      threshold = NA, 
#   #                                      names = c("RNA_end","DNA_end"))
#   
#   cov_data = read.table(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/read_coverage_binned/200k_filter/iMARGI_H1_control/", resolution, ".txt"), header = T)
#   cov_data$log_RNA_end = log(cov_data$RNA_end)
#   cov_data$log_RNA_end[which(is.infinite(cov_data$log_RNA_end))] = 0
#   
#   
#   ##### SPIN states
#   ma_H1 <- fread("/mnt/extraids/SDSC_NFS/wenxingzhao/database/4DN/SPIN/H1_new.SPIN.JAWG.25kb.9_state.bed")
#   STATES <- c("Speckle",  "Interior_Act1","Interior_Act2","Interior_Act3","Interior_Repr1", "Interior_Repr2","Near_Lm1","Near_Lm2","Lamina")
#   color_ <-  c("#8b254a","#c14e4c","#ec7b57","#f2b579","#dbd291","#a8d29f", "#5fbba2","#7d9a98","#54508b")
#   H1_spin <- ma_H1 %>% dplyr::filter(!grepl("NAN",V4)) %>% dplyr::mutate(V2=V2+1) %>% 
#     mutate(spinidx = seq(1,nrow(.))) %>% 
#     GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T,
#                                             seqnames.field = "V1",
#                                             start.field = "V2",
#                                             end.field = "V3")
#   H1_spin$color <- color_[match(H1_spin$V4,STATES)]
#   
#   
#   for (my_chr in c("chr11")){
#     
#     if (figure == "EIGEN"){
#       print(paste0(sample_name, " - ", my_chr))
#       start_coord = 0
#       end_coord = hg38_lengths[my_chr]
#       detail.region <- toGRanges(data.frame(my_chr, start_coord, end_coord))
#       
#       y_text = 0.1 # track label height
#       
#       line_width = 1 # width of the line track
#       label_cex = 0.7
#       
#       # DNA_end area
#       r1_DNA_end = 1
#       r0_DNA_end = r1_DNA_end - 0.2
#       if (is.na(y_max_DNA)) {
#         y_max_DNA_end = round(max(cov_data[which(cov_data$seqnames==my_chr),"DNA_end"]))
#       } else {
#         y_max_DNA_end = y_max_DNA
#       }
#       
#       # gene density area
#       r1_gene_density = r0_DNA_end - 0.12
#       r0_gene_density = r1_gene_density - 0.2
#       y_max_gene_density = round(max(gene_density[which(gene_density$V1==my_chr),"V4"]))
#       
#       # Compartment area
#       r1_comp = r0_gene_density - 0.12
#       r0_comp = r1_comp - 0.2
#       ymin_comp = round(min(compartment[which(compartment$chr==my_chr),4][which(!is.na(compartment[which(compartment$chr==my_chr),4]))]),1)
#       ymax_comp = round(max(compartment[which(compartment$chr==my_chr),3][which(!is.na(compartment[which(compartment$chr==my_chr),3]))]),1)
#       
#       # SPIN states area
#       r1_spin = r0_comp - 0.12
#       r0_spin = r1_spin - 0.08
#       
#       pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_EIGEN_A.pdf"), width = 4, height = 2.8)
#       plot_params <- getDefaultPlotParams(plot.type=1)
#       plot_params$data1inmargin <- 5
#       plot_params$ideogramheight <- 10
#       plot_params$bottommargin <- 30
#       plot_params$topmargin <- 0
#       plot_params$leftmargin <- 0.12
#       plot_params$rightmargin <- 0.02
#       kp <- plotKaryotype(genome=Gr_hg38, plot.type=1, plot.params = plot_params, zoom = detail.region, cex=0.8)
#       kpAddBaseNumbers(kp, tick.dist = tick_dist, tick.len = 6, tick.col="#4d4d4d", cex=0.7,
#                        minor.tick.dist = minor_tick_dist, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
#       #kpAddMainTitle(kp, main="test", col="red")
#       
#       # DNA end track
#       kpDataBackground(kp, data.panel = 1, r0=r0_DNA_end, r1=r1_DNA_end, color = "#ffffff")
#       kpAxis(kp, ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, col="gray50", cex=label_cex, numticks = 3)
#       # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
#       #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
#       kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$DNA_end,
#               col="purple", ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, lwd=line_width)
#       
#       # Compartment barplot
#       kpDataBackground(kp, data.panel = 1, r0=r0_comp, r1=r1_comp, color = "#ffffff")
#       kpAxis(kp, ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, col="gray50", cex=label_cex, numticks = 3)
#       # kpText(kp, chr=compartment$chr, x=mean(c(start(detail.region),end(detail.region))), 
#       #        y=y_text, col="#000000", r0=r1_comp + 0.04, r1=r1_comp + 0.06, labels="A/B compartments", cex=label_cex)
#       kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_pos,
#              col="#E41A1C", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, border = "#E41A1C")
#       kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_neg,
#              col="#0cad01", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, border = "#0cad01")
#       
#       # Gene density track
#       kpDataBackground(kp, data.panel = 1, r0=r0_gene_density, r1=r1_gene_density, color = "#ffffff")
#       kpAxis(kp, ymin=0, ymax=y_max_gene_density, r0=r0_gene_density, r1=r1_gene_density, col="gray50", cex=label_cex, numticks = 3)
#       # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
#       #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
#       kpLines(kp, chr=gene_density$V1, x=rowMeans(gene_density[,2:3]), y=gene_density$V4,
#               col="blue", ymin=0, ymax=y_max_gene_density, r0=r0_gene_density, r1=r1_gene_density, lwd=line_width)
#       
#       ### SPIN states
#       # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
#       #        y=y_text, col="#000000", r0=r1_spin + 0.04, r1=r1_spin + 0.06, labels="SPIN States", cex=label_cex)
#       kpPlotRegions(kp, data=H1_spin, col=H1_spin$color, r0=r0_spin , r1=r1_spin, data.panel=1)
#       
#       
#       dev.off()
#       
#     }
#     
#     else if (track_orientation == "IF"){
#       
#       print(paste0(sample_name, " - ", my_chr))
#       start_coord = 0
#       end_coord = hg38_lengths[my_chr]
#       detail.region <- toGRanges(data.frame(my_chr, start_coord, end_coord))
#       
#       y_text = 0.1 # track label height
#       
#       line_width = 0.6 # width of the line track
#       label_cex = 0.4
#       
#       # Compartment area
#       r1_comp = 0.95
#       r0_comp = r1_comp - 0.3
#       ymin_comp = round(min(compartment[which(compartment$chr==my_chr),4][which(!is.na(compartment[which(compartment$chr==my_chr),4]))]),1)
#       ymax_comp = round(max(compartment[which(compartment$chr==my_chr),3][which(!is.na(compartment[which(compartment$chr==my_chr),3]))]),1)
#       
#       # DNA_end area
#       r1_DNA_end = r0_comp - 0.2
#       r0_DNA_end = r1_DNA_end - 0.3
#       if (is.na(y_max_DNA)) {
#         y_max_DNA_end = round(max(cov_data[which(cov_data$seqnames==my_chr),"DNA_end"]))
#       } else {
#         y_max_DNA_end = y_max_DNA
#       }
#       
#       plot_params <- getDefaultPlotParams(plot.type=1)
#       plot_params$data1inmargin <- 5
#       plot_params$ideogramheight <- 10
#       plot_params$bottommargin <- 30
#       plot_params$topmargin <- 0
#       plot_params$leftmargin <- 0.02
#       plot_params$rightmargin <- 0.02
#       
#       png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_",sample_name,".png"), width = 2.8, height = 1, units = "in", res = 300)
#       kp <- plotKaryotype(genome=Gr_hg38, plot.type=1, plot.params = plot_params, zoom = detail.region, cex=0.8)
#       kpAddBaseNumbers(kp, tick.dist = tick_dist, tick.len = 6, tick.col="#4d4d4d", cex=0.7,
#                        minor.tick.dist = minor_tick_dist, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
#       #kpAddMainTitle(kp, main="test", col="red")
#       
#       # Compartment barplot
#       kpDataBackground(kp, data.panel = 1, r0=r0_comp, r1=r1_comp, color = "#ffffff")
#       kpAxis(kp, ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, col="gray50", cex=label_cex, numticks = 3)
#       # kpText(kp, chr=compartment$chr, x=mean(c(start(detail.region),end(detail.region))), 
#       #        y=y_text, col="#000000", r0=r1_comp + 0.04, r1=r1_comp + 0.06, labels="A/B compartments", cex=label_cex)
#       kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_pos,
#              col="#E41A1C", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, border = "#E41A1C")
#       kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_neg,
#              col="#0cad01", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, border = "#0cad01")
#       
#       # DNA end track
#       kpDataBackground(kp, data.panel = 1, r0=r0_DNA_end, r1=r1_DNA_end, color = "#ffffff")
#       kpAxis(kp, ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, col="gray50", cex=label_cex, numticks = 3)
#       # kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
#       #        y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
#       kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$DNA_end,
#               col="purple", ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, lwd=line_width)
#       
#       dev.off()
#       
#     }
#     
#   }
# }
# 
# ############### Figure EIGEN-A
# sample = all_imargi_samples[1]
# sample_name=gsub("iMARGI_","",sample)
# temp = strsplit(sample,"_")[[1]][3]
# eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_control_merged/H1_", temp, "_500000.txt")
# input_bedpe_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz")
# sample_name=sample_name
# 
# plot_compartment_imargi_tracks(eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_", temp, "_500000.txt"),
#                                input_bedpe_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz"),
#                                sample_name=sample_name,
#                                chromosomes=hg38_chromosomes[1:23],
#                                resolution=500000,
#                                tick_dist=20000000,
#                                minor_tick_dist=5000000,
#                                log_margi = F,
#                                figure = "EIGEN")


######## Figure 3i
# index_to_plot = "P2LL"
# percentile = 0.99
# 
# df = df_upset_H1_loops_full[which(df_upset_H1_loops_full$end2 - df_upset_H1_loops_full$start1 > 0),]
# 
# df1 = df[,paste0(index_to_plot, "_",gsub("H1_","",c("control", "NH4OAc", "FL", "RNase")))]
# df1 = melt(df1)
# df1$variable = gsub(paste0(index_to_plot, "_"), "", df1$variable)
# df1$variable = factor(df1$variable)
# df1 = df1[which(!is.na(df1$value)),]
# df1 = df1[which(df1$value < quantile(df1$value, percentile)),]
# df1$index = "union"
# 
# 
# out_list = list()
# for (j in c("control", "NH4OAc", "FL", "RNase")){
#   temp = data.frame(value = df[which(df[,gsub("H1_","",j)] == 1), paste0(index_to_plot, "_", gsub("H1_","",j))])
#   temp$sample = gsub("H1_","",j)
#   out_list[[j]] = temp
# }
# df2 = data.frame(do.call("rbind", out_list))
# df2$sample = factor(df2$sample)
# df2 = df2[which(df2$value < quantile(df2$value, percentile)),]
# df2$index = "sample"
# df2 = df2[,c(2,1,3)]
# colnames(df2) = c("variable","value","index")
# 
# df_plot = rbind(df1, df2)
# df_plot$variable = factor(df_plot$variable, levels = c("control", "NH4OAc", "FL", "RNase"))
# 
# pdf(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure3i.pdf"), width = 1.52, height = 1.7)
# ggplot(data=df_plot, aes(y=value, x=variable, fill=index)) +
#   stat_boxplot(geom ='errorbar', width = 0.4, size = 0.25, position=position_dodge(0.9)) +
#   geom_boxplot(outlier.shape = NA, lwd = 0.25, position=position_dodge(0.9)) +
#   #stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.x = 0.5, label.y = max(df$value) * 1.05, size = 3, method = "wilcox.test") +
#   xlab("") +
#   ylab(index_to_plot) +
#   ylim(c(0, max(df_plot$value) * 0.9)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=6),
#         axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 6),
#         legend.position = "none")
# dev.off()
  

####### Figure 3a
make_annotation_tracks(chr="chr19", start=55080000, end=55230000, sample_name="H1_MicroC") #AAVS1

####### Figure 3c
library(UpSetR)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_3c.pdf", width = 7, height = 5)
upset(df_upset_H1_loops_10000, 
      sets = c("RNase","FL","NH4OAc","Control"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      #text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
      text.scale = c(1.8, 1.8, 0, 0, 1.6, 1.7),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 3, line.size = 0.5)
dev.off()

#### Figure 3d
all_HiC_samples = paste0("H1_", c("control","FL","NH4OAc","RNase"))

for (my_sample in all_HiC_samples){
  for (replicate in c(1,2)){
    
    my_cutoffs = c(0,30)
    my_bin_size = 10000
    
    p1 <- make_complex_heatmap_HiC(
      chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
      sample=paste0("HiC_", my_sample),
      contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_", replicate, "/", my_bin_size,"/chr12_chr12_",my_bin_size,".txt"),
      sparse=T,
      chr_row="chr12",
      chr_col="chr12",
      bin_size=my_bin_size,
      select_coord=rep(c(52850000,53650000),2),
      contacts_cutoffs=my_cutoffs,
      my_colorbar=c("white","red"),
      Gr_row_annotation=NA,
      Gr_col_annotation=NA,
      Gr_stripes_row=NA,
      Gr_stripes_col=NA,
      compartments_file_path= NA,
      tads_file_path = NA,
      loops_file_path = NA,
      feature_color = "blue")
    
    png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_3d_",my_sample,"_",replicate,".png"), width = 5, height = 5, res = 300, units = "in")
    print(p1)
    dev.off()
  }
}



#### Figure 3e
df1 = summary_apa_5000[which(summary_apa_5000$Index == "P2LL"),]
df1$Index = "Per condition"

df2 = summary_apa_5000_union[which(summary_apa_5000_union$Index == "P2LL"),]
df2$Index = "Union"

df = rbind(df1,df2)
df$variable = factor(df$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_3e.pdf", width = 2, height = 2)
ggplot(df, aes(fill=Index, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "P2LL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "none")
dev.off()

###### Figure 3f
df1 = summary_apa_5000[which(summary_apa_5000$Index == "ZscoreLL"),]
df1$Index = "Per condition"

df2 = summary_apa_5000_union[which(summary_apa_5000_union$Index == "ZscoreLL"),]
df2$Index = "Union"

df = rbind(df1,df2)
df$variable = factor(df$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))

p2<-ggplot(df, aes(fill=Index, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "ZscoreLL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "none")


##### Figure 3g (/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/loop_analysis.r)

# df_10000 = melt(df_p2m_union_10000)
# df_10000$variable = factor(df_10000$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))

df = rbind(reshape2::melt(df_p2m_union_10000), reshape2::melt(df_p2m_union_5000))
df$variable = factor(df$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))

my_comparisons <- list(c("Control","FL"), c("Control","RNase"))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_3g.pdf", width = 2, height = 2)
ggplot(df[which(df$value <= 20),], aes(x=variable, y=value)) +
  stat_boxplot(geom ='errorbar', width = 0.3, size = 0.25, position=position_dodge(0.9)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, position=position_dodge(0.9)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = c(8,10), size = 3.5, method = "wilcox") +
  labs(x="", y="P2M") +
  ylim(0,12) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))
dev.off()


#### Figure 4e
plot_list = list()
for (my_sample in all_imargi_samples[c(1,2,4)]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  if (temp == "RNase"){
    loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_FL_merged") # I want to plot the same emergent loop
  }
  
  my_cutoffs = c(0,60)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
    sparse=T,
    chr_row="chr20",
    chr_col="chr20",
    bin_size=10000,
    select_coord=rep(c(47200000,47400000),2),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = NA, #loop_file_path_temp,
    feature_color = "blue")
  
  plot_list[[my_sample]] = p1
}

plot_list_grob = lapply(plot_list, as.grob)

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_4e.pdf", width = 7.5, height = 2.5)
plot_grid(plotlist=plot_list_grob, ncol = 3)
dev.off()


#### Figure 4f
# temp = data.frame(sample = c(rep("Scramble Ctrl", 3), rep("ZMYND8 gRNA",3)),
#                   mean = c(1.053, 0.223),
#                   std = c(0.376, 0.076))
# 
# pdf("/mnt/extraids/OceanStor-1/rcalandrelli/phase_separation/paper_plots_revision/Figure_4f.pdf", width = 1.6, height = 1.7)
# ggplot(data = temp, aes(x=sample, y=mean)) +
#   geom_bar(stat="identity") +
#   xlab("") +
#   ylab("Normalized ZMYND8 expression level") +
#   geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2, size = 0.4) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=5.5),
#         axis.text.x = element_text(size=5),
#         axis.text.y = element_text(size=5.5),
#         legend.title = element_blank())
# dev.off()


df = data.frame(f = c(21.2925300598145,20.3315877914429,20.2086572647095,22.742413520813,23.403431892395,22.364182472229),
                sample = c(rep("Scramble Ctrl",3), rep("ZMYND8 gRNA",3)))
df$g = c(rep(aggregate(df$f, by = list(df$sample), FUN = mean)[1,2], 6))
df$h = df$f - df$g
df$j = 2^(-1 * df$h)

df1 = summarySE(df, measurevar="j", groupvars=c("sample"))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_4f.pdf", width = 1.8, height = 2)
ggplot() +
  geom_bar(data = df1, aes(x=sample, y=j), stat="identity") +
  geom_errorbar(data = df1, aes(x=sample, ymin=j-se, ymax=j+se), width=.2, size = 0.4) +
  geom_jitter(data=df, aes(x=sample, y=j), size = 1, width = 0.25) +
  xlab("") +
  ylab("Normalized ZMYND8 expression level") +
  ylim(0,1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.title = element_blank())
dev.off()


####### Figure 5a
plot_boxplots_loop_score_index_paper <- function(index_to_plot,
                                                 percentile){
  
  ##### All loops
  plot_list1 = list()
  my_comparisons <- list(c("control","RNase"))
  
  # plot_list1[[1]] <- ggplot(data=df, aes(y=value, x=variable)) +
  #   geom_boxplot() +
  #   stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.x = c(2,3,4), label.y = c(0,0.5,1) + max(df$value) * 1.1, size = 5, method = "wilcox.test") +
  #   xlab("") +
  #   ylab(index_to_plot) +
  #   theme_bw() +
  #   theme(text = element_text(size=12),
  #         axis.text.x = element_text(size = 12),
  #         axis.text.y = element_text(size = 12)) +
  #   ggtitle("All loops")
  
  k = 1
  for (i in c("Convergent CBS","Non-convergent CBS","No CBS")){
    
    df = df_upset_H1_loops_full # [which(df_upset_H1_loops_full$end2 - df_upset_H1_loops_full$start1 > 0),] # why?
    
    if (i == "Convergent CBS"){
      df = df[which(df$CBS_orientation == "+-"),]
    } else if (i == "Non-convergent CBS"){
      df = df[which(df$CBS_orientation %in% c("-+","--","++")),]
    } else {
      df = df[which(df$CBS_orientation == "NA"),]
    }
    
    wilcox.test(df$P2LL_control, df$P2LL_RNase)
    wilcox.test(df$P2LL_control, df$P2LL_RNase, paired = T)
    
    df = df[,paste0(index_to_plot, "_",gsub("H1_","",c("control","RNase")))]
    df = melt(df)
    df$variable = gsub(paste0(index_to_plot, "_"), "", df$variable)
    df$variable = factor(df$variable)
    df = df[which(!is.na(df$value)),]
    df = df[which(df$value < quantile(df$value, percentile)),]
    
    plot_list1[[k]] <- ggplot(data=df, aes(y=value, x=variable)) +
      stat_boxplot(geom ='errorbar', width = 0.2, size = 0.25) +
      geom_boxplot(outlier.shape = NA, lwd = 0.25) +
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.x = 0.5, label.y = max(df$value) * 1.05, size = 3, method = "wilcox.test") +
      xlab("") +
      ylab(index_to_plot) +
      ylim(c(0, max(df$value) * 1.2)) +
      theme_bw() +
      theme(text = element_text(size=8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8))
    k = k + 1
    
    
  }
  
  p <- plot_grid(plotlist = plot_list1, ncol = 3)
  pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/custom/P2LL_boxplots_all_loops_paper.pdf"), width = 4, height = 1.7)
  print(p)
  dev.off()
  
  
  ##### Loops each sample
  plot_list1 = list()
  
  k = 1
  for (i in c("Convergent CBS","Non-convergent CBS","No CBS")){
    
    df = df_upset_H1_loops_full[which(df_upset_H1_loops_full$end2 - df_upset_H1_loops_full$start1 > 100000),]
    
    if (i == "Convergent CBS"){
      df = df[which(df$CBS_orientation == "+-"),]
    } else if (i == "Non-convergent CBS"){
      df = df[which(df$CBS_orientation %in% c("-+","--","++")),]
    } else {
      df = df[which(df$CBS_orientation == "NA"),]
    }
    
    out_list = list()
    for (j in c("control","RNase")){
      temp = data.frame(value = df[which(df[,gsub("H1_","",j)] == 1), paste0(index_to_plot, "_", gsub("H1_","",j))])
      temp$sample = gsub("H1_","",j)
      out_list[[j]] = temp
    }
    df = data.frame(do.call("rbind", out_list))
    df$sample = factor(df$sample)
    df = df[which(df$value < quantile(df$value, percentile)),]
    
    my_comparisons <- list(c("control","RNase"))
    
    plot_list1[[k]] <- ggplot(data=df, aes(y=value, x=sample)) +
      stat_boxplot(geom ='errorbar', width = 0.2, size = 0.25) +
      geom_boxplot(outlier.shape = NA, lwd = 0.25) +
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.x = 0.5, label.y = max(df$value) * 1.05, size = 3, method = "wilcox.test") +
      xlab("") +
      ylab(index_to_plot) +
      ylim(c(0, max(df$value) * 1.2)) +
      theme_bw() +
      theme(text = element_text(size=8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8))
    k = k + 1
  }
  
  p <- plot_grid(plotlist = plot_list1, ncol = 3)
  pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/custom/P2LL_boxplots_loops_each_sample_100000_paper.pdf"), width = 4, height = 1.7)
  print(p)
  dev.off()
  
}


plot_boxplots_loop_score_index_paper(index_to_plot = "P2LL",
                                     percentile = 0.99)




###### Figure 5b
a = nrow(df_upset_H1_loops_full[which(df_upset_H1_loops_full$RNase == 1 & df_upset_H1_loops_full$control == 0 &
                                        df_upset_H1_loops_full$HERV_H_caRNA_associated_reads > 80),])
b = nrow(df_upset_H1_loops_full[which(df_upset_H1_loops_full$RNase == 1 & df_upset_H1_loops_full$control == 0 &
                                        df_upset_H1_loops_full$HERV_H_caRNA_associated_reads <= 80),])
c = nrow(df_upset_H1_loops_full[which(df_upset_H1_loops_full$control == 1 &
                                        df_upset_H1_loops_full$HERV_H_caRNA_associated_reads > 80),])
d = nrow(df_upset_H1_loops_full[which(df_upset_H1_loops_full$control == 1 &
                                        df_upset_H1_loops_full$HERV_H_caRNA_associated_reads <= 80),])
a*d/(b*c)
matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)
chisq.test(as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)))



df = data.frame(x = a*d/(b*c),
                err = sqrt(1/a + 1/b + 1/c + 1/d))

exp(log(a*d/(b*c)) - sqrt(1/a + 1/b + 1/c + 1/d))
exp(log(a*d/(b*c)) + sqrt(1/a + 1/b + 1/c + 1/d))


pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_5b.pdf", width = 2, height = 0.7)
ggplot(df, aes(x=x, y=1)) +
  geom_errorbar(aes(xmin=exp(log(x)-err), xmax=exp(log(x)+err)), width=.3, position=position_dodge(.9)) +
  geom_point(size=2, color="blue") +
  scale_x_continuous(trans='log', breaks = c(1:3), labels = as.character(c(1:3)), limits = c(0.8,2.5)) +
  ylim(c(0.5,1.5)) +
  # xlab("Odds ratio (log scale)") +
  # ylab("") +
  guides(y = "none") +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_bw() +
  theme(text = element_text(size=6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed")
dev.off()






############# Figure 5c
OR = c()
err = c()

thres = 55

a = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation == "+-"),])
b = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$control == 1 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation == "+-"),])
c = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation != "+-"),])
d = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$control == 1 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation != "+-"),])
a*d/(b*c)
matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)
chisq.test(as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)))

OR = c(OR, a*d/(b*c))
err = c(err, sqrt(1/a + 1/b + 1/c + 1/d))


###
a = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation == "+-"),])
b = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$control == 1 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads < thres & df_upset_H1_loops_full_hg19$CBS_orientation == "+-"),])
c = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation != "+-"),])
d = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$control == 1 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads < thres & df_upset_H1_loops_full_hg19$CBS_orientation != "+-"),])
a*d/(b*c)
matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)
chisq.test(as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)))

OR = c(OR, a*d/(b*c))
err = c(err, sqrt(1/a + 1/b + 1/c + 1/d))

###
a = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation == "+-"),])
b = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads < thres & df_upset_H1_loops_full_hg19$CBS_orientation == "+-"),])
c = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads >= thres & df_upset_H1_loops_full_hg19$CBS_orientation != "+-"),])
d = nrow(df_upset_H1_loops_full_hg19[which(df_upset_H1_loops_full_hg19$RNase == 1 & df_upset_H1_loops_full_hg19$control == 0 &
                                             df_upset_H1_loops_full_hg19$HERV_H_caRNA_associated_reads < thres & df_upset_H1_loops_full_hg19$CBS_orientation != "+-"),])
a*d/(b*c)
matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)
chisq.test(as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T)))

OR = c(OR, a*d/(b*c))
err = c(err, sqrt(1/a + 1/b + 1/c + 1/d))


df = data.frame(x = OR,
                y = c(3,2,1),
                err = err)

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_5c.pdf", width = 2, height = 1)
ggplot(df, aes(x=x, y=y)) +
  geom_errorbar(aes(xmin=exp(log(x)-err), xmax=exp(log(x)+err)), width=.3, position=position_dodge(.9)) +
  geom_point(size=2, color="blue") +
  scale_x_continuous(trans='log', breaks = c(1:3), labels = as.character(c(1:3)), limits = c(0.8,2.5)) +
  ylim(c(0.5,3.5)) +
  # xlab("Odds ratio (log scale)") +
  # ylab("") +
  guides(y = "none") +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_bw() +
  theme(text = element_text(size=6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed")
dev.off()





########### Figure S1 RAL vs A/B
figure_s1_RAL_PC1 <- function(sample,
                              eigen_file,
                              resolution){
  
  compartment = read.table(eigen_file, header = T)
  colnames(compartment) = c("chr","coord","eigen")
  #compartment = compartment[which(compartment$chr != "chrY"),]
  
  cov_data = read.table(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/read_coverage_binned/200k_filter/", sample, "/", resolution, ".txt"), header = T)
  
  df = data.frame(compartment$eigen, cov_data$RNA_end, cov_data$DNA_end)
  label_size = 10
  
  if (sample == "iMARGI_H1_control"){
    labels = c("0","25,000","50,000","75,000")
    breaks = c(0,25000,50000,75000)
  } else if (sample == "iMARGI_HFF_control"){
    labels = c("0","10,000","20,000")
    breaks = c(0,10000,20000)
  } else if (sample == "iMARGI_K562_control"){
    labels = c("0","100,000","200,000")
    breaks = c(0,100000,200000)
  }
  
  p<-ggplot(df, aes(x=cov_data.DNA_end, y=compartment.eigen)) + 
    geom_point(color = "purple", size = 0.3) +
    #xlab("DNA end read coverage") +
    xlab("RAL") +
    ylab("PC1") +
    scale_x_continuous(labels = labels, breaks = breaks) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size),
          plot.title = element_text(size = label_size))
  
  png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_S1_PC1_", sample, ".png"), height = 2.5, width = 2.5, units = "in", res = 300)
  print(p)
  dev.off()
  
  write.table(df[,c("cov_data.DNA_end","compartment.eigen")], 
              paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S1c_",gsub("iMARGI_", "", sample), ".txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  return(output)
}


for (i in c("iMARGI_H1_control", "iMARGI_HFF_control", "iMARGI_K562_control")){
  eigen_sample = gsub("iMARGI_", "", i)
  figure_s1_RAL_PC1(sample=i,
                   eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/", eigen_sample, "/", eigen_sample, "_500000.txt"),
                   resolution=500000)
}


################## Figure S1 RAL vs gene density
figure_s1_gene_density <- function(sample,
                           resolution){
  

  cov_data = read.table(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/read_coverage_binned/200k_filter/", sample, "/", resolution, ".txt"), header = T)
  gene_density = read.table(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/gene_density_annotation/gene_density_",resolution,".txt"), stringsAsFactors = F)
  
  print(paste0("RNA end: ", cor(cov_data$RNA_end, gene_density$V4, method = "spearman")))
  print(paste0("DNA end: ", cor(cov_data$DNA_end, gene_density$V4, method = "spearman")))

  
  df = data.frame(gene_density$V4, cov_data$RNA_end, cov_data$DNA_end)
  label_size = 10
  
  ### RNA end
  if (sample == "iMARGI_H1_control"){
    labels = c("0","75,000","150,000")
    breaks = c(0,75000,150000)
  } else if (sample == "iMARGI_HFF_control"){
    labels = c("0","20,000","40,000")
    breaks = c(0,20000,40000)
  } else if (sample == "iMARGI_K562_control"){
    labels = c("0","75,000","150,000")
    breaks = c(0,75000,150000)
  }
  
  df_temp = df
  df_temp = df_temp[which(df_temp$cov_data.RNA_end < quantile(df_temp$cov_data.RNA_end, 0.99)),]
  
  p<-ggplot(df_temp, aes(x=cov_data.RNA_end, y=gene_density.V4)) + 
    geom_point(color = "purple", size = 0.2) +
    xlab("caRNA") +
    ylab("Gene density") +
    scale_x_continuous(labels = labels, breaks = breaks) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size),
          plot.title = element_text(size = label_size))
  
  png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_S1_gene_density_caRNA_", sample, ".png"), height = 2.5, width = 2.5, units = "in", res = 300)
  print(p)
  dev.off()
  
  ### DNA end
  if (sample == "iMARGI_H1_control"){
    labels = c("0","25,000","50,000","75,000")
    breaks = c(0,25000,50000,75000)
  } else if (sample == "iMARGI_HFF_control"){
    labels = c("0","10,000","20,000")
    breaks = c(0,10000,20000)
  } else if (sample == "iMARGI_K562_control"){
    labels = c("0","100,000","200,000")
    breaks = c(0,100000,200000)
  }
  
  p<-ggplot(df, aes(x=cov_data.DNA_end, y=gene_density.V4)) + 
    geom_point(color = "purple", size = 0.2) +
    xlab("RAL") +
    ylab("Gene density") +
    scale_x_continuous(labels = labels, breaks = breaks) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size),
          plot.title = element_text(size = label_size))
  
  png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots_revision/Figure_S1_gene_density_RAL_", sample, ".png"), height = 2.5, width = 2.5, units = "in", res = 300)
  print(p)
  dev.off()
  
  write.table(df[,c("cov_data.DNA_end","gene_density.V4")], 
              paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S1f_",gsub("iMARGI_", "", sample), ".txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  return()
}


for (i in c("iMARGI_H1_control", "iMARGI_HFF_control", "iMARGI_K562_control")){
  figure_s1_gene_density(sample=i,
                 resolution=500000)
}














# plot_list = list()
# 
# my_coord = rep(c(156000000,160000000),2)
# my_coord_zoom = rep(c(157500000,159000000),2)
# 
# my_sample = "H1_control"
# my_cutoffs = c(0,0.995)
# my_bin_size = 10000
# 
# p1 <- make_complex_heatmap_HiC(
#   chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#   sample=paste0("HiC_", my_sample),
#   contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
#   sparse=T,
#   chr_row="chr1",
#   chr_col="chr1",
#   bin_size=my_bin_size,
#   select_coord=my_coord,
#   contacts_cutoffs=my_cutoffs,
#   my_colorbar=c("white","red"),
#   Gr_row_annotation=NA,
#   Gr_col_annotation=NA,
#   Gr_stripes_row=NA,
#   Gr_stripes_col=NA,
#   compartments_file_path= NA,
#   tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
#   loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
#   feature_color = "blue")#dbdeff")
# 
# plot_list[[1]] = p1
# 
# p2 <- make_complex_heatmap_HiC(
#   chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#   sample=paste0("HiC_", my_sample),
#   contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/10000/chr1_chr1_10000.txt"),
#   sparse=T,
#   chr_row="chr1",
#   chr_col="chr1",
#   bin_size=10000,
#   select_coord=my_coord_zoom,
#   contacts_cutoffs=my_cutoffs,
#   my_colorbar=c("white","red"),
#   Gr_row_annotation=NA,
#   Gr_col_annotation=NA,
#   Gr_stripes_row=NA,
#   Gr_stripes_col=NA,
#   compartments_file_path= NA,
#   tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
#   loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
#   feature_color = "blue")#dbdeff")
# 
# plot_list[[2]] = p2
# 
# my_sample = "iMARGI_H1_control"
# my_cutoffs = c(0,0.995)
# 
# p3 <- make_complex_heatmap_margi(
#   chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#   sample=my_sample,
#   contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/",my_bin_size,"_filter1k/chr1_chr1.txt"),
#   sparse=T,
#   chr_row="chr1",
#   chr_col="chr1",
#   bin_size=my_bin_size,
#   select_coord=my_coord,
#   contacts_cutoffs=my_cutoffs,
#   my_colorbar=c("white","red"),
#   Gr_row_annotation=NA,
#   Gr_col_annotation=NA,
#   Gr_stripes_row=NA,
#   Gr_stripes_col=NA,
#   tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label_2.txt"),#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
#   loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
#   feature_color = "blue")#dbdeff")
# 
# plot_list[[3]] = p3
# 
# p4 <- make_complex_heatmap_margi(
#   chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#   sample=my_sample,
#   contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr1_chr1.txt"),
#   sparse=T,
#   chr_row="chr1",
#   chr_col="chr1",
#   bin_size=10000,
#   select_coord=my_coord_zoom,
#   contacts_cutoffs=my_cutoffs,
#   my_colorbar=c("white","red"),
#   Gr_row_annotation=NA,
#   Gr_col_annotation=NA,
#   Gr_stripes_row=NA,
#   Gr_stripes_col=NA,
#   tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label_2.txt"),#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
#   loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
#   feature_color = "blue")#dbdeff")
# 
# plot_list[[4]] = p4
# 
# plot_list_grob = lapply(plot_list, as.grob)
# 
# png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_TAD_IJ_",my_bin_size,".png"), width = 8, height = 4, units = "in", res = 200)
# plot_grid(plot_list_grob[[1]],plot_list_grob[[3]], ncol = 2)
# dev.off()
# 
# png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_TAD_IJ_",my_bin_size,"_zoom.png"), width = 8, height = 4, units = "in", res = 200)
# plot_grid(plot_list_grob[[2]],plot_list_grob[[4]], ncol = 2)
# dev.off()


###### Figure EIGEN H-I
# plot_list = list()

# my_sample = "H1_control"
# my_cutoffs = c(0,0.995)
# my_bin_size = 10000
# 
# p1 <- make_complex_heatmap_HiC(
#   chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#   sample=paste0("HiC_", my_sample),
#   contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
#   sparse=T,
#   chr_row="chr1",
#   chr_col="chr1",
#   bin_size=my_bin_size,
#   select_coord=rep(c(2800000,4800000),2),
#   contacts_cutoffs=my_cutoffs,
#   my_colorbar=c("white","red"),
#   Gr_row_annotation=NA,
#   Gr_col_annotation=NA,
#   Gr_stripes_row=NA,
#   Gr_stripes_col=NA,
#   compartments_file_path= NA,
#   tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
#   loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
#   feature_color = "blue")#dbdeff")
# 
# plot_list[[1]] = p1
# 
# 
# my_sample = all_imargi_samples[1]
# temp = strsplit(my_sample,"_")[[1]][3]
# loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
# 
# p2 <- make_complex_heatmap_margi(
#   chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#   sample=my_sample,
#   contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr1_chr1.txt"),
#   sparse=T,
#   chr_row="chr1",
#   chr_col="chr1",
#   bin_size=10000,
#   select_coord=rep(c(2800000,4800000),2),
#   contacts_cutoffs=c(0,0.995),
#   my_colorbar=c("white","red"),
#   Gr_row_annotation=NA,
#   Gr_col_annotation=NA,
#   Gr_stripes_row=NA,
#   Gr_stripes_col=NA,
#   tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label_2.txt"),
#   loops_file_path = loop_file_path_temp,
#   EP_loops_file = NA,#paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
#   feature_color = "blue")
# 
# plot_list[[2]] = p2

plot_list = list()
my_sample = "H1_control"
my_bin_size = 10000

p3 <- make_complex_heatmap_HiC(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=paste0("HiC_", my_sample),
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr10_chr10_",my_bin_size,".txt"),
  sparse=T,
  chr_row="chr10",
  chr_col="chr10",
  bin_size=my_bin_size,
  select_coord=rep(c(15100000,16100000),2),
  contacts_cutoffs=c(0,0.97),
  log_counts = F,
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  compartments_file_path= NA,
  tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
  loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
  feature_color = "blue")#dbdeff")

plot_list[[1]] = p3


my_sample = all_imargi_samples[1]
temp = strsplit(my_sample,"_")[[1]][3]
loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")

p4 <- make_complex_heatmap_margi(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr10_chr10.txt"),
  sparse=T,
  chr_row="chr10",
  chr_col="chr10",
  bin_size=10000,
  select_coord=rep(c(15100000,16100000),2),
  contacts_cutoffs=c(0,0.995),
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label_2.txt"),
  loops_file_path = loop_file_path_temp,
  EP_loops_file = NA,#paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
  feature_color = "blue")

plot_list[[2]] = p4


plot_list_grob = lapply(plot_list, as.grob)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_EIGEN_H_I.pdf", width = 4, height = 2)
plot_grid(plotlist = plot_list_grob, ncol = 2)
dev.off()





######## Figure STRIPE A-B

plot_list = list()
k = 1
my_sample = all_imargi_samples[1]

temp = strsplit(my_sample,"_")[[1]][3]
loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")

my_cutoffs = c(0,0.995)

p1 <- make_complex_heatmap_margi(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr1_chr1.txt"),
  sparse=T,
  chr_row="chr1",
  chr_col="chr1",
  bin_size=10000,
  select_coord=c(3250000,3750000,3250000,3750000),
  contacts_cutoffs=my_cutoffs,
  my_colorbar=c("white","red"),
  Gr_row_annotation=Gr_annotations_edb_hg38,
  Gr_col_annotation=Gr_annotations_edb_hg38,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = NA,
  loops_file_path = loop_file_path_temp,
  EP_loops_file = NA,#paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
  feature_color = "blue")

plot_list[[k]] = p1
k = k + 1

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_STRIPE_A_heatmap.pdf", width = 8, height = 8)
plot_list[[1]]
dev.off()

p2 <- make_complex_heatmap_margi(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
  sparse=T,
  chr_row="chr20",
  chr_col="chr20",
  bin_size=10000,
  select_coord=c(37150000,37650000,37150000,37650000),
  contacts_cutoffs=my_cutoffs,
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = NA,
  loops_file_path = loop_file_path_temp,
  EP_loops_file = NA,#paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
  feature_color = "blue")

plot_list[[k]] = p2
k = k + 1

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_STRIPE_B_heatmap.pdf", width = 8, height = 8)
plot_list[[2]]
dev.off()


make_annotation_tracks(chr="chr1", start=3250000, end=3750000, sample_name="control", panel = "A")
make_annotation_tracks(chr="chr20", start=37150000, end=37650000, sample_name="control", panel = "B")



######## Figure S2 i-k
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gtable)
library(grid)

# H1_SON<-read.csv('/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/IF/SON_in_H1_statistics.csv')
H1_SON<-read.csv('/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/IF/SON_in_H1_statistics_email.csv')
label_size = 6
statistic_test_size = 2
jitter_size = 0.2

##### NH4OAc
H1_SON_ctr<-H1_SON %>%
  filter(treatment=='0 mM' & quality=='' & Mean>3500)

H1_SON_NH4<-H1_SON %>%
  filter(treatment=='100 mM' & quality=='' & Mean>3500)

H1_NH4_treatment<-rbind(H1_SON_ctr, H1_SON_NH4)


my_comparisons = list(c("0 mM", "100 mM"))

p1<-ggplot(H1_NH4_treatment, aes(x = treatment, y = number_of_foci, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Counts') +
  ylim(c(min(H1_NH4_treatment$number_of_foci)-10, max(H1_NH4_treatment$number_of_foci)+20)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$number_of_foci)+10, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Number of foci')

p2<-ggplot(H1_NH4_treatment, aes(x=treatment, y = integrated_volume, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab(expression(paste('Total volume (', mu, m^3, ')'))) +
  ylim(c(min(H1_NH4_treatment$integrated_volume), max(H1_NH4_treatment$integrated_volume)+2)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$integrated_volume)+1, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Integrated volume')

p3<-ggplot(H1_NH4_treatment, aes(x = treatment, y = normalized_variance, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Variance/Mean') +
  ylim(c(min(H1_NH4_treatment$normalized_variance)-10, max(H1_NH4_treatment$normalized_variance)+600)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$normalized_variance)+300, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Coefficient of dispersion')

p4<-ggplot(H1_NH4_treatment, aes(x = treatment, y = Mean, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Fluorescence (a.u.)') +
  ylim(c(min(H1_NH4_treatment$Mean)-100, max(H1_NH4_treatment$Mean)+800)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$Mean)+400, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Mean nuclear fluorescence')

# png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_H_J.png", width = 8, height = 3, res = 200, units = "in")
# plot_grid(p1,p2,p3, ncol=3)
# dev.off()

df = H1_NH4_treatment[,c("treatment","number_of_foci","integrated_volume","normalized_variance","Mean")]
df$treatment = ifelse(df$treatment=="0 mM","control","NH4OAc")
colnames(df) = c("treatment","number_of_foci","total_volume","normalized_variance","mean_nuclear_fluorescence")
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S2_i.txt", row.names = F, col.names = T, quote = F, sep = "\t")

##### FL
H1_SON_DMSO<-H1_SON %>%
  filter(treatment=='0 uM' & quality=='')

H1_SON_FL<-H1_SON %>%
  filter(treatment=='1 uM' & quality=='')

H1_FL_treatment<-rbind(H1_SON_DMSO, H1_SON_FL)

my_comparisons = list(c("0 uM", "1 uM"))

p5<-ggplot(H1_FL_treatment, aes(x = treatment, y = number_of_foci, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Counts') +
  ylim(c(min(H1_FL_treatment$number_of_foci)-20, max(H1_FL_treatment$number_of_foci)+30)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$number_of_foci)+15, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Number of foci')

p6<-ggplot(H1_FL_treatment, aes(x=treatment, y = integrated_volume, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab(expression(paste('Total volume (', mu, m^3, ')'))) +
  ylim(c(min(H1_FL_treatment$integrated_volume), max(H1_FL_treatment$integrated_volume)+4)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$integrated_volume)+2, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Integrated volume')

p7<-ggplot(H1_FL_treatment, aes(x = treatment, y = normalized_variance, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Variance/Mean') +
  ylim(c(min(H1_FL_treatment$normalized_variance)-10, max(H1_FL_treatment$normalized_variance)+1200)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$normalized_variance)+700, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Coefficient of dispersion')

p8<-ggplot(H1_FL_treatment, aes(x = treatment, y = Mean, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Fluorescence (a.u.)') +
  ylim(c(min(H1_FL_treatment$Mean)-100, max(H1_FL_treatment$Mean)+1600)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$Mean)+1000, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Mean nuclear fluorescence')

# png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_E_G.png", width = 8, height = 3, res = 200, units = "in")
# plot_grid(p1,p2,p3, ncol=3)
# dev.off()

df = H1_FL_treatment[,c("treatment","number_of_foci","integrated_volume","normalized_variance","Mean")]
df$treatment = ifelse(df$treatment=="0 uM","control","FL")
colnames(df) = c("treatment","number_of_foci","total_volume","normalized_variance","mean_nuclear_fluorescence")
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S2_j.txt", row.names = F, col.names = T, quote = F, sep = "\t")




##### RNase
H1_RNase_ctr_SON<-H1_SON %>%
  filter(treatment=='0 ug/mL' & Mean>7500 & Mean<10000)

H1_RNase_treat_SON<-H1_SON %>%
  filter(treatment=='200 ug/mL'  & Mean>6800 & Mean<10000)

H1_RNase_treatment<-rbind(H1_RNase_ctr_SON, H1_RNase_treat_SON)

my_comparisons= list(c("0 ug/mL", "200 ug/mL"))

p9<-ggplot(H1_RNase_treatment, aes(x = treatment, y = number_of_foci, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Counts') +
  ylim(c(min(H1_RNase_treatment$number_of_foci)-20, max(H1_RNase_treatment$number_of_foci)+25)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$number_of_foci)+15, size=statistic_test_size) +
  #stat_compare_means(method='wilcox.test') +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Number of foci')

p10<-ggplot(H1_RNase_treatment, aes(x=treatment, y = integrated_volume, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab(expression(paste('Total volume (', mu, m^3, ')'))) +
  ylim(c(min(H1_RNase_treatment$integrated_volume)-1, max(H1_RNase_treatment$integrated_volume)+3)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$integrated_volume)+2, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Integrated volume')

p11<-ggplot(H1_RNase_treatment, aes(x = treatment, y = normalized_variance, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Variance/Mean') +
  ylim(c(min(H1_RNase_treatment$normalized_variance)-10, max(H1_RNase_treatment$normalized_variance)+800)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$normalized_variance)+500, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Coefficient of dispersion')


p12<-ggplot(H1_RNase_treatment, aes(x = treatment, y = Mean, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Fluorescence (a.u.)') +
  ylim(c(min(H1_RNase_treatment$Mean)-100, max(H1_RNase_treatment$Mean)+500)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$Mean)+300, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Mean nuclear fluorescence')


df = H1_RNase_treatment[,c("treatment","number_of_foci","integrated_volume","normalized_variance","Mean")]
df$treatment = ifelse(df$treatment=="0 ug/mL","control","RNase")
colnames(df) = c("treatment","number_of_foci","total_volume","normalized_variance","mean_nuclear_fluorescence")
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S2_k.txt", row.names = F, col.names = T, quote = F, sep = "\t")



png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_IF_SON.png", width = 4.5, height = 6, res = 300, units = "in")
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol=3, byrow = F)
dev.off()


######## Figure S2 l-m
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gtable)
library(grid)

H1_SC35<-read.csv('/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/IF/SC35_in_H1_statistics_email.csv')

##### NH4OAc
H1_SC35_ctr<-H1_SC35 %>%
  filter(treatment=='control')

H1_SC35_NH4<-H1_SC35 %>%
  filter(treatment=='NH4OAc')

H1_NH4_treatment<-rbind(H1_SC35_ctr, H1_SC35_NH4)


my_comparisons = list(c("control", "NH4OAc"))

p1<-ggplot(H1_NH4_treatment, aes(x = treatment, y = number_of_nuclear_speckles, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Counts') +
  ylim(c(min(H1_NH4_treatment$number_of_nuclear_speckles)-10, max(H1_NH4_treatment$number_of_nuclear_speckles)+20)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$number_of_nuclear_speckles)+10, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Number of foci')

p2<-ggplot(H1_NH4_treatment, aes(x=treatment, y = integrated_volume, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab(expression(paste('Total volume (', mu, m^3, ')'))) +
  ylim(c(min(H1_NH4_treatment$integrated_volume)-1, max(H1_NH4_treatment$integrated_volume)+4)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$integrated_volume)+2, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Integrated volume')

p3<-ggplot(H1_NH4_treatment, aes(x = treatment, y = normalized_variance, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Variance/Mean') +
  ylim(c(min(H1_NH4_treatment$normalized_variance)-10, max(H1_NH4_treatment$normalized_variance)+900)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$normalized_variance)+500, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Coefficient of dispersion')

p4<-ggplot(H1_NH4_treatment, aes(x = treatment, y = Mean, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Fluorescence (a.u.)') +
  ylim(c(min(H1_NH4_treatment$Mean)-100, max(H1_NH4_treatment$Mean)+1000)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_NH4_treatment$Mean)+600, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Mean nuclear fluorescence')

df = H1_NH4_treatment[,c("treatment","number_of_nuclear_speckles","integrated_volume","normalized_variance","Mean")]
colnames(df) = c("treatment","number_of_foci","total_volume","normalized_variance","mean_nuclear_fluorescence")
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S2_l.txt", row.names = F, col.names = T, quote = F, sep = "\t")


# png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_H_J.png", width = 8, height = 3, res = 200, units = "in")
# plot_grid(p1,p2,p3, ncol=3)
# dev.off()

##### FL
H1_SC35_DMSO<-H1_SC35 %>%
  filter(treatment=='DMSO')

H1_SC35_FL<-H1_SC35 %>%
  filter(treatment=='FL')

H1_FL_treatment<-rbind(H1_SC35_DMSO, H1_SC35_FL)

my_comparisons = list(c("DMSO", "FL"))

p5<-ggplot(H1_FL_treatment, aes(x = treatment, y = number_of_nuclear_speckles, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Counts') +
  ylim(c(min(H1_FL_treatment$number_of_nuclear_speckles)-20, max(H1_FL_treatment$number_of_nuclear_speckles)+30)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$number_of_nuclear_speckles)+15, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Number of foci')

p6<-ggplot(H1_FL_treatment, aes(x=treatment, y = integrated_volume, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab(expression(paste('Total volume (', mu, m^3, ')'))) +
  ylim(c(min(H1_FL_treatment$integrated_volume)-1, max(H1_FL_treatment$integrated_volume)+4)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$integrated_volume)+2, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Integrated volume')

p7<-ggplot(H1_FL_treatment, aes(x = treatment, y = normalized_variance, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Variance/Mean') +
  ylim(c(min(H1_FL_treatment$normalized_variance)-10, max(H1_FL_treatment$normalized_variance)+900)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$normalized_variance)+500, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Coefficient of dispersion')

p8<-ggplot(H1_FL_treatment, aes(x = treatment, y = Mean, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Fluorescence (a.u.)') +
  ylim(c(min(H1_FL_treatment$Mean)-100, max(H1_FL_treatment$Mean)+1400)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_FL_treatment$Mean)+800, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Mean nuclear fluorescence')

df = H1_FL_treatment[,c("treatment","number_of_nuclear_speckles","integrated_volume","normalized_variance","Mean")]
df$treatment = ifelse(df$treatment == "DMSO", "control", "FL")
colnames(df) = c("treatment","number_of_foci","total_volume","normalized_variance","mean_nuclear_fluorescence")
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S2_m.txt", row.names = F, col.names = T, quote = F, sep = "\t")


# png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_E_G.png", width = 8, height = 3, res = 200, units = "in")
# plot_grid(p1,p2,p3, ncol=3)
# dev.off()


##### RNase
H1_ctr_SC35<-H1_SC35 %>%
  filter(treatment=='RNase control')

H1_RNase_treat_SC35<-H1_SC35 %>%
  filter(treatment=='RNase treatment')

H1_RNase_treatment<-rbind(H1_ctr_SC35, H1_RNase_treat_SC35)

my_comparisons= list(c("RNase control", "RNase treatment"))

p9<-ggplot(H1_RNase_treatment, aes(x = treatment, y = number_of_nuclear_speckles, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Counts') +
  ylim(c(min(H1_RNase_treatment$number_of_nuclear_speckles)-20, max(H1_RNase_treatment$number_of_nuclear_speckles)+25)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$number_of_nuclear_speckles)+15, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Number of foci')

p10<-ggplot(H1_RNase_treatment, aes(x=treatment, y = integrated_volume, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab(expression(paste('Total volume (', mu, m^3, ')'))) +
  ylim(c(min(H1_RNase_treatment$integrated_volume)-1, max(H1_RNase_treatment$integrated_volume)+4)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$integrated_volume)+2, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Integrated volume')

p11<-ggplot(H1_RNase_treatment, aes(x = treatment, y = normalized_variance, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Variance/Mean') +
  ylim(c(min(H1_RNase_treatment$normalized_variance)-10, max(H1_RNase_treatment$normalized_variance)+900)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$normalized_variance)+500, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Coefficient of dispersion')

p12<-ggplot(H1_RNase_treatment, aes(x = treatment, y = Mean, fill=treatment))+ 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  geom_jitter(width = 0.2, size = jitter_size) +
  xlab('')+
  ylab('Fluorescence (a.u.)') +
  ylim(c(min(H1_RNase_treatment$Mean)-100, max(H1_RNase_treatment$Mean)+1400)) +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y=max(H1_RNase_treatment$Mean)+800, size=statistic_test_size) +
  theme_bw() +
  theme(text = element_text(size=label_size),
        plot.title = element_text(size=label_size),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=label_size),
        axis.title.y = element_text(size=label_size),
        legend.position = "none") +
  ggtitle('Mean nuclear fluorescence')

df = H1_RNase_treatment[,c("treatment","number_of_nuclear_speckles","integrated_volume","normalized_variance","Mean")]
df$treatment = ifelse(df$treatment == "RNase control", "control", "RNase")
colnames(df) = c("treatment","number_of_foci","total_volume","normalized_variance","mean_nuclear_fluorescence")
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S2_n.txt", row.names = F, col.names = T, quote = F, sep = "\t")


png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_IF_SC35.png", width = 4.5, height = 6, res = 300, units = "in")
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol=3, byrow = F)
dev.off()






################### Figure IF-N

temp1 = cbind(seq(from = 50000, to = 20000000, by = 20000000/10),contact_frequency_control,1)
colnames(temp1) = c("distance","contact","Sample")

temp2 = cbind(seq(from = 50000, to = 20000000, by = 20000000/10),contact_frequency_NH4OAc,2)
colnames(temp2) = c("distance","contact","Sample")

temp3 = cbind(seq(from = 50000, to = 20000000, by = 20000000/10),contact_frequency_FL,3)
colnames(temp3) = c("distance","contact","Sample")

temp4 = cbind(seq(from = 50000, to = 20000000, by = 20000000/10),contact_frequency_RNase,4)
colnames(temp4) = c("distance","contact","Sample")

temp_mat = rbind(temp1,temp2,temp3,temp4)
temp_mat = as.data.frame(temp_mat)
temp_mat[which(temp_mat$Sample==1),"Sample"]='Control'
temp_mat[which(temp_mat$Sample==2),"Sample"]='NH4OAc'
temp_mat[which(temp_mat$Sample==3),"Sample"]='FL'
temp_mat[which(temp_mat$Sample==4),"Sample"]='RNase'
temp_mat$Sample = factor(temp_mat$Sample, levels = c("Control", "NH4OAc", "FL", "RNase"))

write.table(temp_mat,"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/HiC_distance_interaction_plot.txt", row.names = F, col.names = T, sep = "\t", quote = F)


library(scales)
#png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_F.png", height = 2.222, width = 2.222, units = "in", res = 300)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_F.pdf", height = 1.6, width = 1.8)
ggplot(temp_mat, aes(x=distance, y=contact, color=Sample)) + 
  geom_line(size=0.2) +
  scale_x_continuous(limits = c(50000,20000000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.5,1000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Genomic Distance (bp)", y="Hi-C Interaction Frequency") +
  scale_color_manual(values = c("black","green","red","magenta")) +
  theme(#legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.line.y = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.ticks = element_line(colour = "black", size = 0.4, linetype = "solid"),
    #axis.ticks.y = element_line(colour = "black", size = 0.3, linetype = "solid"),
    text = element_text(size=7),
    #axis.ticks.length = unit(0.1, "cm")
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  geom_vline(xintercept = 2000000, linetype = "dashed", size = 0.2)
#annotation_logticks()
dev.off()


#png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_3f1.png", height = 1.8, width = 2.2, units = "in", res = 300)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_3f1.pdf", height = 3, width = 2.5)
ggplot(temp_mat[which(temp_mat$distance %in% seq(50000, 2050000, by = 1000000)),], aes(x=distance, y=contact, color=Sample)) + 
  geom_line(size=0.3) +
  scale_x_continuous(limits = c(50000,2100000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(6.5,800), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Genomic Distance (bp)", y="Hi-C Interaction Frequency") +
  scale_color_manual(values = c("black","green","red","magenta")) +
  theme(#legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.line.y = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.ticks = element_line(colour = "black", size = 0.4, linetype = "solid"),
    #axis.ticks.y = element_line(colour = "black", size = 0.3, linetype = "solid"),
    text = element_text(size=7),
    #axis.ticks.length = unit(0.1, "cm")
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)) +
  guides(color = guide_legend(override.aes = list(size = 0.5)))
#annotation_logticks()
dev.off()

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_3f2.pdf", height = 3, width = 2.5)
ggplot(temp_mat[which(temp_mat$distance %in% seq(2050000, 18050000, by = 1000000)),], aes(x=distance, y=contact, color=Sample)) + 
  geom_line(size=0.3) +
  scale_x_continuous(limits = c(2000000,18100000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.6,7.5), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Genomic Distance (bp)", y="Hi-C Interaction Frequency") +
  scale_color_manual(values = c("black","green","red","magenta")) +
  theme(#legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.line.y = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.ticks = element_line(colour = "black", size = 0.4, linetype = "solid"),
    #axis.ticks.y = element_line(colour = "black", size = 0.3, linetype = "solid"),
    text = element_text(size=7),
    #axis.ticks.length = unit(0.1, "cm")
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)) +
  guides(color = guide_legend(override.aes = list(size = 0.5)))
#annotation_logticks()
dev.off()


# pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_3f2.pdf", height = 3, width = 1.7)
# ggplot(temp_mat[which(temp_mat$distance %in% seq(2050000, 10050000, by = 1000000)),], aes(x=distance, y=contact, color=Sample)) + 
#   geom_line(size=0.3) +
#   scale_x_continuous(limits = c(2000000,10100000), trans=log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_continuous(limits = c(0.1,4.1), trans=log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#   labs(x="Genomic Distance (bp)", y="Hi-C Interaction Frequency") +
#   scale_color_manual(values = c("black","green","red","magenta")) +
#   theme(#legend.title = element_blank(),
#     legend.position = "none",
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     # axis.line.x = element_line(colour = "black", size = 0.4, linetype = "solid"),
#     # axis.line.y = element_line(colour = "black", size = 0.4, linetype = "solid"),
#     # axis.ticks = element_line(colour = "black", size = 0.4, linetype = "solid"),
#     #axis.ticks.y = element_line(colour = "black", size = 0.3, linetype = "solid"),
#     text = element_text(size=7),
#     #axis.ticks.length = unit(0.1, "cm")
#     axis.text.x = element_text(size = 7),
#     axis.text.y = element_text(size = 7)) +
#   guides(color = guide_legend(override.aes = list(size = 0.5)))
# #annotation_logticks()
# dev.off()
# 
# 
# pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_3f3.pdf", height = 3, width = 1.7)
# ggplot(temp_mat[which(temp_mat$distance %in% seq(12050000, 18050000, by = 2000000)),], aes(x=distance, y=contact, color=Sample)) + 
#   geom_line(size=0.3) +
#   scale_x_continuous(limits = c(12000000,18100000), trans=log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_continuous(limits = c(0.65,1.1), trans=log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#   labs(x="Genomic Distance (bp)", y="Hi-C Interaction Frequency") +
#   scale_color_manual(values = c("black","green","red","magenta")) +
#   theme(#legend.title = element_blank(),
#     legend.position = "none",
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     # axis.line.x = element_line(colour = "black", size = 0.4, linetype = "solid"),
#     # axis.line.y = element_line(colour = "black", size = 0.4, linetype = "solid"),
#     # axis.ticks = element_line(colour = "black", size = 0.4, linetype = "solid"),
#     #axis.ticks.y = element_line(colour = "black", size = 0.3, linetype = "solid"),
#     text = element_text(size=7),
#     #axis.ticks.length = unit(0.1, "cm")
#     axis.text.x = element_text(size = 7),
#     axis.text.y = element_text(size = 7)) +
#   guides(color = guide_legend(override.aes = list(size = 0.5)))
# #annotation_logticks()
# dev.off()


######
temp1 = cbind(seq(from = 50000, to = 20000000, by = 20000000/20),contact_frequency_control_20,1)
colnames(temp1) = c("distance","contact","Sample")
x = ecdf(temp1[,"contact"])
temp1 = cbind(temp1, x(temp1[,"contact"]))
colnames(temp1) = c("distance","contact","Sample","ecdf")

temp2 = cbind(seq(from = 50000, to = 20000000, by = 20000000/20),contact_frequency_NH4OAc_20,2)
colnames(temp2) = c("distance","contact","Sample")
x = ecdf(temp2[,"contact"])
temp2 = cbind(temp2, x(temp2[,"contact"]))
colnames(temp2) = c("distance","contact","Sample","ecdf")

temp3 = cbind(seq(from = 50000, to = 20000000, by = 20000000/20),contact_frequency_FL_20,3)
colnames(temp3) = c("distance","contact","Sample")
x = ecdf(temp3[,"contact"])
temp3 = cbind(temp3, x(temp3[,"contact"]))
colnames(temp3) = c("distance","contact","Sample","ecdf")

temp4 = cbind(seq(from = 50000, to = 20000000, by = 20000000/20),contact_frequency_RNase_20,4)
colnames(temp4) = c("distance","contact","Sample")
x = ecdf(temp4[,"contact"])
temp4 = cbind(temp4, x(temp4[,"contact"]))
colnames(temp4) = c("distance","contact","Sample","ecdf")

temp_mat = rbind(temp1,temp2,temp3,temp4)
temp_mat = as.data.frame(temp_mat)
temp_mat[which(temp_mat$Sample==1),"Sample"]='Control'
temp_mat[which(temp_mat$Sample==2),"Sample"]='NH4OAc'
temp_mat[which(temp_mat$Sample==3),"Sample"]='FL'
temp_mat[which(temp_mat$Sample==4),"Sample"]='RNase'
temp_mat$Sample = factor(temp_mat$Sample, levels = c("Control", "NH4OAc", "FL", "RNase"))


pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_3g_CDF.pdf", height = 2, width = 2)
ggplot(temp_mat, aes(x=contact, y=ecdf, color=Sample)) + 
  geom_line(size=0.2) +
  labs(x="Hi-C Interaction Frequency", y="CDF") +
  scale_x_continuous(limits = c(0.5,1000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("black","green","red","magenta")) +
  theme(#legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.line.y = element_line(colour = "black", size = 0.4, linetype = "solid"),
    # axis.ticks = element_line(colour = "black", size = 0.4, linetype = "solid"),
    #axis.ticks.y = element_line(colour = "black", size = 0.3, linetype = "solid"),
    text = element_text(size=7),
    #axis.ticks.length = unit(0.1, "cm")
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)) 
#annotation_logticks()
dev.off()






###### Figure IF-F
df1 = melt(df_DRL_50000)
df1$variable = factor(df1$variable, levels = c("Control","NH4OAc", "FL", "RNase"))

my_comparisons = list(c("Control","NH4OAc"),  c("Control","FL"), c("Control","RNase"))

wilcox.test(df1[which(df1$variable == "Control"),2],df1[which(df1$variable == "NH4OAc"),2])
wilcox.test(df1[which(df1$variable == "Control"),2],df1[which(df1$variable == "RNase"),2])
wilcox.test(df1[which(df1$variable == "Control"),2],df1[which(df1$variable == "FL"),2])

png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_F.png", height = 2.222, width = 2.222, units = "in", res = 300)
ggplot(df1, aes(y=value, x=variable)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 3, label.y = c(6,8,10)) +
  scale_y_continuous(breaks = c(-8,-4,0,4,8)) +
  ylim(c(-8,11)) +
  labs(x = "", y = "Log distal-to-local ratio") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black"))
dev.off()

### Figrue IF-H
df2 = rbind(df_local,df_distal)
df2$variable = factor(df2$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))

png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_G.png", height = 2.222, width = 2.222, units = "in", res = 300)
ggplot(df2, aes(y=value, x=variable, fill=label)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  labs(x = "", y = "Number of interactions") + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black"))
dev.off()

########## Figure IF I-J-K
p1 <- make_scatter_plot_switch(full_eigen_data_500000,"NH4OAc","500000")
p2 <- make_scatter_plot_switch(full_eigen_data_500000,"FL","500000")
p3 <- make_scatter_plot_switch(full_eigen_data_500000,"RNase","500000")

png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_I_J_K.png", width = 1.6658, height = 5, units = "in", res = 300)
plot_grid(p1,p2,p3, ncol = 1)
dev.off()



########## Figure IF L-M-N
figure_if_LMN <- function(eigen_file,
                          input_bedpe_file,
                          resolution,
                          log_margi=F,
                          fig_title,
                          force_compartment_control=F,
                          force_imargi_control=F,
                          MALAT1_cov=F){
  
  
  compartment = read.table(eigen_file, header = T)
  colnames(compartment) = c("chr","coord","eigen")
  compartment = compartment[which(compartment$chr != "chrY"),]
  
  
  input_bedpe = data.frame(fread(input_bedpe_file))
  
  input_data = input_bedpe[which(input_bedpe[,1] %in% hg38_chromosomes & 
                                   input_bedpe[,4] %in% hg38_chromosomes),1:6]
  rm(input_bedpe)
  input_data[,2] = input_data[,2] + 1 # 1-based 
  input_data[,5] = input_data[,5] + 1 # 1-based 
  input_data[,1] = as.character(input_data[,1])
  input_data[,4] = as.character(input_data[,4])
  colnames(input_data) = paste0("V",seq(1,6))
  
  Gr_RNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,1])),
    ranges = IRanges(input_data[,2], end = input_data[,3], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  Gr_DNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,4])),
    ranges = IRanges(input_data[,5], end = input_data[,6], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  rm(input_data)
  
  cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
                                       genome_gr = Gr_hg38, 
                                       window_size = resolution,
                                       amplifier = c(1,1), 
                                       threshold = NA, 
                                       names = c("RNA_end","DNA_end"))
  
  
  output = c()
  output = c(output, cor(compartment$eigen, cov_data$RNA_end, method = "pearson"))
  output = c(output, cor(compartment$eigen, cov_data$DNA_end, method = "pearson"))
  output = c(output, cor(compartment$eigen, cov_data$RNA_end, method = "spearman"))
  output = c(output, cor(compartment$eigen, cov_data$DNA_end, method = "spearman"))
  names(output) = c("PCC_RNA","PCC_DNA","SCC_RNA","SCC_DNA")
  
  df_aov = data.frame(ifelse(compartment$eigen < 0, "B", "A"), cov_data$DNA_end)
  colnames(df_aov) = c("eigen","RAL")
  df_aov$eigen = factor(df_aov$eigen)
  res.aov <- aov(RAL ~ eigen, data = df_aov)
  
  output_list = list(output, res.aov)
  
  df = data.frame(compartment$eigen, cov_data$RNA_end, cov_data$DNA_end)
  write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/df_EIGEN_B.txt", row.names = F, col.names = T, sep = "\t", quote = F)
  df = read.table("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/df_EIGEN_B.txt", header = T)
  
  annotations_RNA <- data.frame(
    xpos = c(Inf),
    ypos =  c(-Inf),
    #annotateText = c(paste0("PCC: ",round(output["PCC_RNA"],2),"\nSCC: ",round(output["SCC_RNA"],2))),
    annotateText = c(paste0("SCC: ",round(output["SCC_RNA"],2))),
    hjustvar = c(1) ,
    vjustvar = c(-0.5))
  
  annotations_DNA <- data.frame(
    xpos = c(Inf),
    ypos =  c(-Inf),
    #annotateText = c(paste0("PCC: ",round(output["PCC_DNA"],2),"\nSCC: ",round(output["SCC_DNA"],2))),
    annotateText = c(paste0("SCC: ",round(output["SCC_DNA"],2))),
    hjustvar = c(1) ,
    vjustvar = c(-2))
  
  label_size = 8
  
  y_lab = ifelse(force_compartment_control == T, "Eigenvector control", "Eigenvector")
  
  p1<-ggplot(df, aes(x=cov_data.DNA_end, y=compartment.eigen)) + 
    geom_point(color = "purple", size = 0.2) +
    #xlab("DNA end read coverage") +
    xlab("RAL") +
    ylab("A/B") +
    #scale_x_continuous(labels = c("0","10,000","20,000"), breaks = c(0,10000,20000)) + # NH4OAc
    #scale_x_continuous(labels = c("0","20,000","40,000"), breaks = c(0,20000,40000)) + # FL
    #scale_x_continuous(labels = c("0","2,000","4,000"), breaks = c(0,2000,4000)) + # RNase
    
    # Figure S-RAL
    #scale_x_continuous(labels = c("0","35,000","70,000"), breaks = c(0,35000,70000)) + # control
    #scale_x_continuous(labels = c("0","10,000","20,000"), breaks = c(0,10000,20000)) + # NH4OAc  
    #scale_x_continuous(labels = c("0","20,000","40,000"), breaks = c(0,20000,40000)) + # FL  
    #scale_x_continuous(labels = c("0","2,000","4,000"), breaks = c(0,2000,4000)) + # RNase
    
    # Figure S_EIGEN_C_D
  # scale_x_continuous(labels = c("0","10,000","20,000"), breaks = c(0,10000,20000)) + # HFF
  scale_x_continuous(labels = c("0","100,000","200,000"), breaks = c(0,100000,200000)) + # K562
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size),
          axis.title.x = element_text(size=label_size),
          plot.title = element_text(size = label_size))
  #ggtitle("RAL vs A/B compartment")
  #geom_text(data=annotations_DNA,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=3)
  
  #png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_EIGEN_N.png", height = 1.6658, width = 1.6658, units = "in", res = 300)
  #png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_RAL_D.png", height = 1.6658, width = 1.6658, units = "in", res = 300)
  png(paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/",fig_title,".png"), height = 1.6658, width = 1.6658, units = "in", res = 300)
  print(p1)
  dev.off()
  
  return(output)
}


sample = "iMARGI_H1_NH4OAc"
sample_name=gsub("iMARGI_","",sample)
eigen_sample = strsplit(sample,"_")[[1]][3]

figure_if_LMN(eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
              input_bedpe_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz"),
              resolution=500000,
              log_margi = F,
              fig_title=sample_name,
              force_compartment_control=F,
              force_imargi_control=F)


################# Figure S-RAL

sample = "iMARGI_H1_RNase"
sample_name=gsub("iMARGI_","",sample)
eigen_sample = strsplit(sample,"_")[[1]][3]
figure_if_LMN(eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
              input_bedpe_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
              resolution=500000,
              log_margi = F,
              fig_title=sample_name,
              force_compartment_control=F,
              force_imargi_control=F)


################# Figure S-EIGEN C-D
# scc HFF = 0.7189994
# scc K562 = 0.6975246

figure_if_LMN(eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/HFF_500000.txt"),
              input_bedpe_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_HFF_control/iMARGI_HFF_control.mapq30.200k.final.bedpe.gz"),
              resolution=500000,
              log_margi = F,
              fig_title="Figure_S_EIGEN_C",
              force_compartment_control=F,
              force_imargi_control=F)

figure_if_LMN(eigen_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/compartments/K562_500000.txt"),
              input_bedpe_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_K562_control/iMARGI_K562_control.mapq30.200k.final.bedpe.gz"),
              resolution=500000,
              log_margi = F,
              fig_title="Figure_S_EIGEN_D",
              force_compartment_control=F,
              force_imargi_control=F)





######## Figure IF_V
summary_TAD_label_bio_rep = data.frame(do.call("rbind", summary_TAD_label_bio_rep_list))
summary_TAD_label_bio_rep$sample = gsub("iMARGI_","",rownames(summary_TAD_label_bio_rep))
summary_TAD_label_bio_rep$sample_type = paste0("H1_",sapply(summary_TAD_label_bio_rep$sample, function(x){strsplit(x, "_")[[1]][2]}))
summary_TAD_label_bio_rep[,1:3] = summary_TAD_label_bio_rep[,1:3] / rowSums(summary_TAD_label_bio_rep[,1:3])

summary_TAD_label_sample_mean = aggregate(summary_TAD_label_bio_rep[,1:3],
                                          by = list(sample_type = summary_TAD_label_bio_rep$sample_type),
                                          FUN="mean")
df1 = melt(summary_TAD_label_sample_mean)

summary_TAD_label_sample_sd = aggregate(summary_TAD_label_bio_rep[,1:3],
                                        by = list(sample_type = summary_TAD_label_bio_rep$sample_type),
                                        FUN="sd")[,2:4] / sqrt(c(4,2,2,2))
summary_TAD_label_sample_sd$sample_type = summary_TAD_label_sample_mean$sample_type
df2 = melt(summary_TAD_label_sample_sd)

df = merge(df1,df2,by=c("sample_type","variable"))
colnames(df)[3:4] = c("value", "error")
df$error_pos = NA
df$error_pos[df$variable == "neither"] =  df$value[df$variable == "neither"]
df$error_pos[df$variable == "enrichment"] = df$value[df$variable == "neither"] + df$value[df$variable == "enrichment"]
df$error_pos[df$variable == "depletion"] = df$value[df$variable == "neither"] + df$value[df$variable == "enrichment"] + df$value[df$variable == "depletion"]

png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_IF_V.png", width = 6, height = 5, res = 200, units = "in")
ggplot(df, aes(fill=variable, y=value, x=sample_type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('green','blue',"#808080")) +
  labs(x = "", y = "Percentage of TADs") +
  scale_y_continuous(name="Percentage of TADs", breaks = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75","1")) +
  geom_errorbar(aes(ymin=error_pos-error, ymax=error_pos+error), width=.2, size = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_blank())
dev.off()


# ###### Figure CONC-B-E
# plot_list = list()
# k = 1
# for (my_sample in all_imargi_samples[1:4]){
#   
#   my_cutoffs = c(0,150)
#   
#   p1 <- make_complex_heatmap_margi(
#     chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#     sample=my_sample,
#     contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/50000_filter1k/chr1_chr1.txt"),
#     sparse=T,
#     chr_row="chr1",
#     chr_col="chr1",
#     bin_size=50000,
#     select_coord=my_coord,
#     contacts_cutoffs=my_cutoffs,
#     my_colorbar=c("white","red"),
#     Gr_row_annotation=NA,
#     Gr_col_annotation=NA,
#     Gr_stripes_row=NA,
#     Gr_stripes_col=NA,
#     tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label.txt"),#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
#     loops_file_path = NA,#"/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
#     feature_color = "blue")#dbdeff")
#   
#   plot_list[[k]] = p1
#   k = k + 1
#   
# }
# 
# plot_list_grob = lapply(plot_list, as.grob)
# png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_CONC_B_E.png", width = 16, height = 4, units = "in", res = 200)
# plot_grid(plotlist = plot_list_grob, nrow=1)
# dev.off()


###### Figure LOOPD A
print_loop_numbers <-  function(loop_file){
  print(nrow(loop_file))
  print(nrow(loop_file[which(loop_file$V3- loop_file$V2 == 10000),]))
  print(nrow(loop_file[which(loop_file$V3- loop_file$V2 == 5000),]))
}



loops_control_10000 = loops_control[which(loops_control$V3- loops_control$V2 == 10000),]
loops_NH4OAc_10000 = loops_NH4OAc[which(loops_NH4OAc$V3- loops_NH4OAc$V2 == 10000),]
loops_FL_10000 = loops_FL[which(loops_FL$V3- loops_FL$V2 == 10000),]
loops_RNase_10000 = loops_RNase[which(loops_RNase$V3- loops_RNase$V2 == 10000),]


temp = rbind(loops_control,loops_NH4OAc,loops_FL,loops_RNase)
loops_union = temp[!duplicated(temp[,c(1,2,3,4,5,6)]),]
loops_union$V1 = as.character(loops_union$V1)
loops_union$V4 = as.character(loops_union$V4)

loops_union_5000 = loops_union[which(loops_union$V3- loops_union$V2 == 5000),]
loops_union_10000 = loops_union[which(loops_union$V3- loops_union$V2 == 10000),]

##########
loops_union_5000_ext = loops_union_5000
loops_union_5000$V2 = loops_union_5000$V2 - 5000
loops_union_5000$V3 = loops_union_5000$V3 + 5000
loops_union_5000$V5 = loops_union_5000$V5 - 5000
loops_union_5000$V6 = loops_union_5000$V6 + 5000

loops_union_10000_ext = loops_union_10000
loops_union_10000$V2 = loops_union_10000$V2 - 10000
loops_union_10000$V3 = loops_union_10000$V3 + 10000
loops_union_10000$V5 = loops_union_10000$V5 - 10000
loops_union_10000$V6 = loops_union_10000$V6 + 10000

loops_union_ext = rbind(loops_union_5000_ext,loops_union_10000_ext)

Gr_loop_union_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_ext$V1)),
  ranges = IRanges(loops_union_ext$V2+1, end = loops_union_ext$V3, names = c(1:nrow(loops_union_ext))),
  strand = Rle(strand('*')))

Gr_loop_union_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_ext$V1)),
  ranges = IRanges(loops_union_ext$V5+1, end = loops_union_ext$V6, names = c(1:nrow(loops_union_ext))),
  strand = Rle(strand('*')))


###
Gr_loop_union_10000_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_10000$V1)),
  ranges = IRanges(loops_union_10000$V2+1-10000, end = loops_union_10000$V3+10000, names = c(1:nrow(loops_union_10000))),
  strand = Rle(strand('*')))

Gr_loop_union_10000_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_10000$V1)),
  ranges = IRanges(loops_union_10000$V5+1-10000, end = loops_union_10000$V6+10000, names = c(1:nrow(loops_union_10000))),
  strand = Rle(strand('*')))




library(UpSetR)
png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_LOOPD_A_10000_extended.png", width = 5, height = 5.2, units = "in", res = 400)
upset(df_upset_H1_loops_10000, 
      sets = c("RNase","FL","NH4OAc","Control"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      #text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
      text.scale = c(1.8, 1.8, 0, 0, 1.6, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()


###### Figure LOOPD D
sample = c(rep("Control",length(distance_anchor_control)),
           rep("FL",length(distance_anchor_FL)),
           rep("NH4OAc",length(distance_anchor_NH4OAc)),
           rep("RNase",length(distance_anchor_RNase)))
distances = c(distance_anchor_control,
              distance_anchor_FL, 
              distance_anchor_NH4OAc, 
              distance_anchor_RNase)
df = data.frame(sample, distances)

df$sample = factor(df$sample, levels = c("Control","NH4OAc", "FL","RNase"))

my_comparisons = list(c("Control","FL"), c("Control","RNase"))

wilcox.test(df[which(df$sample == "Control"),2],df[which(df$sample == "FL"),2])
wilcox.test(df[which(df$sample == "Control"),2],df[which(df$sample == "RNase"),2])


png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_LOOPD_D.png", width = 2.222, height = 2.222, units = "in", res = 300)
ggplot(df, aes(x=sample, y=log10(distances))) + 
  geom_boxplot(outlier.size = 0.2, lwd = 0.3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 3, label.y = c(7.3,7.7)) +
  ylim(c(4.5,8)) +
  labs(x="", y="Log10 anchor points distance") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))
dev.off()



##### Figure LOOPD F

png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_LOOPD_F.png", width = 3.5, height = 2, units = "in", res = 300)
ggtern(data=df_p2m_union_10000_rel, aes(x=FL,y=NH4OAc, z=RNase)) +
  geom_point(size=0.1) +
  theme_rgbw() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=8))
dev.off()



##### Figure S3

plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  my_cutoffs = c(0,30)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr12_chr12.txt"),
    sparse=T,
    chr_row="chr12",
    chr_col="chr12",
    bin_size=10000,
    select_coord=rep(c(52700000,53700000),2),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,#Gr_annotations_edb_hg38,
    Gr_col_annotation=NA,#Gr_annotations_edb_hg38,
    annotation_type = "genes",
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
}

plot_list_grob = lapply(plot_list, as.grob)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_LOOPD_F_I.pdf", width = 4, height = 16)
plot_grid(plotlist = plot_list_grob, ncol = 1)
dev.off()

make_annotation_tracks(chr="chr20", start=52700000, end=53700000, sample_name="control", panel = "B")


##### Figure LOOPD G-J
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  p2 <- make_complex_heatmap_margi(
    chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
    sparse=T,
    chr_row="chr20",
    chr_col="chr20",
    bin_size=10000,
    select_coord=rep(c(47000000,47700000),2),
    contacts_cutoffs=c(0,40),
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA, #paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label_2.txt"),
    loops_file_path = loop_file_path_temp,
    feature_color = "blue",
    color_tads = F)
  
  plot_list[[k]] = p2
  k = k + 1
  
  
}

plot_list_grob = lapply(plot_list, as.grob)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_LOOPD_G_J.pdf", width = 16, height = 4)
plot_grid(plotlist = plot_list_grob, ncol = 4)
dev.off()

make_annotation_tracks(chr="chr20", start=47000000, end=47700000, sample_name="control", panel = "B")



###### Figure LOOPD N
plot_list = list()
k = 1
my_sample = all_imargi_samples[1]

temp = strsplit(my_sample,"_")[[1]][3]
loop_file_path_temp = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")

p2 <- make_complex_heatmap_margi(
  chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
  sparse=T,
  chr_row="chr20",
  chr_col="chr20",
  bin_size=10000,
  select_coord=c(47000000,48000000,47000000,48000000),
  contacts_cutoffs=c(0,40),
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label_2.txt"),
  loops_file_path = loop_file_path_temp,
  EP_loops_file = NA,#paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
  feature_color = "blue")

plot_list[[k]] = p2

plot_list_grob = lapply(plot_list, as.grob)
pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/heatmap_LOOPD_N.pdf", width = 6, height = 6)
plot_grid(plotlist = plot_list_grob, ncol = 1)
dev.off()


make_annotation_tracks(chr="chr20", start=47000000, end=48000000, sample_name="control")

make_margi_interaction_track(chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                             chr_plot = "chr20",
                             start_plot=47000000,
                             end_plot=48000000,
                             sample_label = "control")


make_plac_interaction_track(chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                            plac_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.txt",
                            chr_plot = "chr20",
                            start_plot=47000000,
                            end_plot=48000000,
                            read_length = 100)


### Treated samples
make_margi_interaction_track(chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_FL/iMARGI_H1_FL.mapq30.1k.final.bedpe",
                             chr_plot = "chr20",
                             start_plot=47000000,
                             end_plot=48000000,
                             sample_label = "FL")

make_margi_interaction_track(chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_NH4OAc/iMARGI_H1_NH4OAc.mapq30.1k.final.bedpe",
                             chr_plot = "chr20",
                             start_plot=47000000,
                             end_plot=48000000,
                             sample_label = "NH4OAc")

make_margi_interaction_track(chrSize_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_RNase/iMARGI_H1_RNase.mapq30.1k.final.bedpe",
                             chr_plot = "chr20",
                             start_plot=47000000,
                             end_plot=48000000,
                             sample_label = "RNase")


########## Figure S-TAD
sample = c("Control","NH4OAc","FL","RNase")
n_tads = c(length(Gr_tads_control_10000), length(Gr_tads_NH4OAc_10000), length(Gr_tads_FL_10000), length(Gr_tads_RNase_10000))
df = data.frame(sample, n_tads)
df$sample = factor(df$sample, levels = c("Control","NH4OAc","FL","RNase"))


p1<-ggplot(df, aes(y=n_tads, x=sample)) +
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "", y = "Number of TADs") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"))


summary_TAD_label_bio_rep = data.frame(do.call("rbind", summary_TAD_label_bio_rep_list))
summary_TAD_label_bio_rep$sample = gsub("iMARGI_","",rownames(summary_TAD_label_bio_rep))
summary_TAD_label_bio_rep$sample_type = paste0("H1_",sapply(summary_TAD_label_bio_rep$sample, function(x){strsplit(x, "_")[[1]][2]}))
summary_TAD_label_bio_rep[,1:3] = summary_TAD_label_bio_rep[,1:3] / rowSums(summary_TAD_label_bio_rep[,1:3])

summary_TAD_label_sample_mean = aggregate(summary_TAD_label_bio_rep[,1:3],
                                          by = list(sample_type = summary_TAD_label_bio_rep$sample_type),
                                          FUN="mean")
df1 = melt(summary_TAD_label_sample_mean)

summary_TAD_label_sample_sd = aggregate(summary_TAD_label_bio_rep[,1:3],
                                        by = list(sample_type = summary_TAD_label_bio_rep$sample_type),
                                        FUN="sd")[,2:4] / sqrt(c(4,2,2,2)) 
summary_TAD_label_sample_sd$sample_type = summary_TAD_label_sample_mean$sample_type
df2 = melt(summary_TAD_label_sample_sd)

df = merge(df1,df2,by=c("sample_type","variable"))
colnames(df)[3:4] = c("value", "error")
df$variable = factor(df$variable, levels = c("neither","depletion","enrichment"))

df$error_pos = NA
df$error_pos[df$variable == "enrichment"] =  df$value[df$variable == "enrichment"]
df$error_pos[df$variable == "depletion"] = df$value[df$variable == "depletion"] + df$value[df$variable == "enrichment"]
df$error_pos[df$variable == "neither"] = df$value[df$variable == "neither"] + df$value[df$variable == "enrichment"] + df$value[df$variable == "depletion"]
df$sample_type = gsub("H1_","",df$sample_type)
df$sample_type = gsub("control","Control",df$sample_type)
df$sample_type = factor(df$sample_type, levels = c("Control","NH4OAc","FL","RNase"))

p2<-ggplot(df, aes(fill=variable, y=value, x=sample_type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#808080",'green','blue')) +
  labs(x = "", y = "Percentage of TADs") +
  scale_y_continuous(name="Percentage of TADs", breaks = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75","1")) +
  geom_errorbar(aes(ymin=error_pos-error, ymax=error_pos+error), width=.2, size = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1),
        axis.text.y = element_text(size=10),
        legend.position = "none")



pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_TAD.pdf", width = 6, height = 3)
plot_grid(p2,p1)
dev.off()


######## Figure S_LOOPD
df1 = summary_apa_5000[which(summary_apa_5000$Index == "P2LL"),]
df1$Index = "All"

df2 = summary_apa_5000_union[which(summary_apa_5000_union$Index == "P2LL"),]
df2$Index = "Union"

df = rbind(df1,df2)
df$variable = factor(df$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))


p1<-ggplot(df, aes(fill=Index, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "P2LL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "none")


df1 = summary_apa_5000[which(summary_apa_5000$Index == "ZscoreLL"),]
df1$Index = "All"

df2 = summary_apa_5000_union[which(summary_apa_5000_union$Index == "ZscoreLL"),]
df2$Index = "Union"

df = rbind(df1,df2)
df$variable = factor(df$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))


p2<-ggplot(df, aes(fill=Index, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "ZscoreLL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "none")

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_LOOPD_AB.pdf", width = 4.5, height = 2.5)
plot_grid(p1,p2)
dev.off()

df = df_p2m_union_5000[,c(1,3,2,4)]
my_comparisons = list(c("Control","FL"), c("Control","RNase"))
df = melt(df)
df = df[which(df$value <= 20),]

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_S_LOOPD_C.pdf", width = 2.5, height = 2.5)
ggplot(df, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
  labs(x="", y="P2M") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 3, label.y = c(18,21)) +
  ylim(0,22) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10)) 
dev.off()


######### Figure S-RNAdomain B



######### Figure S4d
temp = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/C_stripe_rnaTreatment_upset.csv")
colnames(temp)[2:5] = c("Control","FL","NH4OAc","RNase")

library(UpSetR)
png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_RNAdomain_C.png", width = 5, height = 5, units = "in", res = 400)
upset(temp, 
      sets = c("RNase","FL","NH4OAc","Control"),
      sets.x.label = "Number of caRNA Domains",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      #text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
      text.scale = c(1.8, 1.8, 0, 0, 1.6, 0),
      mainbar.y.label = "Shared caRNA domains",
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()

######### Figure S-RNAdomain D
temp = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/D_stripe_cellLine_heightWidth_boxPlot.csv", header = F)


######### Figure S4a-c
temp1 = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/E_stripe_H1_heightGeneLength_scatter.csv")
temp2 = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/F_stripe_HFF_heightGeneLength_scatter.csv")
temp3 = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/G_stripe_K562_heightGeneLength_scatter.csv")

dot_size = 0.5

p1<-ggplot(temp1, aes(x=stripeHeight/10^6, y=geneLength/10^6)) + 
  geom_point(size = dot_size) +
  geom_smooth(method=lm, se = F,
              color="red", size = dot_size) +
  labs(x="Rectangular block's height (Mb)", y=paste0("Length of the longest", "\n", "overlapping gene (Mb)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none",
        title = element_text(size = 8)) +
  ggtitle("H1")

p2<-ggplot(temp2, aes(x=stripeHeight/10^6, y=geneLength/10^6)) + 
  geom_point(size = dot_size) +
  geom_smooth(method=lm, se = F,
              color="red", size = dot_size) +
  labs(x="Rectangular block's height (Mb)", y=paste0("Length of the longest", "\n", "overlapping gene (Mb)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none",
        title = element_text(size = 8)) +
  ggtitle("HFF")

p3<-ggplot(temp3, aes(x=stripeHeight/10^6, y=geneLength/10^6)) + 
  geom_point(size = dot_size) +
  geom_smooth(method=lm, se = F,
              color="red", size = dot_size) +
  labs(x="Rectangular block's height (Mb)", y=paste0("Length of the longest", "\n", "overlapping gene (Mb)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none",
        title = element_text(size = 8)) +
  ggtitle("K562")


png("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_RNAdomain_E_G.png", height = 2.2, width = 7, res = 200, units = "in")
plot_grid(p1,p2,p3, nrow=1)
dev.off()


######### Figure S4e
temp = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/H_stripe_rnaTreatment_height_boxPlot.csv", header = F)
temp = t(temp)[-1,]
temp = data.frame(temp)
colnames(temp) = c("Control", "NH4OAc","FL","RNase")

temp$Control = as.numeric(as.character(temp$Control))
temp$NH4OAc = as.numeric(as.character(temp$NH4OAc))
temp$FL = as.numeric(as.character(temp$FL))
temp$RNase = as.numeric(as.character(temp$RNase))

df = melt(temp)
df$variable = factor(df$variable, levels = c("Control","NH4OAc","FL","RNase"))
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S4e.txt", row.names = F, col.names = T, quote = F, sep = "\t")

my_comparisons = list(c("Control","NH4OAc"),  c("Control","FL"), c("Control","RNase"))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_RNAdomain_H.pdf", height = 2.3, width = 2.3)
ggplot(df, aes(x=variable, y=value)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3) +
  ylim(c(0,700)) +
  labs(x="", y="Rectangular block's height (kb)") +
  stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.x = 1.5, label.y = c(500, 600, 700), size=0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.position = "none")
dev.off()


######### Figure S4f
temp = read.csv("/mnt/extraids/SDSC_NFS/qizhijie/iMARGI_project/stripeCalling/code/results/Figure_S_RNAdomain_rawData/I_stripe_rnaTreatment_width_boxPlot.csv", header = F)
temp = t(temp)[-1,]
temp = data.frame(temp)
colnames(temp) = c("Control", "NH4OAc","FL","RNase")

temp$Control = as.numeric(as.character(temp$Control))
temp$NH4OAc = as.numeric(as.character(temp$NH4OAc))
temp$FL = as.numeric(as.character(temp$FL))
temp$RNase = as.numeric(as.character(temp$RNase))

df = melt(temp)
df$variable = factor(df$variable, levels = c("Control","NH4OAc","FL","RNase"))
write.table(df, "/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plot_source_data/Figure_S4f.txt", row.names = F, col.names = T, quote = F, sep = "\t")

my_comparisons = list(c("Control","NH4OAc"),  c("Control","FL"), c("Control","RNase"))

pdf("/mnt/extraids/SDSC_NFS/rcalandrelli/phase_separation/paper_plots/Figure_RNAdomain_I.pdf", height = 2.3, width = 2.3)
ggplot(df, aes(x=variable, y=value)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, lwd = 0.3) +
  labs(x="", y="Rectangular block's width (kb)") +
  ylim(c(0,3000)) +
  #stat_compare_means(comparisons = my_comparisons, method='wilcox.test', label='p.signif', label.y = c(15000, 16000, 17000), size=0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=7),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.position = "none")
dev.off()







