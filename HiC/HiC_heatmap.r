library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
require(dplyr)
library(igraph)

options(scipen=999)

##################################################################################################
make_complex_heatmap_HiC <- function(
  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample="test", # sample name
  contact_matrix_file, # contact matrix 
  sparse=T, # T if the contact matrix is sparse, F if dense
  oe=F, # T if the contact matrix is oe, F if not
  chr_row, # chr on the rows (RNA end)
  chr_col, # chr on the columns (DNA end)
  bin_size, # bin size or resolution
  select_coord=NA, # to select a portion of the matrix input c(start_row_coord, end_row_coord, start_col_coord, end_col_coord)
  contacts_cutoffs=c(0,NA), # minimum and maximum contact value cutoffs. Leave NA to do not have any max cutoff
  log_counts = F,
  my_colorbar=c("white","red"), # colorbar colors
  Gr_row_annotation=NA, # Granges object for annotation track on the rows
  Gr_col_annotation=NA, # Granges object for annotation track on the columns
  Gr_stripes_row=NA, # set True to add horizontal stripes over the heatmap corresponding to the row annotations
  Gr_stripes_col=NA, # set True to add vertical stripes over the heatmap corresponding to the column annotations
  compartments_file_path=NA,
  tads_file_path=NA,
  loops_file_path=NA,
  EP_loops_file=NA,
  feature_color = "#dbdeff" # set color of the features on the heatmap
){
  
  edge_to_matrix <- function(df,attrname="V3"){
    # df is source target wrights 
    g=graph.data.frame(df,directed = FALSE)
    matrix <- as_adjacency_matrix(g,type = "both", names = TRUE, sparse = FALSE, attr = attrname)
    return(matrix)
  }
  
  options(scipen=999)
  ########## Load genome data
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ##########  Tile genome and make genomic windows for selected chromosomes
  genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = bin_size, cut.last.tile.in.chrom = T)
  chr_row_genome_window = genome_window[seqnames(genome_window)==chr_row]
  chr_col_genome_window = genome_window[seqnames(genome_window)==chr_col]
  
  ########## Load contact matrix
  if (sparse == T){
    if (oe == F){
      temp = data.frame(fread(contact_matrix_file))
      contact_matrix = matrix(0, nrow = length(chr_row_genome_window), ncol = length(chr_col_genome_window))
      if (chr_row == chr_col){
        contact_matrix[as.matrix(temp[,c("V1","V2")] / bin_size + 1)] = temp$V3
        contact_matrix[as.matrix(temp[,c("V2","V1")] / bin_size + 1)] = temp$V3
      } else {
        contact_matrix[as.matrix(temp[,c("V1","V2")] / bin_size + 1)] = temp$V3
      }
      rownames(contact_matrix) = seq(1,nrow(contact_matrix)) # ???
      contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
      
    } else {
      contact_matrix <- data.frame(fread(contact_matrix_file)) %>% edge_to_matrix()
      contact_matrix[is.na(contact_matrix)] = 0
      contact_matrix <- cor(contact_matrix, method = "pearson") # pearson correlation matrix of the OE
      contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = -2
    }
    
  } else {
    contact_matrix = data.frame(fread(contact_matrix_file))
    contact_matrix = contact_matrix[,-1]
    rownames(contact_matrix) = seq(1,nrow(contact_matrix)) # ???
    contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
  }
  
  
  ########## Load compartment file if declared
  if (!is.na(compartments_file_path)){
    compartment = read.table(compartments_file_path, header = T)
    compartment = compartment[which(compartment$chr == chr_col),]
    compartment = cbind(compartment,compartment[,3])
    compartment[,4][which(compartment[,4]<0)] = "#0cad01"
    compartment[,4][which(compartment[,4]>=0)] = "#E41A1C"
    colnames(compartment)[4] = "color"
    
    compartment_ha = HeatmapAnnotation(eigen = anno_barplot(compartment$eigen,
                                                            bar_width = 1,
                                                            height = unit(1, "cm"),
                                                            gp = gpar(fill = compartment$color)),
                                       which = "column")
  }
  
  
  ########## Load TAD file if declared
  if (!is.na(tads_file_path)){
    
    ##### Using TAD file from HiC analysis without iMARGI enrichment and depletion information
    tads_file = read.table(paste0(tads_file_path,"/TADs_",chr_row,"/10000_blocks.bedpe"))[,c(1:3)]
    Gr_tads = GRanges(
      seqnames = Rle(paste0("chr",tads_file$V1)),
      ranges = IRanges(as.numeric(tads_file$V2), end = as.numeric(tads_file$V3), names = c(1:nrow(tads_file))),
      strand = Rle(strand('*')))
    Gr_tads = Gr_tads[seqnames(Gr_tads) == chr_row]
    
    overlaps = findOverlaps(Gr_tads, chr_row_genome_window)
    
    start_tad_bin = c()
    end_tad_bin = c()
    for (i in 1:length(Gr_tads)){
      temp_tad = subjectHits(overlaps[which(queryHits(overlaps)==i),])
      start_tad_bin = c(start_tad_bin,temp_tad[1])
      end_tad_bin = c(end_tad_bin,temp_tad[length(temp_tad)])
    }
    
    df_tad_bin = data.frame(start_tad_bin, end_tad_bin)
    for (i in 1:nrow(df_tad_bin)){
      start_tad = df_tad_bin[i,1]
      end_tad = df_tad_bin[i,2]
      contact_matrix[start_tad:end_tad,start_tad] = -4
      contact_matrix[start_tad:end_tad,end_tad] = -4
      contact_matrix[start_tad,start_tad:end_tad] = -4
      contact_matrix[end_tad,start_tad:end_tad] = -4
      }
    
    
    ##### Using TAD file from iMARGI analysis with iMARGI enrichment and depletion information
    # tads_file = read.table(tads_file_path, header = T)
    # Gr_tads = GRanges(
    #   seqnames = Rle(tads_file$seqnames),
    #   ranges = IRanges(as.numeric(tads_file$start), end = as.numeric(tads_file$end), names = c(1:nrow(tads_file))),
    #   strand = Rle(strand('*')),
    #   label = as.character(tads_file$label))
    # Gr_tads = Gr_tads[seqnames(Gr_tads) == chr_row]
    # 
    # overlaps = findOverlaps(Gr_tads, chr_row_genome_window)
    # 
    # start_tad_bin = c()
    # end_tad_bin = c()
    # for (i in 1:length(Gr_tads)){
    #   temp_tad = subjectHits(overlaps[which(queryHits(overlaps)==i),])
    #   start_tad_bin = c(start_tad_bin,temp_tad[1])
    #   end_tad_bin = c(end_tad_bin,temp_tad[length(temp_tad)])
    # }
    # 
    # label_bin = mcols(Gr_tads)[,"label"]
    # df_tad_bin = data.frame(start_tad_bin, end_tad_bin, label_bin)
    # 
    # for (i in 1:nrow(df_tad_bin)){
    #   start_tad = df_tad_bin[i,1]
    #   end_tad = df_tad_bin[i,2]
    #   if (df_tad_bin[i,3] == "enrichment"){
    #     contact_matrix[start_tad:end_tad,start_tad] = -2
    #     contact_matrix[start_tad:end_tad,end_tad] = -2
    #     contact_matrix[start_tad,start_tad:end_tad] = -2
    #     contact_matrix[end_tad,start_tad:end_tad] = -2
    #   }
    #   else if (df_tad_bin[i,3] == "depletion"){
    #     contact_matrix[start_tad:end_tad,start_tad] = -3
    #     contact_matrix[start_tad:end_tad,end_tad] = -3
    #     contact_matrix[start_tad,start_tad:end_tad] = -3
    #     contact_matrix[end_tad,start_tad:end_tad] = -3
    #   } else {
    #     contact_matrix[start_tad:end_tad,start_tad] = -4
    #     contact_matrix[start_tad:end_tad,end_tad] = -4
    #     contact_matrix[start_tad,start_tad:end_tad] = -4
    #     contact_matrix[end_tad,start_tad:end_tad] = -4
    #   }
    # }
    
  }
  
  ##########  Load loop file if declared
  if (!is.na(loops_file_path)){
    loops_file = read.table(paste0(loops_file_path,"/merged_loops.bedpe"))[,c(1:6)]
    
    Gr_loops_RNA = GRanges(
      seqnames = Rle(paste0("chr",loops_file$V1)),
      ranges = IRanges(as.numeric(loops_file$V2), end = as.numeric(loops_file$V3), names = c(1:nrow(loops_file))),
      strand = Rle(strand('*')))
    Gr_loops_RNA = Gr_loops_RNA[seqnames(Gr_loops_RNA) == chr_row]
    
    Gr_loops_DNA = GRanges(
      seqnames = Rle(paste0("chr",loops_file$V4)),
      ranges = IRanges(as.numeric(loops_file$V5), end = as.numeric(loops_file$V6), names = c(1:nrow(loops_file))),
      strand = Rle(strand('*')))
    Gr_loops_DNA = Gr_loops_DNA[seqnames(Gr_loops_DNA) == chr_row]
    
    overlaps_RNA = findOverlaps(Gr_loops_RNA, chr_row_genome_window)
    overlaps_DNA = findOverlaps(Gr_loops_DNA, chr_col_genome_window)
    
    start_loop_bin_row = c()
    end_loop_bin_row = c()
    for (i in 1:length(Gr_loops_RNA)){
      temp_loop = subjectHits(overlaps_RNA[which(queryHits(overlaps_RNA)==i),])
      start_loop_bin_row = c(start_loop_bin_row, temp_loop[1])
      end_loop_bin_row = c(end_loop_bin_row, temp_loop[length(temp_loop)])
    }
    
    start_loop_bin_col = c()
    end_loop_bin_col = c()
    for (i in 1:length(Gr_loops_DNA)){
      temp_loop = subjectHits(overlaps_DNA[which(queryHits(overlaps_DNA)==i),])
      start_loop_bin_col = c(start_loop_bin_col, temp_loop[1])
      end_loop_bin_col = c(end_loop_bin_col, temp_loop[length(temp_loop)])
    }
    
    df_loop_bin = data.frame(start_loop_bin_row, end_loop_bin_row, start_loop_bin_col, end_loop_bin_col)
    
    for (i in 1:nrow(df_loop_bin)){
      start_loop_row = df_loop_bin[i,1]
      end_loop_row = ifelse(df_loop_bin[i,2] > df_loop_bin[i,1], df_loop_bin[i,2] - 1, df_loop_bin[i,2])
      start_loop_col = df_loop_bin[i,3]
      end_loop_col = ifelse(df_loop_bin[i,4] > df_loop_bin[i,3], df_loop_bin[i,4] - 1, df_loop_bin[i,4])
      
      contact_matrix[start_loop_row:end_loop_row,start_loop_col] = -3
      contact_matrix[start_loop_row:end_loop_row,end_loop_col] = -3
      contact_matrix[start_loop_row,start_loop_col:end_loop_col] = -3
      contact_matrix[end_loop_row,start_loop_col:end_loop_col] = -3
      
      # Plot also loops under the diagonal
      contact_matrix[start_loop_col:end_loop_col,start_loop_row] = -3
      contact_matrix[start_loop_col:end_loop_col,end_loop_row] = -3
      contact_matrix[start_loop_col,start_loop_row:end_loop_row] = -3
      contact_matrix[end_loop_col,start_loop_row:end_loop_row] = -3
      
    }
  }
  
  ##########  Load enhancer-promoter loop file if declared
  if (!is.na(EP_loops_file)){
    loops_file = read.table(EP_loops_file)[,c(1:6,25)]
    loops_file = loops_file[which(!is.na(loops_file$V25)),]
    
    Gr_loops_RNA = GRanges(
      seqnames = Rle(paste0("chr",loops_file$V1)),
      ranges = IRanges(as.numeric(loops_file$V2), end = as.numeric(loops_file$V3), names = c(1:nrow(loops_file))),
      strand = Rle(strand('*')),
      type=loops_file$V25)
    Gr_loops_RNA = Gr_loops_RNA[seqnames(Gr_loops_RNA) == chr_row]
    
    Gr_loops_DNA = GRanges(
      seqnames = Rle(paste0("chr",loops_file$V4)),
      ranges = IRanges(as.numeric(loops_file$V5), end = as.numeric(loops_file$V6), names = c(1:nrow(loops_file))),
      strand = Rle(strand('*')),
      type=loops_file$V25)
    Gr_loops_DNA = Gr_loops_DNA[seqnames(Gr_loops_DNA) == chr_row]
    
    overlaps_RNA = findOverlaps(Gr_loops_RNA, chr_row_genome_window)
    overlaps_DNA = findOverlaps(Gr_loops_DNA, chr_col_genome_window)
    
    start_loop_bin_row = c()
    end_loop_bin_row = c()
    for (i in 1:length(Gr_loops_RNA)){
      temp_loop = subjectHits(overlaps_RNA[which(queryHits(overlaps_RNA)==i),])
      start_loop_bin_row = c(start_loop_bin_row, temp_loop[1])
      end_loop_bin_row = c(end_loop_bin_row, temp_loop[length(temp_loop)])
    }
    
    start_loop_bin_col = c()
    end_loop_bin_col = c()
    for (i in 1:length(Gr_loops_DNA)){
      temp_loop = subjectHits(overlaps_DNA[which(queryHits(overlaps_DNA)==i),])
      start_loop_bin_col = c(start_loop_bin_col, temp_loop[1])
      end_loop_bin_col = c(end_loop_bin_col, temp_loop[length(temp_loop)])
    }
    
    df_loop_bin = data.frame(start_loop_bin_row, end_loop_bin_row, start_loop_bin_col, end_loop_bin_col, mcols(Gr_loops_RNA)["type"])
    
    for (i in 1:nrow(df_loop_bin)){
      start_loop_row = df_loop_bin[i,1]
      end_loop_row = ifelse(df_loop_bin[i,2] > df_loop_bin[i,1], df_loop_bin[i,2] - 1, df_loop_bin[i,2])
      start_loop_col = df_loop_bin[i,3]
      end_loop_col = ifelse(df_loop_bin[i,4] > df_loop_bin[i,3], df_loop_bin[i,4] - 1, df_loop_bin[i,4])
      
      value = ifelse(df_loop_bin[i,5] == "enhancer_promoter", -2, 
                     ifelse(df_loop_bin[i,5] == "promoter_enhancer", -4, -5))
      
      contact_matrix[start_loop_row:end_loop_row,start_loop_col] = value
      contact_matrix[start_loop_row:end_loop_row,end_loop_col] = value
      contact_matrix[start_loop_row,start_loop_col:end_loop_col] = value
      contact_matrix[end_loop_row,start_loop_col:end_loop_col] = value
      
      # Plot also loops under the diagonal
      value = ifelse(df_loop_bin[i,5] == "enhancer_promoter", -4, 
                     ifelse(df_loop_bin[i,5] == "promoter_enhancer", -2, -5))
      
      contact_matrix[start_loop_col:end_loop_col,start_loop_row] = value
      contact_matrix[start_loop_col:end_loop_col,end_loop_row] = value
      contact_matrix[start_loop_col,start_loop_row:end_loop_row] = value
      contact_matrix[end_loop_col,start_loop_row:end_loop_row] = value
      
    }
  }
  
  
  
  ########## Select part of the contact matrix
  if (length(select_coord) == 4){
    coord_bin = ceiling(select_coord/bin_size)
    
    coord_bin[1] = coord_bin[1] + 1
    coord_bin[3] = coord_bin[3] + 1
    
    if (is.na(coord_bin[2])){
      coord_bin[2] = ceiling(chrSize[which(chrSize$V1 == chr_row),2] / bin_size)
    }
    
    if (is.na(coord_bin[4])){
      coord_bin[4] = ceiling(chrSize[which(chrSize$V1 == chr_col),2] / bin_size)
    }
    
    contact_matrix = apply(as.matrix(contact_matrix[coord_bin[1]:coord_bin[2],coord_bin[3]:coord_bin[4]]), 2, as.numeric)
    rownames(contact_matrix) = seq(1,nrow(contact_matrix)) # ???
  }
  
  if (log_counts == T){
    contact_matrix[which(contact_matrix > 0, arr.ind = T)] = log(contact_matrix[which(contact_matrix > 0, arr.ind = T)])
  }
  
  
  ##########  Check if making annotation tracks
  if (!is.na(Gr_row_annotation[1])){
    row_annotation_bin = unique(queryHits(findOverlaps(chr_row_genome_window,Gr_row_annotation)))
    if (length(select_coord) == 4){
      row_annotation_bin = row_annotation_bin[which(row_annotation_bin>=coord_bin[1] & row_annotation_bin<= coord_bin[2])] - coord_bin[1]
    }
    bar_annotation = rep(0,nrow(contact_matrix))
    bar_annotation[row_annotation_bin] = 1
    # same as rowAnnotation without the "which" parameter
    row_ha = HeatmapAnnotation(SE = anno_barplot(x=bar_annotation,
                                                 bar_width = 1,
                                                 width = unit(0.5, "cm"), 
                                                 axis = F),
                               which = "row")
  }
  
  if (!is.na(Gr_col_annotation[1])){
    col_annotation_bin = unique(queryHits(findOverlaps(chr_col_genome_window,Gr_col_annotation)))
    # Check and subset the annotation if only a part of the matrix was selected 
    if (length(select_coord) == 4){
      col_annotation_bin = col_annotation_bin[which(col_annotation_bin>=coord_bin[3] & col_annotation_bin<= coord_bin[4])] - coord_bin[3]
    }
    bar_annotation = rep(0,nrow(contact_matrix))
    bar_annotation[col_annotation_bin] = 1
    # same as columnAnnotation without the "which" parameter
    column_ha = HeatmapAnnotation(SE = anno_barplot(bar_annotation, 
                                                    bar_width = 1,
                                                    height = unit(0.5, "cm"), 
                                                    axis = F),
                                  which = "column")
  }
  
  ########## Check if making stripes
  zero_indexes = as.matrix(which(contact_matrix==0, arr.ind = T))
  if (!is.na(Gr_stripes_row[1])){
    row_stripe_bin = unique(queryHits(findOverlaps(chr_row_genome_window,Gr_stripes_row)))
    if (length(select_coord) == 4){
      row_stripe_bin = row_stripe_bin[which(row_stripe_bin>=coord_bin[1] & row_stripe_bin<= coord_bin[2])] - coord_bin[1]
    }
    contact_matrix[zero_indexes[which(zero_indexes[,1] %in% row_stripe_bin),]] = -1
  }
  if (!is.na(Gr_stripes_col[1])){
    col_stripe_bin = unique(queryHits(findOverlaps(chr_col_genome_window,Gr_stripes_col)))
    if (length(select_coord) == 4){
      col_stripe_bin = col_stripe_bin[which(col_stripe_bin>=coord_bin[1] & col_stripe_bin<= coord_bin[2])] - coord_bin[1]
    }
    contact_matrix[zero_indexes[which(zero_indexes[,2] %in% col_stripe_bin),]] = -1
  }
  
  ##########  Set colorbar
  if (oe == F){
    min_cutoff = as.numeric(contacts_cutoffs[1])
    if (is.na(contacts_cutoffs[2])){
      max_cutoff = max(contact_matrix)
    } else if (contacts_cutoffs[2] < 1){
      max_cutoff = round(as.numeric(quantile(contact_matrix[contact_matrix > 0], contacts_cutoffs[2])))
    }
    
    col_fun1 = colorRamp2(c(min_cutoff, max_cutoff), my_colorbar)
    col_fun = function(x) {
      ifelse(x == -1, "#dbdeff", 
             ifelse(x == -2, "green",  
                    ifelse(x == -3, "blue", # loops in general
                           ifelse(x == -4, "#b3b3b3", # neither TADs
                                  col_fun1(x)))))
    }
    attr(col_fun, "breaks") = c(min_cutoff, max_cutoff)
  } 
  else {
    col_fun1 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    col_fun = function(x) {
      ifelse(x == -2, "gray", col_fun1(x))
    }
    attr(col_fun, "breaks") = c(-1, 0, 1)
    
  }

  
  # Set axis labels
  # coord_labels_padding_left = 1
  # coord_labels_padding_right = 30
  
  # coord_labels_rows = rep("",nrow(contact_matrix))
  # coord_labels_columns = rep("",ncol(contact_matrix))
  # 
  # coord_labels_rows[1] = ifelse(is.na(select_coord), 0, select_coord[1])
  # coord_labels_rows[nrow(contact_matrix)] =  ifelse(is.na(select_coord), end(chr_row_genome_window[length(chr_row_genome_window)]), select_coord[2])
  # coord_labels_rows[round(nrow(contact_matrix)/2)] = paste0(chr_row,"\n(RNA)")
  # 
  # coord_labels_columns[1] = ifelse(is.na(select_coord), 0, select_coord[3])
  # coord_labels_columns[ncol(contact_matrix)] =  ifelse(is.na(select_coord), end(chr_col_genome_window[length(chr_col_genome_window)]), select_coord[4])
  # coord_labels_columns[round(ncol(contact_matrix)/2)] = paste0(chr_col,"\n(DNA)")
  
  # Make figure main title
  if(length(select_coord) != 1){
    if (is.na(select_coord[2])){
      select_coord[2] = chrSize[which(chrSize$V1 == chr_row),2]
    }
    
    if (is.na(select_coord[4])){
      select_coord[4] = chrSize[which(chrSize$V1 == chr_col),2]
    }
  }
  
  title_1 = paste0(sample, "\n")
  title_2 = ifelse(length(select_coord) == 1, 
                   paste0(chr_row,":0-",end(chr_row_genome_window[length(chr_row_genome_window)])), 
                   paste0(chr_row,":",select_coord[1],"-",select_coord[2]))
  if (chr_row != chr_col){
    title_3 = ifelse(length(select_coord) == 1, 
                     paste0("\n",chr_col,":0-",end(chr_col_genome_window[length(chr_col_genome_window)])),
                     paste0("\n",chr_col,":",select_coord[3],"-",select_coord[4]))
  } else {
    title_3 = ""
  }
  
  main_title = paste0(title_1, title_2, title_3)
  
  coord_fontsize = 12
  
  if (is.na(Gr_row_annotation[1]) & is.na(Gr_col_annotation[1])){
    p<-Heatmap(contact_matrix, 
               cluster_rows = FALSE, cluster_columns = FALSE, border = T,
               col = col_fun,
               show_row_names = F,
               show_column_names = F,
               #right_annotation = row_ha,
               #bottom_annotation = column_ha,
               na_col = feature_color,
               # row_labels = coord_labels_rows, row_names_gp = gpar(fontsize = coord_fontsize),
               # column_labels = coord_labels_columns, column_names_gp = gpar(fontsize = coord_fontsize),
               # row_names_side = "left",
               # column_names_side = "top",
               column_title = main_title,
               heatmap_legend_param = list(title = ""))
    # draw(p, heatmap_legend_list = list(Legend(at = c(-2, -1), 
    #                                            labels = c("Enrich", "Deplet", "Neither"), 
    #                                            legend_gp = gpar(fill = c("green", "blue", "yellow")))))
  }
  else if (!is.na(Gr_row_annotation[1]) & is.na(Gr_col_annotation[1])){
    p<-Heatmap(contact_matrix, 
               cluster_rows = FALSE, cluster_columns = FALSE, border = T,
               col = col_fun,
               show_row_names = F,
               show_column_names = F,
               right_annotation = row_ha,
               #bottom_annotation = column_ha,
               na_col = feature_color,
               row_names_side = "left",
               column_names_side = "top",
               # row_labels = coord_labels_rows, row_names_gp = gpar(fontsize = coord_fontsize),
               # column_labels = coord_labels_columns, column_names_gp = gpar(fontsize = coord_fontsize),
               column_title = main_title,
               heatmap_legend_param = list(title = ""))
  }
  else if (is.na(Gr_row_annotation[1]) & !is.na(Gr_col_annotation[1])){
    p<-Heatmap(contact_matrix, 
               cluster_rows = FALSE, cluster_columns = FALSE, border = T,
               col = col_fun,
               show_row_names = F,
               show_column_names = F,
               #right_annotation = row_ha,
               bottom_annotation = column_ha,
               na_col = feature_color,
               row_names_side = "left",
               column_names_side = "top",
               # row_labels = coord_labels_rows, row_names_gp = gpar(fontsize = coord_fontsize),
               # column_labels = coord_labels_columns, column_names_gp = gpar(fontsize = coord_fontsize),
               column_title = main_title,
               heatmap_legend_param = list(title = ""))
  }
  else if(!is.na(Gr_row_annotation[1]) & !is.na(Gr_col_annotation[1])){
    p<-Heatmap(contact_matrix, 
               cluster_rows = FALSE, cluster_columns = FALSE, border = T,
               col = col_fun,
               show_row_names = F,
               show_column_names = F,
               right_annotation = row_ha,
               bottom_annotation = column_ha,
               na_col = feature_color,
               row_names_side = "left",
               column_names_side = "top",
               # row_labels = coord_labels_rows, row_names_gp = gpar(fontsize = coord_fontsize),
               # column_labels = coord_labels_columns, column_names_gp = gpar(fontsize = coord_fontsize),
               column_title = main_title,
               heatmap_legend_param = list(title = ""))
  }
  
  if (!is.na(compartments_file_path)){
    p<-Heatmap(contact_matrix, 
               cluster_rows = FALSE, cluster_columns = FALSE, border = T,
               col = col_fun,
               show_row_names = F,
               show_column_names = F,
               #right_annotation = row_ha,
               bottom_annotation = compartment_ha,
               na_col = feature_color,
               row_names_side = "left",
               column_names_side = "top",
               # row_labels = coord_labels_rows, row_names_gp = gpar(fontsize = coord_fontsize),
               # column_labels = coord_labels_columns, column_names_gp = gpar(fontsize = coord_fontsize),
               column_title = main_title,
               heatmap_legend_param = list(title = ""))
  }
  
  
  return(p)
  
}

############ Test function
my_bin_size = 10000
my_sample = "H1_control"

chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes"
sample="test"
contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr20_chr20_",my_bin_size,".txt")
sparse=T
chr_row="chr20"
chr_col="chr20"
bin_size=my_bin_size
select_coord=c(47000000,48000000,47000000,48000000)
contacts_cutoffs=c(0,NA)
my_colorbar=c("white","red")
Gr_row_annotation=NA
Gr_col_annotation=NA
Gr_stripes_row=NA
Gr_stripes_col=NA
compartments_file_path= NA # "/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt"
#tads_file_path = "/dataOS/rcalandrelli/phase_separation/HiC/tads/H1_control_merged/"
tads_file_path = NA#"/dataOS/rcalandrelli/phase_separation/HiC/tads/H1_control_merged/"
loops_file_path = NA #"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged"
feature_color = "#dbdeff"



############ TADs heatmap plots
all_HiC_samples = paste0("H1_", c("control","FL","NH4OAc","RNase"))

plot_list = list()
k = 1
for (my_sample in all_HiC_samples){
  
  my_cutoffs = c(0,NA)
  my_bin_size = 50000
  
  p1 <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=my_bin_size,
    select_coord=c(32000000,44000000,32000000,44000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    compartments_file_path= NA,
    tads_file_path = paste0("/dataOS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
    loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
    feature_color = "blue")#dbdeff")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=my_bin_size,
    select_coord=c(72000000,84000000,72000000,84000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    compartments_file_path= NA,
    tads_file_path = paste0("/dataOS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
    loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
    feature_color = "blue")#dbdeff")
  
  plot_list[[k]] = p2
  k = k + 1
  
  p3 <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=my_bin_size,
    select_coord=rep(c(153000000,161000000,),2),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    compartments_file_path= NA,
    tads_file_path = paste0("/dataOS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
    loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
    feature_color = "blue")#dbdeff")
  
  plot_list[[k]] = p3
  k = k + 1
  
}

plot_list_grob = lapply(plot_list, as.grob)
png("/dataOS/rcalandrelli/phase_separation/HiC/tads/heatmaps/chr1_HiC_tad_heatmaps_FULL.png", width = 12, height = 16, units = "in", res = 200)
#png("/dataOS/rcalandrelli/phase_separation/HiC/tads/heatmaps/chr1_HiC_tad_heatmaps_98.png", width = 12, height = 16, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 3, labels = c("A",rep("",2), "B",rep("",2), "C",rep("",2), "D",rep("",2)))
dev.off()




############ Loops heatmap plots zoomed over patterns (only control)
plot_list = list()
k = 1
for (my_sample in all_HiC_samples[1]){
  
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_", my_sample, "_merged")
  
  my_cutoffs = c(0,0.99)
  my_bin_size = 10000
  
  p1 <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=10000,
    select_coord=c(3250000,3750000,3250000,3750000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/", my_sample, "_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=10000,
    select_coord=c(63200000,64700000,63200000,64700000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/", my_sample, "_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p2
  k = k + 1
}

#plot_list_grob = lapply(plot_list, as.grob)

pdf("/dataOS/rcalandrelli/phase_separation/HiC/loops/heatmaps/heatmap_loops_control_chr1_zoom.pdf", width = 8, height = 8)
plot_list[[1]]
dev.off()

pdf("/dataOS/rcalandrelli/phase_separation/HiC/loops/heatmaps/heatmap_loops_control_chr1_zoom_99.pdf", width = 8, height = 8)
p2
dev.off()


############ Loop domains
plot_list = list()
k = 1
for (my_sample in all_HiC_samples[1]){
  
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_", my_sample, "_merged")
  
  my_cutoffs = c(0,0.99)
  my_bin_size = 10000
  
  p1 <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/",my_sample,"_merged/", my_bin_size,"/chr1_chr1_",my_bin_size,".txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=10000,
    select_coord=c(44500000,45500000,44500000,45500000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = paste0("/dataOS/rcalandrelli/phase_separation/HiC/tads/", my_sample, "_merged"),
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/", my_sample, "_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1

}

pdf("/dataOS/rcalandrelli/phase_separation/HiC/loops/heatmaps/heatmap_loop_domains_control_chr1_zoom_99.pdf", width = 8, height = 8)
p1
dev.off()


###### oe maps

i = "chr21"
plot_list = list()
k = 1

for (my_sample in c(all_HiC_samples,"H1_RNaseCtrl")){
  
  my_bin_size = 1000000
  
  p <- make_complex_heatmap_HiC(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=paste0("HiC_", my_sample),
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/oe_KR/",my_sample,"_merged/", my_bin_size,"/",i,"_",i,"_",my_bin_size,".txt"),
    sparse=T,
    oe = T,
    chr_row=i,
    chr_col=i,
    bin_size=my_bin_size,
    select_coord=NA,
    contacts_cutoffs=NA,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = NA,
    EP_loops_file = NA,
    feature_color = "blue")
  
  plot_list[[k]] = p
  k = k + 1
}

plot_list_grob = lapply(plot_list, as.grob)
pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/heatmaps/",my_bin_size,"/",i,"_oe.pdf"), width = 20, height = 4)
plot_grid(plotlist = plot_list_grob, nrow = 1)
dev.off()
  


i = "chr21"
my_bin_size = 1000000
my_sample = "H1_RNaseCtrl"

p <- make_complex_heatmap_HiC(
  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=paste0("HiC_", my_sample),
  contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/oe_KR/",my_sample,"_merged/", my_bin_size,"/",i,"_",i,"_",my_bin_size,".txt"),
  sparse=T,
  oe = T,
  chr_row=i,
  chr_col=i,
  bin_size=my_bin_size,
  select_coord=NA,
  contacts_cutoffs=NA,
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = NA,
  loops_file_path = NA,
  EP_loops_file = NA,
  feature_color = "blue")


pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/heatmaps/",my_bin_size,"/",i,"_RNaseCtrl_oe.pdf"), width = 4, height = 4)
p
dev.off()











###################################################################
#### Interaction frequency curve
extract_contact_frequency <- function(directory, 
                                      resolution, # resolution of the contact matrix used
                                      max_distance, # maximum distance for calculation of interaction frequency
                                      n_points # number of points in the curve
){
  contact_frequency_list = list()
  start_bin_distance = 1
  end_bin_distance = max_distance / resolution
  step_bin_distance = end_bin_distance / n_points
  
  options(scipen=999)
  ########## Load genome data
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  for (i in hg38_chromosomes){
    
    ##########  Tile genome and make genomic windows for selected chromosomes
    genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = resolution, cut.last.tile.in.chrom = T)
    chr_row_genome_window = genome_window[seqnames(genome_window)==i]
    chr_col_genome_window = genome_window[seqnames(genome_window)==i]
    
    ########## Load contact matrix
    temp = data.frame(fread(paste0(directory, "/", as.character(resolution), "/", i,"_",i,"_",as.character(resolution),".txt")))
    contact_matrix = matrix(0, nrow = length(chr_row_genome_window), ncol = length(chr_col_genome_window))
    contact_matrix[as.matrix(temp[,c("V1","V2")] / resolution + 1)] = temp$V3
    contact_matrix[as.matrix(temp[,c("V2","V1")] / resolution + 1)] = temp$V3

    rownames(contact_matrix) = seq(1,nrow(contact_matrix)) # ???
    contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
    
    #contact_matrix = as.matrix(read.table(paste0(directory,i,"_",i,"_",as.character(resolution),".txt"), stringsAsFactors = F))
    contact_frequency = c()
    for (bin_distance in seq(start_bin_distance, end_bin_distance, step_bin_distance)){
      k = 1
      contact_frequency_bin_distance = 0
      while(k + bin_distance < nrow(contact_matrix)){
        contact_frequency_bin_distance = contact_frequency_bin_distance + contact_matrix[k,k+bin_distance]
        k = k + 1
      }
      contact_frequency = c(contact_frequency, contact_frequency_bin_distance / k)
    }
    contact_frequency_list[[i]] = contact_frequency
  }
  return(contact_frequency_list)
}

### Control
temp = extract_contact_frequency(directory="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/H1_control_merged/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=100)
# contact_frequency_control = colMeans(do.call("rbind", temp))
contact_frequency_control_20 = colMeans(do.call("rbind", temp))
contact_frequency_control_100 = colMeans(do.call("rbind", temp))


### NH4OAc
temp = extract_contact_frequency(directory="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/H1_NH4OAc_merged/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=100)
# contact_frequency_NH4OAc = colMeans(do.call("rbind", temp))
contact_frequency_NH4OAc_20 = colMeans(do.call("rbind", temp))
contact_frequency_NH4OAc_100 = colMeans(do.call("rbind", temp))

### FL
temp = extract_contact_frequency(directory="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/H1_FL_merged/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=100)
# contact_frequency_FL = colMeans(do.call("rbind", temp))
contact_frequency_FL_20 = colMeans(do.call("rbind", temp))
contact_frequency_FL_100 = colMeans(do.call("rbind", temp))


### RNase
temp = extract_contact_frequency(directory="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/H1_RNase_merged/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=100)
# contact_frequency_RNase = colMeans(do.call("rbind", temp))
contact_frequency_RNase_20 = colMeans(do.call("rbind", temp))
contact_frequency_RNase_100 = colMeans(do.call("rbind", temp))



### Plot
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

write.table(temp_mat,"/dataOS/rcalandrelli/phase_separation/HiC/HiC_distance_interaction_plot.txt", row.names = F, col.names = T, sep = "\t", quote = F)


library(scales)
png("/dataOS/rcalandrelli/phase_separation/HiC/HiC_distance_interaction_plot.png", height = 5, width = 5, units = "in", res = 200)
ggplot(temp_mat, aes(x=distance, y=contact, color=Sample)) + 
  geom_line(size=0.8) +
  scale_x_continuous(limits = c(50000,20000000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.5,1000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Genomic Distance (bp)", y="Interaction Frequency") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks()
dev.off()




