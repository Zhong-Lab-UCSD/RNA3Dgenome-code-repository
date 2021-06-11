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
library(GenomicFeatures)

##################################################################################################
make_complex_heatmap_margi <- function(
  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample="test", # sample name
  contact_matrix_file, # contact matrix 
  sparse=T, # T if the contact matrix is sparse, F if dense
  chr_row, # chr on the rows (RNA end)
  chr_col, # chr on the columns (DNA end)
  bin_size, # bin size or resolution
  select_coord=NA, # to select a portion of the matrix input c(start_row_coord, end_row_coord, start_col_coord, end_col_coord)
  contacts_cutoffs=c(0,NA), # minimum and maximum contact value cutoffs. Leave NA to do not have any max cutoff
  log_contacts=F, # to plot the log of the contacts
  my_colorbar=c("white","red"), # colorbar colors
  Gr_row_annotation=NA, # Granges object for annotation track on the rows
  Gr_col_annotation=NA, # Granges object for annotation track on the columns
  annotation_type="", # set to "genes" to select only protein coding genesß
  Gr_stripes_row=NA, # set True to add horizontal stripes over the heatmap corresponding to the row annotations
  Gr_stripes_col=NA, # set True to add vertical stripes over the heatmap corresponding to the column annotations
  compartments_file_path=NA,
  tads_file_path=NA,
  loops_file_path=NA,
  EP_loops_file=NA,
  feature_color = "#dbdeff",
  color_tads=T, # set color of the features on the heatmap
  normalized=F,
  norm_raw_ratio=F,
  RNA_domain_coord = NA
){
  
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
  
  ##########  Load contact matrix
  if (sparse == T){
    temp = data.frame(fread(contact_matrix_file))
    contact_matrix = matrix(0, nrow = length(chr_row_genome_window), ncol = length(chr_col_genome_window))
    contact_matrix[as.matrix(temp[,c("RNA_end_bins","DNA_end_bins")])] = temp$x
  } else {
    contact_matrix = data.frame(fread(contact_matrix_file))
    contact_matrix = contact_matrix[,-1]
  }
  rownames(contact_matrix) = seq(1,nrow(contact_matrix)) # ???
  
  # Normalize contact matrix by colSums. Each entry is divided by the sum of RNAs attached to that DNA locus.
  if (normalized == T){
    contact_matrix = sweep(contact_matrix, 2, colSums(contact_matrix), FUN = '/') * 100
    contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
  }
  
  if (norm_raw_ratio == T){
    contact_matrix_norm = sweep(contact_matrix, 2, colSums(contact_matrix), FUN = '/')
    contact_matrix_norm[which(is.na(contact_matrix_norm), arr.ind = T)] = 0
    contact_matrix_norm_scale = (contact_matrix_norm - min(contact_matrix_norm)) / (max(contact_matrix_norm)-min(contact_matrix_norm))

    contact_matrix_scale = (contact_matrix - min(contact_matrix)) / (max(contact_matrix)-min(contact_matrix))
    
    contact_matrix_ratio = contact_matrix_norm_scale / contact_matrix_scale
    contact_matrix_ratio[which(is.na(contact_matrix_ratio), arr.ind = T)] = 0
    
    contact_matrix_ratio = log2(contact_matrix_ratio)
    contact_matrix_ratio[which(is.infinite(contact_matrix_ratio), arr.ind = T)] = 0
    contact_matrix = contact_matrix_ratio
  }
  
  ##########  Load compartment file if declared
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
  
  
  ##########  Load TAD file if declared
  if (!is.na(tads_file_path)){
    
    # tads_file = read.table(paste0(tads_file_path,"/TADs_",chr_row,"/50000_blocks.bedpe"))[,c(1:3)]
    # Gr_tads = GRanges(
    #   seqnames = Rle(paste0("chr",tads_file$V1)),
    #   ranges = IRanges(as.numeric(tads_file$V2), end = as.numeric(tads_file$V3), names = c(1:nrow(tads_file))),
    #   strand = Rle(strand('*')))
    # Gr_tads = Gr_tads[seqnames(Gr_tads) == chr_row]
    
    tads_file = read.table(tads_file_path, header = T)
    Gr_tads = GRanges(
      seqnames = Rle(tads_file$seqnames),
      ranges = IRanges(as.numeric(tads_file$start), end = as.numeric(tads_file$end), names = c(1:nrow(tads_file))),
      strand = Rle(strand('*')),
      label = as.character(tads_file$label))
    Gr_tads = Gr_tads[seqnames(Gr_tads) == chr_row]
    
    overlaps = findOverlaps(Gr_tads, chr_row_genome_window)
    
    start_tad_bin = c()
    end_tad_bin = c()
    for (i in 1:length(Gr_tads)){
      temp_tad = subjectHits(overlaps[which(queryHits(overlaps)==i),])
      start_tad_bin = c(start_tad_bin,temp_tad[1])
      end_tad_bin = c(end_tad_bin,temp_tad[length(temp_tad)])
    }
    
    
    # df_tad_bin = data.frame(start_tad_bin, end_tad_bin)
    # for (i in 1:nrow(df_tad_bin)){
    #   start_tad = df_tad_bin[i,1]
    #   end_tad = df_tad_bin[i,2]
    #   contact_matrix[start_tad:end_tad,start_tad] = -2
    #   contact_matrix[start_tad:end_tad,end_tad] = -2
    #   contact_matrix[start_tad,start_tad:end_tad] = -2
    #   contact_matrix[end_tad,start_tad:end_tad] = -2
    #   } 
    
    label_bin = mcols(Gr_tads)[,"label"]
    df_tad_bin = data.frame(start_tad_bin, end_tad_bin, label_bin)

    for (i in 1:nrow(df_tad_bin)){
      # In this way start boundary overlaps the square border, while the end boundary is included within the square. This makes the visualization good
      # reducing the square overlaps
      start_tad = df_tad_bin[i,1]
      end_tad = df_tad_bin[i,2]
      
      if (color_tads == T){
        if (df_tad_bin[i,3] == "enrichment"){
          contact_matrix[start_tad:end_tad,start_tad] = -3
          contact_matrix[start_tad:end_tad,end_tad] = -3
          contact_matrix[start_tad,start_tad:end_tad] = -3
          contact_matrix[end_tad,start_tad:end_tad] = -3
        }
        else if (df_tad_bin[i,3] == "depletion"){
          contact_matrix[start_tad:end_tad,start_tad] = -2
          contact_matrix[start_tad:end_tad,end_tad] = -2
          contact_matrix[start_tad,start_tad:end_tad] = -2
          contact_matrix[end_tad,start_tad:end_tad] = -2
        } else {
          contact_matrix[start_tad:end_tad,start_tad] = -6
          contact_matrix[start_tad:end_tad,end_tad] = -6
          contact_matrix[start_tad,start_tad:end_tad] = -6
          contact_matrix[end_tad,start_tad:end_tad] = -6
        }
      } else {
        if (df_tad_bin[i,3] == "enrichment"){
          contact_matrix[start_tad:end_tad,start_tad] = -6
          contact_matrix[start_tad:end_tad,end_tad] = -6
          contact_matrix[start_tad,start_tad:end_tad] = -6
          contact_matrix[end_tad,start_tad:end_tad] = -6
        }
        else if (df_tad_bin[i,3] == "depletion"){
          contact_matrix[start_tad:end_tad,start_tad] = -6
          contact_matrix[start_tad:end_tad,end_tad] = -6
          contact_matrix[start_tad,start_tad:end_tad] = -6
          contact_matrix[end_tad,start_tad:end_tad] = -6
        } else {
          contact_matrix[start_tad:end_tad,start_tad] = -6
          contact_matrix[start_tad:end_tad,end_tad] = -6
          contact_matrix[start_tad,start_tad:end_tad] = -6
          contact_matrix[end_tad,start_tad:end_tad] = -6
        }
      }
    }
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
      # In this way the dot is exactly over the loop  
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
  
  if (length(nrow(RNA_domain_coord)) > 0){
    RNA_domain_coord = RNA_domain_coord[which(RNA_domain_coord$V1 == chr_row),]
    for (i in seq_len(nrow(RNA_domain_coord))){
      RNA_domain_coord_bin = as.numeric(ceiling(RNA_domain_coord[i,c(2,3,5,6)]/bin_size))
      contact_matrix[RNA_domain_coord_bin[1]:RNA_domain_coord_bin[2], RNA_domain_coord_bin[3]] = -3
      contact_matrix[RNA_domain_coord_bin[1]:RNA_domain_coord_bin[2], RNA_domain_coord_bin[4]] = -3
      contact_matrix[RNA_domain_coord_bin[1], RNA_domain_coord_bin[3]:RNA_domain_coord_bin[4]] = -3
      contact_matrix[RNA_domain_coord_bin[2], RNA_domain_coord_bin[3]:RNA_domain_coord_bin[4]] = -3
    }
  } else if (sum(!is.na(RNA_domain_coord) & length(nrow(RNA_domain_coord)) == 0) == 4) {
      RNA_domain_coord_bin = ceiling(RNA_domain_coord/bin_size)
      contact_matrix[RNA_domain_coord_bin[1]:RNA_domain_coord_bin[2], RNA_domain_coord_bin[3]] = -3
      contact_matrix[RNA_domain_coord_bin[1]:RNA_domain_coord_bin[2], RNA_domain_coord_bin[4]] = -3
      contact_matrix[RNA_domain_coord_bin[1], RNA_domain_coord_bin[3]:RNA_domain_coord_bin[4]] = -3
      contact_matrix[RNA_domain_coord_bin[2], RNA_domain_coord_bin[3]:RNA_domain_coord_bin[4]] = -3
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
  
  
  ##########  Check if making annotation tracks
  if (!is.na(Gr_row_annotation[1])){
    
    if (annotation_type == "genes"){
      Gr_row_annotation = Gr_row_annotation[mcols(Gr_row_annotation)["GENEBIOTYPE"][,1] == "protein_coding"]
    }
    
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
    
    if (annotation_type == "genes"){
      Gr_col_annotation = Gr_col_annotation[mcols(Gr_col_annotation)["GENEBIOTYPE"][,1] == "protein_coding"]
    }
    
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
  
  
  #### Log contacts
  if (log_contacts == T){
    contact_matrix[which(contact_matrix>0, arr.ind = T)] = log(contact_matrix[which(contact_matrix>0, arr.ind = T)])
  }
  
  
  ##########  Set colorbar
  
  if (norm_raw_ratio == F){
    min_cutoff = as.numeric(contacts_cutoffs[1])
    if (is.na(contacts_cutoffs[2])){
      max_cutoff = max(contact_matrix)
    } else if (contacts_cutoffs[2] < 1){
      max_cutoff = round(as.numeric(quantile(contact_matrix[contact_matrix > 0], contacts_cutoffs[2])))
    } else if (contacts_cutoffs[2] > 1){
      max_cutoff = contacts_cutoffs[2]
    }
    
    ### Custom part to plot Figure STRIPE B
    # contact_matrix[23,13] = 8
    # contact_matrix[13,23] = 0
    # 
    # contact_matrix[39,10] = -2
    # contact_matrix[10,39] = -2
    # 
    # contact_matrix[37,12] = -4
    # contact_matrix[12,37] = -4
    # 
    # contact_matrix[38,25] = -5
    # contact_matrix[25,38] = -5
    
    col_fun1 = colorRamp2(c(min_cutoff, max_cutoff), my_colorbar)
    col_fun = function(x) {
      ifelse(x == -1, "#dbdeff",
             ifelse(x == -2, "green",  # depleted TADs / enhancer_promoter loops
                    ifelse(x == -3, "blue", # enriched TADs / loops in general
                           ifelse(x == -4, "#ff00ff", # promoter_enhancer loops
                                  ifelse(x == -5, "#000000", # both promoter_enhancer and vice versa loops
                                         ifelse(x == -6, "#b3b3b3", # neither tads
                                                col_fun1(x)))))))
    }
    
    attr(col_fun, "breaks") = c(min_cutoff, max_cutoff)
 
   } else {
    
    col_fun = colorRamp2(c(min(contact_matrix), 0, max(contact_matrix)), c("blue","white","red"))
    # col_fun = function(x) {
    #   col_fun1(x)
    # }
    # 
    # attr(col_fun, "breaks") = c(min(contact_matrix), 0, max(contact_matrix))
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

chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes"
sample="test"
contact_matrix_file="/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_control/500000_filter1k/chr11_chr11.txt"
sparse=T
chr_row="chr11"
chr_col="chr11"
bin_size=500000
select_coord=NA
contacts_cutoffs=c(0,NA)
my_colorbar=c("white","red")
Gr_row_annotation=NA
Gr_col_annotation=NA
Gr_stripes_row=NA
Gr_stripes_col=NA
compartments_file_path="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt"
tads_file_path = NA #"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged/")
loops_file_path = NA #"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged"
feature_color = "#dbdeff"


my_sample = "iMARGI_H1_control_1"
all_imargi_samples = list.dirs(path = "/dataOS/rcalandrelli/phase_separation/MARGI/data", full.names = F, recursive = F)



############ Compartments heatmap plots with the same colorbar scale as control (0, 800)
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = gsub("iMARGI_","",my_sample)
  my_cutoffs = c(0,400)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/500000_filter1k/chr1_chr1.txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=500000,
    select_coord=NA,
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    compartments_file_path=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/",temp,"_500000.txt"),
    tads_file_path = NA,
    loops_file_path = NA,
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
}

plot_list_grob = lapply(plot_list, as.grob)
#png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/heatmaps/chr11_margi_compartment_heatmaps_800.png", width = 12, height = 12, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/heatmaps/chr1_margi_compartment_heatmaps_400.png", width = 12, height = 12, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()



############ Compartments heatmap plots
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = gsub("iMARGI_","",my_sample)
  my_cutoffs = c(0,NA)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/500000_filter1k/chr11_chr11.txt"),
    sparse=T,
    chr_row="chr11",
    chr_col="chr11",
    bin_size=500000,
    select_coord=NA,
    log_contacts = T,
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    compartments_file_path=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/",temp,"_500000.txt"),
    tads_file_path = NA,
    loops_file_path = NA,
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
}

plot_list_grob = lapply(plot_list, as.grob)

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/heatmaps/chr11_margi_compartment_heatmaps_995.png", width = 12, height = 12, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/heatmaps/chr11_margi_compartment_heatmaps_99.png", width = 12, height = 12, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/heatmaps/chr11_margi_compartment_heatmaps_999.png", width = 12, height = 12, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/heatmaps/chr11_margi_compartment_heatmaps_log.png", width = 12, height = 12, units = "in", res = 200)

plot_grid(plotlist = plot_list_grob, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()




############ TADs heatmap plots
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  my_cutoffs = c(0,0.995)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/100000_filter1k/chr1_chr1.txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=100000,
    select_coord=NA,
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/MARGI/tads/iMARGI_H1_control_1/10000_tads_enrich_label.txt",#"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/dataOS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
    loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
    feature_color = "blue")#dbdeff")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/100000_filter1k/chr1_chr1.txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=100000,
    select_coord=c(72000000,84000000,72000000,84000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label.txt"),#"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/dataOS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
    loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
    feature_color = "blue")#dbdeff")
  
  plot_list[[k]] = p2
  k = k + 1
  
  p3 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/50000_filter1k/chr1_chr1.txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=50000,
    select_coord=rep(c(153000000,161000000,),2),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/",my_sample,"/10000_tads_RNA_enrich_label.txt"),#"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/dataOS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
    loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
    feature_color = "blue")#dbdeff")

  plot_list[[k]] = p3
  k = k + 1

}

plot_list_grob = lapply(plot_list, as.grob)
#png("/dataOS/rcalandrelli/phase_separation/MARGI/tads/heatmaps/chr1_margi_tad_heatmaps_FULL.png", width = 12, height = 16, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/tads/heatmaps/chr1_margi_tad_heatmaps_995.png", width = 12, height = 16, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 3, labels = c("A",rep("",2), "B",rep("",2), "C",rep("",2), "D",rep("",2)))
dev.off()


############ Loops heatmap plots
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  my_cutoffs = c(0,0.995)
  
  my_chr = "chr1"
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/",my_chr,"_",my_chr,".txt"),
    sparse=T,
    chr_row=my_chr,
    chr_col=my_chr,
    bin_size=10000,
    select_coord=rep(c(3000000,4000000),2),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1

  p2 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
    sparse=T,
    chr_row="chr20",
    chr_col="chr20",
    bin_size=10000,
    select_coord=c(37000000,38000000,37000000,38000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    feature_color = "blue")

  plot_list[[k]] = p2
  k = k + 1
  
  # p3 <- make_complex_heatmap_margi(
  #   chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  #   sample=my_sample,
  #   contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr3_chr3.txt"),
  #   sparse=T,
  #   chr_row="chr3",
  #   chr_col="chr3",
  #   bin_size=10000,
  #   select_coord=c(47500000,48500000,47500000,48500000),
  #   contacts_cutoffs=my_cutoffs,
  #   my_colorbar=c("white","red"),
  #   Gr_row_annotation=NA,
  #   Gr_col_annotation=NA,
  #   Gr_stripes_row=NA,
  #   Gr_stripes_col=NA,
  #   tads_file_path = NA,
  #   loops_file_path = loop_file_path_temp,
  #   feature_color = "blue")
  # 
  # plot_list[[k]] = p3
  # k = k + 1
  
}

plot_list_grob = lapply(plot_list, as.grob)
#png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_loops_FULL.png", width = 8, height = 16, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_loops_995.png", width = 8, height = 16, units = "in", res = 200)
#plot_grid(plotlist = plot_list_grob, ncol = 3, labels = c("A","B","C","D","","","E","","","F","",""))
plot_grid(plotlist = plot_list_grob, ncol = 2, labels = c("A","B","C","","D","","E",""))
dev.off()


#### test to find the best looking example
temp_loops = loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio > 3.7 & 
                                      !is.na(loops_control_EP_ratio$density_ratio)),]
rownames(temp_loops) = c(1:nrow(temp_loops))

for (i in 235:nrow(temp_loops)){

  print(i)
  
  my_start_row = round_any(temp_loops[i,2], 1000000)
  if (my_start_row > temp_loops[i,2]){
    my_start_row = my_start_row - 1000000
  }
  my_end_row = my_start_row + 1000000
  
  my_start_col = round_any(temp_loops[i,5], 1000000)
  if (my_start_col > temp_loops[i,5]){
    my_start_col = my_start_col - 1000000
  }
  my_end_col = my_start_col + 1000000
  
  my_start = min(my_start_row, my_start_col)
  my_end = max(my_end_row, my_end_col)
  
  plot_list = list()
  k = 1
  
  for (my_sample in all_imargi_samples[1:4]){
    
    temp = strsplit(my_sample,"_")[[1]][3]
    loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
    
    my_cutoffs = c(0,0.995)
    
    my_chr = paste0("chr", temp_loops[i,1])
    
    p1 <- make_complex_heatmap_margi(
      chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
      sample=my_sample,
      contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/",my_chr,"_",my_chr,".txt"),
      sparse=T,
      chr_row=my_chr,
      chr_col=my_chr,
      bin_size=10000,
      #select_coord=c(my_start, my_end, my_start, my_end),
      select_coord=rep(c(46500000, 47500000),2),
      contacts_cutoffs=my_cutoffs,
      my_colorbar=c("white","red"),
      Gr_row_annotation=NA,
      Gr_col_annotation=NA,
      Gr_stripes_row=NA,
      Gr_stripes_col=NA,
      tads_file_path = NA,
      loops_file_path = loop_file_path_temp,
      feature_color = "blue")
    
    plot_list[[k]] = p1
    k = k + 1
    
  }
  
  plot_list_grob = lapply(plot_list, as.grob)
  a <- plot_grid(plotlist = plot_list_grob, ncol = 2)
  png(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/pattern_1_all_plots/", i, "_heatmaps_loops_995.png"), width = 8, height = 8, units = "in", res = 200)
  print(a)
  dev.off()

}






############ Loops heatmap plots to match Figure S-loop-examples
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  my_cutoffs = c(0,0.995)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr12_chr12.txt"),
    sparse=T,
    chr_row="chr12",
    chr_col="chr12",
    bin_size=10000,
    select_coord=c(52500000,54000000,52500000,54000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = NA,
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
}
  
plot_list_grob = lapply(plot_list, as.grob)
png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_margi_S-loop-examples1.png", width = 8, height = 8, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 2, labels = c("A","B","C","D"))
dev.off()


############ Loops (with marker enhancer-promoter) heatmap plots
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  my_cutoffs = c(0,0.995)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr1_chr1.txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=10000,
    select_coord=c(3000000,4000000,3000000,4000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
    sparse=T,
    chr_row="chr20",
    chr_col="chr20",
    bin_size=10000,
    select_coord=c(37000000,38000000,37000000,38000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p2
  k = k + 1
  
}

plot_list_grob = lapply(plot_list, as.grob)
#png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_loops_FULL.png", width = 8, height = 16, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_EP_loops_995.png", width = 8, height = 16, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 2, labels = c("A","B","C","","D","","E",""))
dev.off()




############ Loops heatmap plots (without superimposed loops)
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1:4]){
  
  temp = strsplit(my_sample,"_")[[1]][3]

  my_cutoffs = c(0,0.995)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr1_chr1.txt"),
    sparse=T,
    chr_row="chr1",
    chr_col="chr1",
    bin_size=10000,
    select_coord=c(3000000,4000000,3000000,4000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = NA,
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr20_chr20.txt"),
    sparse=T,
    chr_row="chr20",
    chr_col="chr20",
    bin_size=10000,
    select_coord=c(37000000,38000000,37000000,38000000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = NA,
    feature_color = "blue")
  
  plot_list[[k]] = p2
  k = k + 1
  
  # p3 <- make_complex_heatmap_margi(
  #   chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  #   sample=my_sample,
  #   contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr3_chr3.txt"),
  #   sparse=T,
  #   chr_row="chr3",
  #   chr_col="chr3",
  #   bin_size=10000,
  #   select_coord=c(47500000,48500000,47500000,48500000),
  #   contacts_cutoffs=my_cutoffs,
  #   my_colorbar=c("white","red"),
  #   Gr_row_annotation=NA,
  #   Gr_col_annotation=NA,
  #   Gr_stripes_row=NA,
  #   Gr_stripes_col=NA,
  #   tads_file_path = NA,
  #   loops_file_path = NA,
  #   feature_color = "blue")
  # 
  # plot_list[[k]] = p3
  # k = k + 1
  
}

plot_list_grob = lapply(plot_list, as.grob)
#png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_loops_FULL_noLoops.png", width = 8, height = 16, units = "in", res = 200)
png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmaps_loops_995_noLoops.png", width = 8, height = 16, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 2, labels = c("A","B","C","","D","","E",""))
dev.off()



temp = "control"
my_cutoffs = c(0,0.995)
my_sample = "iMARGI_H1_control"

p1 <- make_complex_heatmap_margi(
  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr1_chr1.txt"),
  sparse=T,
  chr_row="chr1",
  chr_col="chr1",
  bin_size=10000,
  select_coord=c(10000000,20000000,10000000,20000000),
  contacts_cutoffs=my_cutoffs,
  my_colorbar=c("white","red"),
  Gr_row_annotation=NA,
  Gr_col_annotation=NA,
  Gr_stripes_row=NA,
  Gr_stripes_col=NA,
  tads_file_path = NA,
  loops_file_path = NA,
  feature_color = "blue")





############ Loops heatmap plots zoomed over patterns (only control)
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  my_cutoffs = c(0,0.995)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k_RNA_from_gene/chr1_chr1.txt"),
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
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k_RNA_from_gene/chr20_chr20.txt"),
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
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p2
  k = k + 1
}

#plot_list_grob = lapply(plot_list, as.grob)

pdf("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmap_loops_control_chr1_zoom_995_RNA_from_gene.pdf", width = 8, height = 8)
plot_list[[1]]
dev.off()

pdf("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmap_loops_control_chr20_zoom_995_RNA_from_gene.pdf", width = 8, height = 8)
plot_list[[2]]
dev.off()
  

################# More examples with pattern 1
plot_list = list()
k = 1
for (my_sample in all_imargi_samples[1]){
  
  temp = strsplit(my_sample,"_")[[1]][3]
  loop_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_",temp,"_merged")
  
  my_cutoffs = c(0,0.995)
  
  p1 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr10_chr10.txt"),
    sparse=T,
    chr_row="chr10",
    chr_col="chr10",
    bin_size=10000,
    select_coord=c(69600000,70100000,69600000,70100000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p1
  k = k + 1
  
  p2 <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr9_chr9.txt"),
    sparse=T,
    chr_row="chr9",
    chr_col="chr9",
    bin_size=10000,
    select_coord=c(96300000,96800000,96300000,96800000),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = loop_file_path_temp,
    EP_loops_file = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_",temp,"_EP.bedpe"),
    feature_color = "blue")
  
  plot_list[[k]] = p2
  k = k + 1
}

pdf("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmap_loops_control_chr10_zoom_995.pdf", width = 8, height = 8)
p1
dev.off()

pdf("/dataOS/rcalandrelli/phase_separation/MARGI/loops/heatmaps/heatmap_loops_control_chr9_zoom_995.pdf", width = 8, height = 8)
p2
dev.off()


############ Figure S-RNAdomain
my_sample = "iMARGI_HFF_control"
loop_file_path_temp = "/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_HFF/"
my_cutoffs = c(0,0.995)
HFF_stripes = read.table("/dataOS/qizhijie/iMARGI_project/stripeCalling/stripe_cellLine_narrow/HFF.raw.2dBed")
  
temp_start = 55000000
temp_end = 75000000
  
temp_stripes = HFF_stripes[which(HFF_stripes$V1 == "chr2" &
                                 HFF_stripes$V2 > temp_start &
                                   HFF_stripes$V5 > temp_s)]
  
p <- make_complex_heatmap_margi(
    chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
    sample=my_sample,
    contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/50000_filter200k/chr2_chr2.txt"),
    sparse=T,
    chr_row="chr2",
    chr_col="chr2",
    bin_size=50000,
    select_coord=rep(c(temp_start,temp_end),2),
    contacts_cutoffs=my_cutoffs,
    my_colorbar=c("white","red"),
    Gr_row_annotation=NA,
    Gr_col_annotation=NA,
    Gr_stripes_row=NA,
    Gr_stripes_col=NA,
    tads_file_path = NA,
    loops_file_path = NA,
    EP_loops_file = NA,
    feature_color = "blue",
    RNA_domain_coord = HFF_stripes)
  

pdf("/dataOS/rcalandrelli/phase_separation/paper_plots/Figure_S_RNAdomain_A1.pdf", width = 8, height = 8)
print(p)
dev.off()

p <- make_complex_heatmap_margi(
  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  sample=my_sample,
  contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/10000_filter1k/chr6_chr6.txt"),
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

pdf("/dataOS/rcalandrelli/phase_separation/paper_plots/Figure_S_RNAdomain_A2.pdf", width = 8, height = 8)
print(p)
dev.off()


############################### Normalized maps


for (my_chr in hg38_chromosomes){
  
  my_cutoffs = c(0,0.995)
  bin_size = 100000
  plot_list = list()
  k = 1
  
  for (my_sample in all_imargi_samples[1:4]){

    p1 <- make_complex_heatmap_margi(
      chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
      sample=my_sample,
      contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/",bin_size,"_filter1k/",my_chr,"_",my_chr,".txt"),
      sparse=T,
      chr_row=my_chr,
      chr_col=my_chr,
      bin_size=bin_size,
      select_coord=c(0,NA),#rep(c(37000000,38000000),2),
      contacts_cutoffs=my_cutoffs,
      my_colorbar=c("white","red"),
      Gr_row_annotation=NA,
      Gr_col_annotation=NA,
      Gr_stripes_row=NA,
      Gr_stripes_col=NA,
      tads_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/MARGI/tads/iMARGI_H1_control_1/10000_tads_enrich_label.txt",#"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/dataOS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
      loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
      feature_color = "blue")#dbdeff")
    
    plot_list[[k]] = p1
    k = k + 1
    
    p2 <- make_complex_heatmap_margi(
      chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
      sample=my_sample,
      contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/",bin_size,"_filter1k/",my_chr,"_",my_chr,".txt"),
      sparse=T,
      chr_row=my_chr,
      chr_col=my_chr,
      bin_size=bin_size,
      select_coord=c(0,NA),#rep(c(37000000,38000000),2),
      contacts_cutoffs=my_cutoffs,
      my_colorbar=c("white","red"),
      Gr_row_annotation=NA,
      Gr_col_annotation=NA,
      Gr_stripes_row=NA,
      Gr_stripes_col=NA,
      tads_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/MARGI/tads/iMARGI_H1_control_1/10000_tads_enrich_label.txt",#"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/dataOS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
      loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
      feature_color = "blue",
      normalized = T)
    
    plot_list[[k]] = p2
    k = k + 1
    
    p3 <- make_complex_heatmap_margi(
      chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
      sample=my_sample,
      contact_matrix_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",my_sample,"/",bin_size,"_filter1k/",my_chr,"_",my_chr,".txt"),
      sparse=T,
      chr_row=my_chr,
      chr_col=my_chr,
      bin_size=bin_size,
      select_coord=c(0,NA),#rep(c(37000000,38000000),2),
      contacts_cutoffs=c(0,NA),
      my_colorbar=c("white","red"),
      Gr_row_annotation=NA,
      Gr_col_annotation=NA,
      Gr_stripes_row=NA,
      Gr_stripes_col=NA,
      tads_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/MARGI/tads/iMARGI_H1_control_1/10000_tads_enrich_label.txt",#"/dataOS/rcalandrelli/phase_separation/TADs/H1_control_merged",#/dataOS/rcalandrelli/phase_separation/MARGI/tads/tads_enrich_label.txt",
      loops_file_path = NA,#"/dataOS/rcalandrelli/phase_separation/loops/hiccups_H1_control_merged",
      feature_color = "blue",
      normalized = F,
      norm_raw_ratio = T)
    
    plot_list[[k]] = p3
    k = k + 1
    
  }
  
  plot_list_grob = lapply(plot_list, as.grob)
  p <- plot_grid(plotlist = plot_list_grob, ncol = 3)
  png(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/normalized_heatmaps/",my_chr,"_100k_995.png"), width = 12, height = 16, units = "in", res = 200)
  print(p)
  dev.off()
 
}

plot_list_grob = lapply(plot_list, as.grob)
png("/dataOS/rcalandrelli/phase_separation/MARGI/normalized_heatmaps/chr20_10k_995.png", width = 8, height = 16, units = "in", res = 200)
plot_grid(plotlist = plot_list_grob, ncol = 2)
dev.off()







######################################## Annotations
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


chrSize = read.table("/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes")
hg38_chromosomes = as.character(chrSize$V1)
hg38_lengths = as.numeric(chrSize$V2)
names(hg38_lengths) = hg38_chromosomes

# annotations_orgHsDb <- AnnotationDbi::select(org.Hs.eg.db, # database
#                                              keys = keys(org.Hs.eg.db),  # data to use for retrieval
#                                              columns = c("SYMBOL", "ENTREZID","GENENAME"), # information to retrieve for given data
#                                              keytype = "ENTREZID") # type of data given in 'keys' argument
# 
# annotations_edb_hg38 <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
#                                               keys = keys(EnsDb.Hsapiens.v86),
#                                               columns = c("SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND", "SYMBOL", "ENTREZID","GENEBIOTYPE"),
#                                               keytype = "GENEID")
# 
# annotations_edb_hg38$SEQNAME = paste0("chr", annotations_edb_hg38$SEQNAME)
# annotations_edb_hg38 = annotations_edb_hg38[which(annotations_edb_hg38$SEQNAME %in% hg38_chromosomes),]
# 
Gr_annotations_edb_hg38 = GRanges(
  seqnames = Rle(annotations_edb_hg38$SEQNAME),
  ranges = IRanges(annotations_edb_hg38$GENESEQSTART, end = annotations_edb_hg38$GENESEQEND, names = c(1:nrow(annotations_edb_hg38))),
  strand = Rle(strand('*')),
  SYMBOL = annotations_edb_hg38$SYMBOL,
  GENEBIOTYPE = annotations_edb_hg38$GENEBIOTYPE)
# 
# overlaps = countOverlaps(Gr_annotations_edb_hg38, chr1.region)
# data.frame(Gr_annotations_edb_hg38[overlaps>0])

TxDb.Hsapiens.encode.GRCh38.84.knownGene = makeTxDbFromGFF('/dataOS/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf')
TxDb.Hsapiens.gencode.24 = makeTxDbFromGFF('/home/frankyan/research/refGenome/standard4DN/GRCh38/gencode.v24.primary_assembly.annotation.gtf')

H3K27ac_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF423TVA_H3K27ac.bw"
H3K4me1_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF584AVI_H3K4me1.bw"
H3K4me3_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF422PZQ_H3K36me3.bw"
CTCF_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF473IZV_CTCF.bw"
RAD21_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF913JGA_RAD21.bw"
HDAC2_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF640QBB_HDAC2.bw"
H2AFZ_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF775ZWT_H2AFZ.bw"
YY1_chipseq <- "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/ENCFF406PYH_YY1.bw"

chipseq_tracks_list = list(H3K27ac_chipseq, H3K4me3_chipseq, H3K4me1_chipseq, CTCF_chipseq, RAD21_chipseq, YY1_chipseq)
chipseq_tracks = c("H3K27ac", "H3K4me3", "H3K4me1", "CTCF", "RAD21", "YY1")
names(chipseq_tracks_list) = chipseq_tracks

make_annotation_tracks <- function(chr,
                                   start,
                                   end,
                                   sample_name,
                                   panel){
  
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
  
  
  a = which(mcols(genes.data[["genes"]])["gene_id"] %in% protein_coding_genes)[[1]]
  genes.data[["genes"]] = genes.data[["genes"]][a]
  
  genes.data <- mergeTranscripts(genes.data)
  
  gene_track_border = 0.75
  step = 0.8 / length(chipseq_tracks_list) # height of each track
  track_borders = seq(0.8,0, by=-step)
  
  #png(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_name, "_", chr, "_", start, "_", end, "_tracks.png"), width = 8.5, height = 6, units = "in", res = 500)
  
  if (panel == "A"){
    pdf(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_name, "_", chr,"_", start, "_", end, "_tracks.pdf"), width = 8.5, height = 6)
  } 
  else if (panel == "B"){
    pdf(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_name, "_", chr,"_", start, "_", end, "_tracks.pdf"), width = 8.5, height = 6)
  }
  
  plot_params <- getDefaultPlotParams(plot.type=1)
  plot_params$data1inmargin <- 5
  plot_params$ideogramheight <- 10
  plot_params$bottommargin <- 30
  plot_params$topmargin <- 0
  plot_params$leftmargin <- 0.15
  plot_params$rightmargin <- 0.02
  
  kp <- plotKaryotype(zoom = chr.region, plot.params = plot_params)
  
  if (panel == "A"){
    kpPlotGenes(kp, data=genes.data, r1=1, r0=0.8, gene.name.position = "left", gene.name.cex = 1, mark.height = 0.3)
  } 
  else if (panel == "B"){
    kpPlotGenes(kp, data=genes.data, r1=1, r0=0.7, gene.name.position = "left", gene.name.cex = 1.2, mark.height = 0.4) # Figure LOOPD G-J
  }
  
  # for(i in 1:length(chipseq_tracks_list)){
  #   kp <- kpPlotBigWig(kp, data=chipseq_tracks_list[[i]], ymax="visible.region", r1 = track_borders[i]-step/3, r0 = track_borders[i+1])
  #   computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  #   kpAxis(kp, ymin=0, ymax=computed.ymax, r1 = track_borders[i]-step/3, r0 = track_borders[i+1], cex=0.8)
  #   kpAddLabels(kp, labels = chipseq_tracks[i], r1=track_borders[i]-step/2, r0=track_borders[i]-step/2, cex=1, srt = 0, offset=5, label.margin=0.06)
  # }
  
  dev.off()
  
}


make_annotation_tracks(chr="chr1", start=3000000, end=4000000, sample_name="control")
make_annotation_tracks(chr="chr20", start=37000000, end=38000000, sample_name="control")


make_annotation_tracks(chr="chr20", start=37150000, end=37650000, sample_name="control")

make_annotation_tracks(chr="chr9", start=96300000, end=96800000, sample_name="control")
make_annotation_tracks(chr="chr10", start=69600000, end=70100000, sample_name="control")


make_annotation_tracks(chr="chr6", start=15500000, end=17500000, sample_name="HFF", panel="A")


######################################## Interactions
library(chicane)

hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY')) # UCSC

TxDb.Hsapiens.gencode.24 = makeTxDbFromGFF('/home/frankyan/research/refGenome/standard4DN/GRCh38/gencode.v24.primary_assembly.annotation.gtf')
gencode.24_genes = genes(TxDb.Hsapiens.gencode.24)
gencode.24_genes = gencode.24_genes[seqnames(gencode.24_genes) %in% hg38_chromosomes]

make_margi_interaction_track <- function(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                   margi_bedpe_file,
                                   chr_plot,
                                   start_plot,
                                   end_plot,
                                   only_RNA_starting_from_gene=F,
                                   sample_label
){
  
  ########## Load genome data
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ### Tile genome and make genomic windows for selected chromosomes
  genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = 10000, cut.last.tile.in.chrom = T)
  #genome_window <- genome_window[seqnames(genome_window) == chr_plot]
  
  ########## Loading input margi data for selecting interactions
  margi_bedpe = data.frame(fread(margi_bedpe_file))
  margi_bedpe = margi_bedpe[which(margi_bedpe$X.chrom1 == margi_bedpe$chrom2),]
  margi_bedpe = margi_bedpe[which(margi_bedpe$X.chrom1 == chr_plot),]
  
  Gr_margi_RNA = GRanges(
    seqnames = Rle(margi_bedpe$X.chrom1),
    ranges = IRanges(margi_bedpe$start1+1, end = margi_bedpe$end1, names = c(1:nrow(margi_bedpe))),
    strand = Rle(strand('*')))
  
  Gr_margi_DNA = GRanges(
    seqnames = Rle(margi_bedpe$chrom2),
    ranges = IRanges(margi_bedpe$start2+1, end = margi_bedpe$end2, names = c(1:nrow(margi_bedpe))),
    strand = Rle(strand('*')))
  
  if (only_RNA_starting_from_gene == T){
    
    # Make Granges with start coordinate of each RNA end read
    Gr_margi_RNA_start = GRanges(
      seqnames = Rle(margi_bedpe$X.chrom1),
      ranges = IRanges(margi_bedpe$start1+1, end = margi_bedpe$start1+2, names = c(1:nrow(margi_bedpe))),
      strand = Rle(strand('*')))
    
    # Overlap with genes
    overlaps = countOverlaps(Gr_margi_RNA_start, gencode.24_genes, ignore.strand = T)
    
    Gr_margi_RNA = Gr_margi_RNA[overlaps > 0]
    Gr_margi_DNA = Gr_margi_DNA[overlaps > 0]
    
    names(Gr_margi_RNA) = seq(1:length(Gr_margi_RNA))
    names(Gr_margi_DNA) = seq(1:length(Gr_margi_DNA))
    
  }
  
  ########## Selecting ranges corresponding to interaction region based also on pattern type
  interaction_region = GRanges(paste0(chr_plot,":",start_plot,"-",end_plot))
      
  overlaps_RNA = countOverlaps(Gr_margi_RNA, interaction_region, ignore.strand = T)
  Gr_margi_RNA_overlaps = Gr_margi_RNA[overlaps_RNA > 0]
      
  overlaps_DNA = countOverlaps(Gr_margi_DNA, interaction_region, ignore.strand = T)
  Gr_margi_DNA_overlaps = Gr_margi_DNA[overlaps_DNA > 0]
      
  temp = intersect(names(Gr_margi_RNA_overlaps), names(Gr_margi_DNA_overlaps))
      
  Gr_margi_RNA_overlaps = Gr_margi_RNA_overlaps[temp]
  Gr_margi_DNA_overlaps = Gr_margi_DNA_overlaps[temp]
  
  
  # test
  # interaction_region_RNA = GRanges(paste0(chr_plot,":",3740000,"-",3750000))
  # overlaps_RNA = countOverlaps(Gr_margi_RNA, interaction_region_RNA, ignore.strand = T)
  # Gr_margi_RNA_overlaps = Gr_margi_RNA[overlaps_RNA > 0]
  # 
  # 
  # temp = countOverlaps(Gr_margi_DNA_overlaps, GRanges(paste0(chr_plot,":",3260001,"-",3270000)))
  # temp1 = Gr_margi_DNA_overlaps[temp>0]
  # temp2 = data.frame(Gr_margi_RNA_overlaps[temp>0])
      
  
  ##########  Converting from read pairs to bin size interactions
  overlaps = data.frame(findOverlaps(genome_window, Gr_margi_RNA_overlaps))
  overlaps = overlaps[order(overlaps$subjectHits),]
  temp = overlaps[!duplicated(overlaps$subjectHits),"queryHits"] # to remove cases where a read may overlap two genomic windows
  anchor_RNA = genome_window[temp]
  
  overlaps = data.frame(findOverlaps(genome_window, Gr_margi_DNA_overlaps))
  overlaps = overlaps[order(overlaps$subjectHits),]
  temp = overlaps[!duplicated(overlaps$subjectHits),"queryHits"] # to remove cases where a read may overlap two genomic windows
  anchor_DNA = genome_window[temp]
  
  # test
  # interaction_region_RNA = GRanges(paste0(chr_plot,":",3740001,"-",3750000))
  # overlaps_RNA = countOverlaps(anchor_RNA, interaction_region_RNA, ignore.strand = T)
  # anchor_RNA_overlaps = anchor_RNA[overlaps_RNA > 0]
  # 
  # anchor_DNA_overlaps = anchor_DNA[overlaps_RNA > 0]
  # 
  # temp = countOverlaps(anchor_DNA, GRanges(paste0(chr_plot,":",3260001,"-",3270000)))
  # temp1 = anchor_DNA[temp>0]
  # temp2=data.frame(anchor_RNA[names(temp1)])
  
  
  ########## Making interaction object
  interaction.object <- GenomicInteractions(
    anchor1 = anchor_RNA,
    anchor2 = anchor_DNA
  )
  
  df.interaction.object = data.frame(interaction.object)
  df.interaction.object_agg = ddply(df.interaction.object,.(seqnames1,start1,end1,seqnames2,start2,end2),nrow)
  
  anchor_RNA_agg = GRanges(
    seqnames = Rle(df.interaction.object_agg$seqnames1),
    ranges = IRanges(df.interaction.object_agg$start1, end = df.interaction.object_agg$end1, names = c(1:nrow(df.interaction.object_agg))),
    strand = Rle(strand('*')))
  
  anchor_DNA_agg = GRanges(
    seqnames = Rle(df.interaction.object_agg$seqnames2),
    ranges = IRanges(df.interaction.object_agg$start2, end = df.interaction.object_agg$end2, names = c(1:nrow(df.interaction.object_agg))),
    strand = Rle(strand('*')))
  
  interaction.object.agg <- GenomicInteractions(
    anchor1 = anchor_RNA_agg,
    anchor2 = anchor_DNA_agg,
    count = df.interaction.object_agg$V1
  )
  
  ##### Making plot
  ideogram.track <- IdeogramTrack(genome = 'hg38', chromosome = chr_plot)
  genome.axis.track <- GenomeAxisTrack()
  
  interaction.track <- InteractionTrack(
    interaction.object.agg,
    name = ''
  )
  
  if (only_RNA_starting_from_gene == F){
    out_file = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_label,"_",chr_plot,"_",start_plot,"_",end_plot,"_interactions.png")
  } else {
    out_file = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_label,"_",chr_plot,"_",start_plot,"_",end_plot,"_interactions_RNA_from_gene.png")
  }
  
  png(out_file, width = 8, height = 2, units = "in", res = 500)
  plotTracks(
    list(ideogram.track, interaction.track),
    sizes = c(0.1,0.8),
    from = start_plot,
    to = end_plot,
    col="red",
    col.anchors.line = "black",
    col.anchors.fill = "white",
    cex = 1.5,
    fontcolor = "black",
    background.title = "white"
  )
  dev.off()
  
}

make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                       margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                       chr_plot = "chr1",
                       start_plot=3250000,
                       end_plot=3750000)

make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                       margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                       chr_plot = "chr20",
                       start_plot=37150000,
                       end_plot=37650000)

make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                             chr_plot = "chr1",
                             start_plot=3250000,
                             end_plot=3750000,
                             only_RNA_starting_from_gene = T)

make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                             chr_plot = "chr20",
                             start_plot=37150000,
                             end_plot=37650000,
                             only_RNA_starting_from_gene = T)


make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                             chr_plot = "chr9",
                             start_plot=96300000,
                             end_plot=96800000)

make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe",
                             chr_plot = "chr10",
                             start_plot=69600000,
                             end_plot=70100000)


make_margi_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             margi_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_HFF_control/iMARGI_HFF_control.mapq30.1k.final.bedpe.gz",
                             chr_plot = "chr2",
                             start_plot=149300000,
                             end_plot=151600000,
                             sample_label = "HFF")


##### Loops
make_loop_interaction_track <- function(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                        loop_file,
                                         chr_plot,
                                         start_plot,
                                         end_plot,
                                         sample_label
){
  
  ########## Load genome data
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ### Tile genome and make genomic windows for selected chromosomes
  genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = 10000, cut.last.tile.in.chrom = T)
  #genome_window <- genome_window[seqnames(genome_window) == chr_plot]
  
  ########## Loading input margi data for selecting interactions
  loop_bedpe = read.table(loop_file)[,1:6]
  loop_bedpe$V1 = paste0("chr",loop_bedpe$V1)
  loop_bedpe$V4 = paste0("chr",loop_bedpe$V4)
  loop_bedpe = loop_bedpe[which(loop_bedpe$V1 == chr_plot),]
  
  loop_bedpe = loop_bedpe[which(loop_bedpe$V2 > start_plot &
                                  loop_bedpe$V6 < end_plot),]
  
  

  interaction_region = GRanges(paste0(chr_plot,":",start_plot,"-",end_plot))
  
  Gr_anchor_1 = GRanges(
    seqnames = Rle(loop_bedpe$V1),
    ranges = IRanges(loop_bedpe$V2+1, end = loop_bedpe$V3, names = c(1:nrow(loop_bedpe))),
    strand = Rle(strand('*')))
  Gr_anchor_2 = GRanges(
    seqnames = Rle(loop_bedpe$V4),
    ranges = IRanges(loop_bedpe$V5+1, end = loop_bedpe$V6, names = c(1:nrow(loop_bedpe))),
    strand = Rle(strand('*')))
  
  
  ##########  Converting from read pairs to bin size interactions
  overlaps = data.frame(findOverlaps(genome_window, Gr_anchor_1))
  overlaps = overlaps[order(overlaps$subjectHits),]
  temp = overlaps[!duplicated(overlaps$subjectHits),"queryHits"] # to remove cases where a read may overlap two genomic windows
  anchor_RNA = genome_window[temp]
  
  overlaps = data.frame(findOverlaps(genome_window, Gr_anchor_2))
  overlaps = overlaps[order(overlaps$subjectHits),]
  temp = overlaps[!duplicated(overlaps$subjectHits),"queryHits"] # to remove cases where a read may overlap two genomic windows
  anchor_DNA = genome_window[temp]
  
  # test
  # interaction_region_RNA = GRanges(paste0(chr_plot,":",3740001,"-",3750000))
  # overlaps_RNA = countOverlaps(anchor_RNA, interaction_region_RNA, ignore.strand = T)
  # anchor_RNA_overlaps = anchor_RNA[overlaps_RNA > 0]
  # 
  # anchor_DNA_overlaps = anchor_DNA[overlaps_RNA > 0]
  # 
  # temp = countOverlaps(anchor_DNA, GRanges(paste0(chr_plot,":",3260001,"-",3270000)))
  # temp1 = anchor_DNA[temp>0]
  # temp2=data.frame(anchor_RNA[names(temp1)])
  
  
  ########## Making interaction object
  interaction.object <- GenomicInteractions(
    anchor1 = anchor_RNA,
    anchor2 = anchor_DNA
  )
  
  df.interaction.object = data.frame(interaction.object)
  df.interaction.object_agg = ddply(df.interaction.object,.(seqnames1,start1,end1,seqnames2,start2,end2),nrow)
  
  anchor_RNA_agg = GRanges(
    seqnames = Rle(df.interaction.object_agg$seqnames1),
    ranges = IRanges(df.interaction.object_agg$start1, end = df.interaction.object_agg$end1, names = c(1:nrow(df.interaction.object_agg))),
    strand = Rle(strand('*')))
  
  anchor_DNA_agg = GRanges(
    seqnames = Rle(df.interaction.object_agg$seqnames2),
    ranges = IRanges(df.interaction.object_agg$start2, end = df.interaction.object_agg$end2, names = c(1:nrow(df.interaction.object_agg))),
    strand = Rle(strand('*')))
  
  interaction.object.agg <- GenomicInteractions(
    anchor1 = anchor_RNA_agg,
    anchor2 = anchor_DNA_agg,
    count = df.interaction.object_agg$V1
  )
  
  ##### Making plot
  ideogram.track <- IdeogramTrack(genome = 'hg38', chromosome = chr_plot)
  genome.axis.track <- GenomeAxisTrack()
  
  interaction.track <- InteractionTrack(
    interaction.object.agg,
    name = ''
  )
  
  out_file = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/tracks/",sample_label,"_",chr_plot,"_",start_plot,"_",end_plot,"_loops.pdf")

  pdf(out_file, width = 8, height = 2)
  plotTracks(
    list(ideogram.track, interaction.track),
    sizes = c(0.1,0.8),
    from = start_plot,
    to = end_plot,
    col="blue",
    col.anchors.line = "black",
    col.anchors.fill = "white",
    cex = 1.5,
    fontcolor = "black",
    background.title = "white"
  )
  dev.off()
  
}

make_loop_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_HFF/merged_loops.bedpe",
                             chr_plot = "chr2",
                             start_plot=149300000,
                             end_plot=151600000,
                             sample_label = "HFF")


########### PLAC-seq

make_plac_interaction_track <- function(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                        plac_file,
                                        chr_plot,
                                        start_plot,
                                        end_plot,
                                        read_length # read length PLAC-seq
){
  
  ########## Load genome data
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ### Tile genome and make genomic windows for selected chromosomes
  genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = 10000, cut.last.tile.in.chrom = T)

  ########## Loading input margi data for selecting interactions
  plac_data = data.frame(fread(plac_file))
  #plac_data = plac_data[which(plac_data$X.chrom1 == plac_data$chrom2),]
  plac_data = plac_data[which(plac_data$V1 == chr_plot),]
  
  Gr_plac_start = GRanges(
    seqnames = Rle(plac_data$V1),
    ranges = IRanges(plac_data$V2, end = plac_data$V2+read_length, names = c(1:nrow(plac_data))),
    strand = Rle(strand(plac_data$V5)))
  
  Gr_plac_end = GRanges(
    seqnames = Rle(plac_data$V3),
    ranges = IRanges(plac_data$V4, end = plac_data$V4+read_length, names = c(1:nrow(plac_data))),
    strand = Rle(strand(plac_data$V6)))
  
  ########## Selecting ranges corresponding to interaction region based also on pattern type
  interaction_region = GRanges(paste0(chr_plot,":",start_plot,"-",end_plot))
  
  overlaps_start = countOverlaps(Gr_plac_start, interaction_region, ignore.strand = T)
  Gr_plac_start_overlaps = Gr_plac_start[overlaps_start > 0]
  
  overlaps_end = countOverlaps(Gr_plac_end, interaction_region, ignore.strand = T)
  Gr_plac_end_overlaps = Gr_plac_end[overlaps_end > 0]
  
  temp = intersect(names(Gr_plac_start_overlaps), names(Gr_plac_end_overlaps))
  
  Gr_plac_start_overlaps = Gr_plac_start_overlaps[temp]
  Gr_plac_end_overlaps = Gr_plac_end_overlaps[temp]
  
  
  ##########  Converting from read pairs to bin size interactions
  overlaps = data.frame(findOverlaps(genome_window, Gr_plac_start_overlaps))
  overlaps = overlaps[order(overlaps$subjectHits),]
  temp = overlaps[!duplicated(overlaps$subjectHits),"queryHits"] # to remove cases where a read may overlap two genomic windows
  anchor_start = genome_window[temp]
  
  overlaps = data.frame(findOverlaps(genome_window, Gr_plac_end_overlaps))
  overlaps = overlaps[order(overlaps$subjectHits),]
  temp = overlaps[!duplicated(overlaps$subjectHits),"queryHits"] # to remove cases where a read may overlap two genomic windows
  anchor_end = genome_window[temp]
  
  
  ########## Making interaction object
  interaction.object <- GenomicInteractions(
    anchor1 = anchor_start,
    anchor2 = anchor_end
  )
  
  df.interaction.object = data.frame(interaction.object)
  df.interaction.object_agg = ddply(df.interaction.object,.(seqnames1,start1,end1,seqnames2,start2,end2),nrow)
  
  anchor_start_agg = GRanges(
    seqnames = Rle(df.interaction.object_agg$seqnames1),
    ranges = IRanges(df.interaction.object_agg$start1, end = df.interaction.object_agg$end1, names = c(1:nrow(df.interaction.object_agg))),
    strand = Rle(strand('*')))
  
  anchor_end_agg = GRanges(
    seqnames = Rle(df.interaction.object_agg$seqnames2),
    ranges = IRanges(df.interaction.object_agg$start2, end = df.interaction.object_agg$end2, names = c(1:nrow(df.interaction.object_agg))),
    strand = Rle(strand('*')))
  
  interaction.object.agg <- GenomicInteractions(
    anchor1 = anchor_start_agg,
    anchor2 = anchor_end_agg,
    count = df.interaction.object_agg$V1
  )
  
  ##### Making plot
  ideogram.track <- IdeogramTrack(genome = 'hg38', chromosome = chr_plot)
  genome.axis.track <- GenomeAxisTrack()
  
  interaction.track <- InteractionTrack(
    interaction.object.agg,
    name = ''
  )
  
  png(paste0("/dataOS/rcalandrelli/phase_separation/PLACseq/tracks/H1_",chr_plot,"_",start_plot,"_",end_plot,"_interactions.png"), width = 8, height = 2, units = "in", res = 500)
  plotTracks(
    list(ideogram.track, interaction.track),
    sizes = c(0.1,0.8),
    from = start_plot,
    to = end_plot,
    col="red",
    col.anchors.line = "black",
    col.anchors.fill = "white",
    col.interactions = "blue",
    cex = 1.5,
    fontcolor = "black",
    background.title = "white"
  )
  dev.off()
  
}


make_plac_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             plac_file="/dataOS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.txt",
                             chr_plot = "chr1",
                             start_plot=3250000,
                             end_plot=3750000,
                            read_length = 100)

make_plac_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                             plac_file="/dataOS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.txt",
                             chr_plot = "chr20",
                             start_plot=37150000,
                             end_plot=37650000,
                            read_length = 100)

make_plac_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                            plac_file="/dataOS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.txt",
                            chr_plot = "chr9",
                            start_plot=96300000,
                            end_plot=96800000,
                            read_length = 100)

make_plac_interaction_track(chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                            plac_file="/dataOS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.txt",
                            chr_plot = "chr10",
                            start_plot=69600000,
                            end_plot=70100000,
                            read_length = 100)




###################################################################
#### Interaction frequency curve
extract_contact_frequency_margi <- function(directory, 
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
    temp = data.frame(fread(paste0(directory, "/", as.character(resolution), "_filter1k/", i,"_",i,".txt")))
    contact_matrix = matrix(0, nrow = length(chr_row_genome_window), ncol = length(chr_col_genome_window))
    contact_matrix[as.matrix(temp[,c("RNA_end_bins","DNA_end_bins")])] = temp$x
    
    rownames(contact_matrix) = seq(1,nrow(contact_matrix)) # ???
    contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
    
    #contact_matrix = as.matrix(read.table(paste0(directory,i,"_",i,"_",as.character(resolution),".txt"), stringsAsFactors = F))
    contact_frequency = c()
    for (bin_distance in seq(start_bin_distance, end_bin_distance, step_bin_distance)){
      k = 1
      contact_frequency_bin_distance = 0
      while(k + bin_distance < nrow(contact_matrix)){
        contact_frequency_bin_distance = contact_frequency_bin_distance + contact_matrix[k,k+bin_distance] + contact_matrix[k+bin_distance,k]
        k = k + 1
      }
      contact_frequency = c(contact_frequency, contact_frequency_bin_distance / k)
    }
    contact_frequency_list[[i]] = contact_frequency
  }
  return(contact_frequency_list)
}

### Control
temp = extract_contact_frequency_margi(directory="/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_control/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=10)
contact_frequency_control_margi = colMeans(do.call("rbind", temp))

### NH4OAc
temp = extract_contact_frequency_margi(directory="/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_NH4OAc/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=10)
contact_frequency_NH4OAc_margi = colMeans(do.call("rbind", temp))

### FL
temp = extract_contact_frequency_margi(directory="/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_FL/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=10)
contact_frequency_FL_margi = colMeans(do.call("rbind", temp))

### RNase
temp = extract_contact_frequency_margi(directory="/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_RNase/",
                                 resolution=50000,
                                 max_distance=20000000,
                                 n_points=10)
contact_frequency_RNase_margi = colMeans(do.call("rbind", temp))




  