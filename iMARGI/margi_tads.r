####################### TAD analysis on iMARGI data
library(gtools)
library(GenomicRanges)
library(GenomicAlignments)
library(data.table)

options(scipen=999)

########################################## 
# Function to calculate the intensity of each genomic window at a defined resolution (100kb) based on the number of contacts over a pre-defined flanking
# region (+-1Mb) compared to the average intensity across all the genomic windows. Output: 1 if the bin is enriched in intensity, 0 if it's depleted.
# Dataframe with this information is saved as output
compute_genome_imargi_intensity <- function(contact_matrix_path,
                                            chromosomes = as.character(chrSize$V1),
                                            resolution=100000,
                                            flanking_region = 1000000,
                                            Gr_genome=Gr_hg38,
                                            output_path){
  
  shift_bin = flanking_region/resolution
  genome_window_chr_list = list()
  
  # Function to calculate the intensity of bin i
  compute_bin_intensity <- function(i){
    #print(i)
    upstream = ifelse(i - shift_bin > 0, i - shift_bin, 1) # most upstream bin of flanking region
    downstream = ifelse(i + shift_bin <= length(Gr_genome_window_chr), i + shift_bin, length(Gr_genome_window_chr)) # most downstream bin of flanking region
    
    # Calculate intensity for the RNA end bin (on the rows)
    contact_matrix_row_bin = contact_matrix_sparse[which(contact_matrix_sparse$RNA_end_bins == i),]
    contact_matrix_row_bin_flanking = contact_matrix_row_bin[which(contact_matrix_row_bin$DNA_end_bins %in% seq(upstream,downstream)),]
    row_intensity = ifelse(nrow(contact_matrix_row_bin_flanking) > 0, sum(contact_matrix_row_bin_flanking$x), 0)
    
    # Calculate intensity for the DNA end bin (on the columns)
    # contact_matrix_col_bin = contact_matrix_sparse[which(contact_matrix_sparse$DNA_end_bins == i),]
    # contact_matrix_col_bin_flanking = contact_matrix_col_bin[which(contact_matrix_col_bin$RNA_end_bins %in% seq(upstream,downstream)),]
    # col_intensity = ifelse(nrow(contact_matrix_col_bin_flanking) > 0, sum(contact_matrix_col_bin_flanking$x), 0)
    
    intensity_bin = row_intensity #+ col_intensity #- contact_matrix_sparse[which(contact_matrix_sparse$RNA_end_bins == i & contact_matrix_sparse$DNA_end_bins == i),"x"]
    return(intensity_bin)
  }
  
  for (chromosome in chromosomes){
    # Load contact matrix
    contact_matrix_sparse = data.frame(fread(paste0(contact_matrix_path,"/",chromosome,"_",chromosome,".txt"), header = T))
    # Make genomic windows for specific chromosome
    Gr_genome_window_chr = tileGenome(seqinfo(Gr_genome)[chromosome], tilewidth = resolution, cut.last.tile.in.chrom = T)
    # Calculate intensities per each bin
    mcols(Gr_genome_window_chr)["intensity"] = sapply(1:length(Gr_genome_window_chr), compute_bin_intensity)
    genome_window_chr_list[[chromosome]] = data.frame(Gr_genome_window_chr)
  }
  
  genome_window_chr = do.call("rbind", genome_window_chr_list)
  avg_intensity = mean(genome_window_chr$intensity)
  genome_window_chr$enrichment = ifelse(genome_window_chr$intensity>avg_intensity, 1, 0)
  write.table(genome_window_chr, paste0(output_path,"/genome_RNA_window_intensity_", resolution, ".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  
}

########################################## 
# Function to see if a TAD is: 1) enriched of high intensity bins (more than 80% of the bins are enriched in intensity) 2) depleted of high intensity bins
# (more than 80% of the bins are depleted in intensity) 3) neither in the rest of the cases. Output is a data frame with all the TADs and the
# corresponding label
library(angsdr)

label_tad_enrichment_depletion <- function(intensity_file_path,
                                           chromosomes = as.character(chrSize$V1),
                                           resolution=100000,
                                           tad_file_path,
                                           tad_resolution=10000,
                                           thres=0.8,
                                           output_path){
  
  Gr_window_intensity = makeGRangesFromDataFrame(data.frame(fread(paste0(intensity_file_path,"/genome_RNA_window_intensity_", resolution, ".txt"), header = T)),
                           keep.extra.columns=T,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("seqnames", "seqname",
                                            "chromosome", "chrom",
                                            "chr", "chromosome_name",
                                            "seqid"),
                           start.field="start",
                           end.field=c("end", "stop"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)
  
  tads_file_list = list()
  for (chromosome in chromosomes){
    if (file.exists(paste0(tad_file_path,"/TADs_",chromosome,"/", tad_resolution, "_blocks.bedpe"))){
      tads_file = read.table(paste0(tad_file_path,"/TADs_",chromosome,"/", tad_resolution, "_blocks.bedpe"))[,c(1:3)]
      tads_file[,1] = paste0("chr", tads_file[,1])
      tads_file_list[[chromosome]] = tads_file
    }
  }
  
  tads_all = do.call("rbind", tads_file_list)
  
  Gr_tads <- GRanges(
    seqnames = Rle(tads_all$V1),
    ranges = IRanges(as.numeric(tads_all$V2), end = as.numeric(tads_all$V3), names = c(1:nrow(tads_all))),
    strand = Rle(strand('*')))
  Gr_tads <- sort(Gr_tads)
  
  overlaps = findOverlaps(Gr_tads, Gr_window_intensity)
  
  tads_df = data.frame(Gr_tads)
  tads_df$index = seq(1,nrow(tads_df))
  
  # Weizhong's method: 80% of the overlapping tads have to be enriched
  # compute_tad_label_1 <- function(x){
  #   tad_overlap = overlaps[which(queryHits(overlaps)==x)] # x["index"]
  #   bin_tad_overlap_enrich = mcols(Gr_window_intensity[subjectHits(tad_overlap)])[,"enrichment"]
  #   
  #   tad_label = ifelse(sum(bin_tad_overlap_enrich==1) >= round(length(bin_tad_overlap_enrich)*thres),
  #                      "enrichment",
  #                      ifelse(sum(bin_tad_overlap_enrich==0) >= round(length(bin_tad_overlap_enrich)*thres),
  #                             "depletion",
  #                             "neither"))
  #   return(tad_label)
  # }
  # 
  # tads_df$label = sapply(tads_df$index, compute_tad_label)
  # write.table(tads_df, paste0(output_path,"/", tad_resolution, "_tads_RNA_enrich_label.txt"), row.names = F, col.names = T, sep="\t", quote = F)
  
  # 80% of the TAD bins overlapping enriched tads
  compute_tad_label_2 <- function(x){
    tad_overlap = overlaps[which(queryHits(overlaps)==x)] # x["index"]
    
    tad_bin_intensity = c()
    for (i in subjectHits(tad_overlap)){
      gr.over <- pintersect(Gr_tads[x], Gr_window_intensity[i])
      gr.counts <- tapply(gr.over, x, FUN=function(x) sum(width(x)))
      tad_bin_intensity = c(tad_bin_intensity, rep(mcols(Gr_window_intensity[i])["enrichment"][1,1], as.integer(gr.counts/10000)))
    }
    
    perc_thres = as.integer(length(tad_bin_intensity) * thres)
    
    tad_label = ifelse(sum(tad_bin_intensity==1) >= perc_thres, "enrichment",
                       ifelse(sum(tad_bin_intensity==0) >= perc_thres, "depletion", "neither"))
    return(tad_label)
  }
  
  tads_df$label = sapply(tads_df$index, compute_tad_label_2)
  write.table(tads_df, paste0(output_path,"/", tad_resolution, "_tads_RNA_enrich_label_2.txt"), row.names = F, col.names = T, sep="\t", quote = F)
}


########################################## 
# Function which permutes the intensity labels across bins and then re-labels each TADs as above. The output is a dataframe with all the TADs and then
# n_permutation additional columns with the labels associated to the TADs after each permutation.
label_tad_enrichment_depletion_permutation <- function(intensity_file_path,
                                                      chromosomes,
                                                      resolution=100000,
                                                      tad_file_path,
                                                      tad_resolution=10000,
                                                      thres=0.8,
                                                      output_path,
                                                      n_permutations = 100){
  
  tads_file_list = list()
  for (chromosome in chromosomes){
    if (file.exists(paste0(tad_file_path,"/TADs_",chromosome,"/", tad_resolution, "_blocks.bedpe"))){
      tads_file = read.table(paste0(tad_file_path,"/TADs_",chromosome,"/", tad_resolution, "_blocks.bedpe"))[,c(1:3)]
      tads_file[,1] = paste0("chr", tads_file[,1])
      tads_file_list[[chromosome]] = tads_file
    }
  }
  
  tads_all = do.call("rbind", tads_file_list)
  
  Gr_tads <- GRanges(
    seqnames = Rle(tads_all$V1),
    ranges = IRanges(as.numeric(tads_all$V2), end = as.numeric(tads_all$V3), names = c(1:nrow(tads_all))),
    strand = Rle(strand('*')))
  Gr_tads = sort(Gr_tads)
  
  tads_df = data.frame(Gr_tads)
  tads_df$index = seq(1,nrow(tads_df))
  
  Gr_window_intensity = makeGRangesFromDataFrame(data.frame(fread(paste0(intensity_file_path,"/genome_RNA_window_intensity_", resolution, ".txt"), header = T)),
                                                 keep.extra.columns=T,
                                                 ignore.strand=FALSE,
                                                 seqinfo=NULL,
                                                 seqnames.field=c("seqnames", "seqname",
                                                                  "chromosome", "chrom",
                                                                  "chr", "chromosome_name",
                                                                  "seqid"),
                                                 start.field="start",
                                                 end.field=c("end", "stop"),
                                                 strand.field="strand",
                                                 starts.in.df.are.0based=FALSE)
  
  label_list = list()
  for (i in 1:n_permutations){
    print(paste0("permutation ", i))

    mcols(Gr_window_intensity)["enrichment_permuted"] = permute(mcols(Gr_window_intensity)[,"enrichment"]) # permute enrichment values
    
    overlaps = findOverlaps(Gr_tads, Gr_window_intensity)
    
    # compute_tad_label <- function(x){
    #   tad_overlap = overlaps[which(queryHits(overlaps)==x)]
    #   bin_tad_overlap_enrich = mcols(Gr_window_intensity[subjectHits(tad_overlap)])[,"enrichment_permuted"]
    #   
    #   tad_label = ifelse(sum(bin_tad_overlap_enrich==1) >= round(length(bin_tad_overlap_enrich)*thres),
    #                      "enrichment",
    #                      ifelse(sum(bin_tad_overlap_enrich==0) >= round(length(bin_tad_overlap_enrich)*thres),
    #                             "depletion",
    #                             "neither"))
    #   return(tad_label)
    # }
    
    # 80% of the TAD bins overlapping enriched tads
    compute_tad_label_2 <- function(x){
      tad_overlap = overlaps[which(queryHits(overlaps)==x)] # x["index"]
      
      tad_bin_intensity = c()
      for (i in subjectHits(tad_overlap)){
        gr.over <- pintersect(Gr_tads[x], Gr_window_intensity[i])
        gr.counts <- tapply(gr.over, x, FUN=function(x) sum(width(x)))
        tad_bin_intensity = c(tad_bin_intensity, rep(mcols(Gr_window_intensity[i])["enrichment"][1,1], as.integer(gr.counts/10000)))
      }
      
      perc_thres = as.integer(length(tad_bin_intensity) * thres)
      
      tad_label = ifelse(sum(tad_bin_intensity==1) >= perc_thres, "enrichment",
                         ifelse(sum(tad_bin_intensity==0) >= perc_thres, "depletion", "neither"))
      return(tad_label)
    }
    
    label_list[[i]] = sapply(tads_df$index, compute_tad_label_2)
  }
  
  all_permutation_labels = do.call("cbind", label_list)
  tads_df_permutation = cbind(tads_df, all_permutation_labels)
  colnames(tads_df_permutation)[7:ncol(tads_df_permutation)] = paste0("perm_",1:n_permutations)
  write.table(tads_df_permutation, paste0(output_path,"/", tad_resolution, "_tads_RNA_enrich_label_permutations_",n_permutations,"_2.txt"), row.names = F, col.names = T, sep="\t", quote = F)

}


########################################## Run the functions
chrSize_file = "/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes"
chrSize = read.table(chrSize_file)

Gr_hg38 <- GRanges(
  seqnames = Rle(as.character(chrSize$V1)),
  ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
  strand = Rle(strand('*')))
seqlengths(Gr_hg38) <- as.numeric(chrSize$V2)

all_imargi_samples = c("iMARGI_H1_control", "iMARGI_H1_FL", "iMARGI_H1_NH4OAc", "iMARGI_H1_RNase")
all_imargi_bio_replicates = c(paste0("iMARGI_H1_control_", c(1,3,4,5)),
                              paste0("iMARGI_H1_FL_", c(1,2)),
                              paste0("iMARGI_H1_NH4OAc_", c(1,2)),
                              paste0("iMARGI_H1_RNase_", c(1,2)))

for (sample in all_imargi_samples){
  sample_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample)
  if (!dir.exists(sample_dir)){
    dir.create(sample_dir)
  }
}

for (sample in all_imargi_samples){
  compute_genome_imargi_intensity(contact_matrix_path=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/",sample,"/10000_filter1k"),
                                  chromosomes = as.character(chrSize$V1),
                                  resolution=100000,
                                  flanking_region = 1000000,
                                  Gr_genome=Gr_hg38,
                                  output_path=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample))
}


for (sample in all_imargi_samples){
  temp = strsplit(sample,"_")[[1]][3]
  tad_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/tads/H1_",temp,"_merged")

  label_tad_enrichment_depletion(intensity_file_path=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample),
                                 chromosomes = as.character(chrSize$V1),
                                 resolution=100000,
                                 tad_file_path=tad_file_path_temp,
                                 tad_resolution=10000,
                                 thres=0.8,
                                 output_path=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample))
}


for (sample in all_imargi_samples){
  temp = strsplit(sample,"_")[[1]][3]
  tad_file_path_temp = paste0("/dataOS/rcalandrelli/phase_separation/HiC/tads/H1_",temp,"_merged")

  label_tad_enrichment_depletion_permutation(intensity_file_path=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample),
                                 chromosomes = as.character(chrSize$V1),
                                 resolution=100000,
                                 tad_file_path=tad_file_path_temp,
                                 thres=0.8,
                                 output_path=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample),
                                 n_permutations = 100)
}



######################## Summarizing TAD enrichment results
summary_TAD_label_list = list()
for (sample in all_imargi_samples){
  temp_data = read.table(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample, "/10000_tads_RNA_enrich_label_2.txt"), header = T)
  summary_TAD_label_list[[sample]] = table(as.factor(temp_data$label))
}

summary_TAD_label_bio_rep_list = list()
for (sample in all_imargi_bio_replicates[2:11]){
  temp_data = read.table(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample, "/10000_tads_RNA_enrich_label_2.txt"), header = T)
  summary_TAD_label_bio_rep_list[[sample]] = table(as.factor(temp_data$label))
}

### Plot by merged sample
summary_TAD_label = data.frame(do.call("rbind", summary_TAD_label_list))
summary_TAD_label = summary_TAD_label / rowSums(summary_TAD_label)
summary_TAD_label$sample = gsub("iMARGI_","",rownames(summary_TAD_label))
df = melt(summary_TAD_label)
png("/dataOS/rcalandrelli/phase_separation/MARGI/tads/TAD_label_barplot.png", width = 6, height = 5, res = 200, units = "in")
p1<-ggplot(df, aes(x=sample, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('blue','green',"#ff00ff")) +
  xlab("") +
  ylab("Percentage of TADs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_blank())
dev.off()

### Plot by biological replciates with error bars
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
df$sample_type = gsub("H1_","",df$sample_type)
df$sample_type = gsub("control","Control",df$sample_type)
df$sample_type = factor(df$sample_type, levels = c("Control","NH4OAc","FL","RNase"))

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

png("/dataOS/rcalandrelli/phase_separation/MARGI/tads/TAD_label_barplot_merged_bio_rep.png", width = 12, height = 5, res = 200, units = "in")
plot_grid(p1,p2, labels = c("A", "B"))
dev.off()


### Save data table

# Merged data
summary_TAD_label = data.frame(do.call("rbind", summary_TAD_label_list))
summary_TAD_label$depletion_perc = summary_TAD_label$depletion / rowSums(summary_TAD_label)
summary_TAD_label$enrichment_perc = summary_TAD_label$enrichment / rowSums(summary_TAD_label)
summary_TAD_label$neither_perc = summary_TAD_label$neither / rowSums(summary_TAD_label)
write.table(summary_TAD_label, "/dataOS/rcalandrelli/phase_separation/MARGI/tads/summary_TAD_label.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# Biological replicates
summary_TAD_label_bio_rep = data.frame(do.call("rbind", summary_TAD_label_bio_rep_list))

summary_TAD_label_bio_rep$depletion_perc = summary_TAD_label_bio_rep$depletion / rowSums(summary_TAD_label_bio_rep)
summary_TAD_label_bio_rep$enrichment_perc = summary_TAD_label_bio_rep$enrichment / rowSums(summary_TAD_label_bio_rep)
summary_TAD_label_bio_rep$neither_perc = summary_TAD_label_bio_rep$neither / rowSums(summary_TAD_label_bio_rep)
write.table(summary_TAD_label_bio_rep, "/dataOS/rcalandrelli/phase_separation/MARGI/tads/summary_TAD_label_bio_replicates.txt", row.names = T, col.names = T, sep = "\t", quote = F)


######################## Summarizing permutation results
summary_tads_permutation <- function(tads_df_permutation_file,
                                     tads_df_file,
                                     n_permutations,
                                     sample){
  ### Real data
  tads_df = data.frame(fread(tads_df_file))
  real_enrich = as.numeric(table(as.factor(tads_df$label))["enrichment"])
  real_deplet = as.numeric(table(as.factor(tads_df$label))["depletion"])
  real_neither = as.numeric(table(as.factor(tads_df$label))["neither"])

  ### Permutation data
  tads_df_permutation = data.frame(fread(tads_df_permutation_file))[,paste0("perm_",seq(1,n_permutations))]
  summary_permutation_list = list()
  for (i in 1:ncol(tads_df_permutation)){
    summary_permutation_list[[i]] = table(as.factor(tads_df_permutation[,i]))
  }
  summary_permutation = as.data.frame(do.call("rbind", summary_permutation_list))
  
  ### Plot 
  label_cex = 10
  p1<-ggplot(summary_permutation, aes(x=enrichment)) +
    geom_histogram(binwidth = 10) +
    xlab("Number of TADs") +
    ylab("Number of\npermutations") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_cex),
          axis.text.x = element_text(size = label_cex, color = "black"),
          axis.text.y = element_text(size = label_cex, color = "black"),
          axis.ticks.length = unit(0.2, "cm")) +
    geom_vline(xintercept = real_enrich, color = "green", size=1.5) +
      ggtitle(sample)
    
    p2<-ggplot(summary_permutation, aes(x=depletion)) +
      geom_histogram(binwidth = 10) +
      xlab("Number of TADs") +
      ylab("Number of\npermutations") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_cex),
          axis.text.x = element_text(size = label_cex, color = "black"),
          axis.text.y = element_text(size = label_cex, color = "black"),
          axis.ticks.length = unit(0.2, "cm")) +
      geom_vline(xintercept = real_deplet, color = "blue", size=1.5) +
      ggtitle(sample)
    
    p3<-ggplot(summary_permutation, aes(x=neither)) +
      geom_histogram(binwidth = 10) +
      xlab("Number of TADs") +
      ylab("Number of\npermutations") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_cex),
          axis.text.x = element_text(size = label_cex, color = "black"),
          axis.text.y = element_text(size = label_cex, color = "black"),
          axis.ticks.length = unit(0.2, "cm")) +
      geom_vline(xintercept = real_neither, color = "#ff00ff", size=1.5) +
      ggtitle(sample)
    
  p <- plot_grid(p1,p2,p3, nrow=1)
  return(p)
  
}

permutation_plot_list = list()
for (sample in all_imargi_samples[1:4]){
  permutation_plot_list[[sample]] = summary_tads_permutation(tads_df_permutation_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample, "/10000_tads_RNA_enrich_label_permutations_100.txt"),
                           tads_df_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/tads/", sample, "/10000_tads_RNA_enrich_label.txt"),
                           n_permutations=100,
                           sample=sample)
}

png("/dataOS/rcalandrelli/phase_separation/MARGI/tads/permutation_barplots.png", width = 8, height = 8, res = 200, units = "in")
plot_grid(plotlist = permutation_plot_list, nrow = 4)
dev.off()







