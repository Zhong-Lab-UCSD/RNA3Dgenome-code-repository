########## TAD and TAD boundaries analysis for RNA-genome interaction project
library(GenomicRanges)
library(GenomicAlignments)
library(reshape2)

chrSize = read.table("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes")

hg38_chromosomes = as.character(chrSize$V1)
hg38_lengths = as.numeric(chrSize$V2)
names(hg38_lengths) = hg38_chromosomes

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/"

##### TAD analysis global assessment through Measure of Concordance

create_Gr_TADs <- function(TAD_directory, # with trailing slash in the end
                           sample,
                           resolution){
  tads = matrix(nrow=0, ncol=9)
  for (i in hg38_chromosomes){
    if (file.exists(paste0(TAD_directory,sample,"/TADs_",i,"/",resolution,"_blocks.bedpe"))){
      temp_file = read.table(paste0(TAD_directory,sample,"/TADs_",i,"/",resolution,"_blocks.bedpe"), stringsAsFactors = F)
      temp_file = temp_file[,c(1:3,11:16)]
      temp_file[,1] = paste0("chr",temp_file[,1])
      tads = rbind(tads,temp_file)
    }
  }
  colnames(tads) = c("chr","start","end","color","corner_score","Uvar","Lvar","Usign","Lsign")
  Gr_tads <- GRanges(
    seqnames = Rle(tads[,1]),
    ranges = IRanges(as.numeric(tads[,2]), end = as.numeric(tads[,3]), names = c(1:nrow(tads))),
    strand = Rle(strand("*")),
    color = as.character(tads[,4]),
    corner_score = as.numeric(tads[,5]),
    Uvar = as.numeric(tads[,6]),
    Lvar = as.numeric(tads[,7]),
    Usign = as.numeric(tads[,8]),
    Lsign = as.numeric(tads[,9]))
  return(Gr_tads)
}

# 50k resolution
# Gr_tads_control_50000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/",
#                                           "H1_control_merged",
#                                           50000)
# Gr_tads_NH4OAc_50000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/",
#                                          "H1_NH4OAc_merged",
#                                          50000)
# Gr_tads_FL_50000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/",
#                                          "H1_FL_merged",
#                                          50000)
# Gr_tads_RNaseTreat_50000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/",
#                                          "H1_RNaseTreat_merged",
#                                          50000)

# 10k resolution
Gr_tads_control_10000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/",
                                       "H1_control_merged",
                                       10000)
Gr_tads_NH4OAc_10000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/",
                                      "H1_NH4OAc_merged",
                                      10000)
Gr_tads_FL_10000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/",
                                  "H1_FL_merged",
                                  10000)
Gr_tads_RNase_10000 = create_Gr_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/",
                                          "H1_RNase_merged",
                                          10000)


# sample = rep(c("Ctrl","NH4OAc","FL","RNase"),2)
# resolution = c(rep("50k",4), rep("10k",4))
# n_tads = c(length(Gr_tads_control_50000), length(Gr_tads_NH4OAc_50000), length(Gr_tads_FL_50000), length(Gr_tads_RNaseTreat_50000),
#            length(Gr_tads_control_10000), length(Gr_tads_NH4OAc_10000), length(Gr_tads_FL_10000), length(Gr_tads_RNase_10000))
# df = data.frame(sample, resolution, n_tads)
# df$sample = factor(df$sample, levels = c("Ctrl","NH4OAc","FL","RNase"))
# 
# png(paste0(directory,"result/barplot_n_tads.png"), width = 8, height = 8, units = "in", res = 200)
# ggplot(df, aes(fill=resolution, y=n_tads, x=sample)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   #scale_fill_manual(values = c("#ffeb00","blue")) +
#   labs(x = "", y = "Number of TADs") + 
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=24),
#         axis.text.x = element_text(size = 24, color = "black"),
#         axis.text.y = element_text(size = 24, color = "black"),
#         axis.ticks.length = unit(0.2, "cm"))
# dev.off()

sample = c("Control","NH4OAc","FL","RNase")
n_tads = c(length(Gr_tads_control_10000), length(Gr_tads_NH4OAc_10000), length(Gr_tads_FL_10000), length(Gr_tads_RNase_10000))
df = data.frame(sample, n_tads)
df$sample = factor(df$sample, levels = c("Control","NH4OAc","FL","RNase"))

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/barplot_n_tads.png", width = 8, height = 8, units = "in", res = 200)
ggplot(df, aes(y=n_tads, x=sample)) +
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "", y = "Number of TADs") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()

# Plot TAD length distribution
sample = c(rep("Control",length(Gr_tads_control_10000)),
           rep("FL",length(Gr_tads_FL_10000)),
           rep("NH4OAc",length(Gr_tads_NH4OAc_10000)),
           rep("RNase",length(Gr_tads_RNase_10000)))
length_tads = c(end(ranges(Gr_tads_control_10000)) - start(ranges(Gr_tads_control_10000)),
                end(ranges(Gr_tads_FL_10000)) - start(ranges(Gr_tads_FL_10000)), 
                end(ranges(Gr_tads_NH4OAc_10000)) - start(ranges(Gr_tads_NH4OAc_10000)), 
                end(ranges(Gr_tads_RNase_10000)) - start(ranges(Gr_tads_RNase_10000)))
df = data.frame(sample, length_tads)

png("/dataOS/rcalandrelli/phase_separation/HiC/tads/result/TAD_size_boxplot_log.png", width = 8, height = 8, units = "in", res = 200)
ggplot(df, aes(y=log10(length_tads), x=sample)) +
  geom_boxplot() +
  labs(x = "", y = "Log10(TAD size)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()

wilcox.test(end(ranges(Gr_tads_control_10000)) - start(ranges(Gr_tads_control_10000)), end(ranges(Gr_tads_NH4OAc_10000)) - start(ranges(Gr_tads_NH4OAc_10000)))
wilcox.test(end(ranges(Gr_tads_control_10000)) - start(ranges(Gr_tads_control_10000)), end(ranges(Gr_tads_FL_10000)) - start(ranges(Gr_tads_FL_10000)))
wilcox.test(end(ranges(Gr_tads_control_10000)) - start(ranges(Gr_tads_control_10000)), end(ranges(Gr_tads_RNase_10000)) - start(ranges(Gr_tads_RNase_10000)))


# Plot corner score distribution
df1 = data.frame(Gr_tads_control_10000)[,c("seqnames","corner_score")]
df1[,1] = "Control"

df2 = data.frame(Gr_tads_NH4OAc_10000)[,c("seqnames","corner_score")]
df2[,1] = "NH4OAc"

df3 = data.frame(Gr_tads_FL_10000)[,c("seqnames","corner_score")]
df3[,1] = "FL"

df4 = data.frame(Gr_tads_RNase_10000)[,c("seqnames","corner_score")]
df4[,1] = "RNase"

df = rbind(df1,df2,df3,df4)
colnames(df)[1] = "sample"
df$sample = factor(df$sample, levels = c("Control", "NH4OAc", "FL", "RNase"))

png("/dataOS/rcalandrelli/phase_separation/TADs/result/corner_score_boxplot.png", width = 8, height = 8, units = "in", res = 200)
ggplot(df, aes(x=sample, y=corner_score)) + 
  geom_boxplot() +
  labs(x="", y="TAD corner score") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()

wilcox.test(df1$corner_score,
            df2$corner_score,
            paired = F)

wilcox.test(df1$corner_score,
            df3$corner_score,
            paired = F)

wilcox.test(df1$corner_score,
            df4$corner_score,
            paired = F)


##### Measure of Concordance analysis

calculate_MoC <- function(Gr1,Gr2){
  hits <- findOverlaps(Gr1,Gr2)
  overlaps <- pintersect(Gr1[queryHits(hits)], Gr2[subjectHits(hits)])
  F_ij <- as.numeric(width(overlaps))
  P_i <- as.numeric(width(Gr1[queryHits(hits)]))
  Q_j <- as.numeric(width(Gr2[subjectHits(hits)]))
  
  MoC <- 1/(sqrt(length(Gr1) * length(Gr2)) - 1) * (sum(F_ij^2 / (P_i * Q_j)) - 1)
  return(MoC)
}

# MoC_control_NH4OAc_50000 <- calculate_MoC(Gr_tads_control_50000,Gr_tads_NH4OAc_50000)
# MoC_control_FL_50000 <- calculate_MoC(Gr_tads_control_50000,Gr_tads_FL_50000)
# MoC_control_RNaseTreat_50000 <- calculate_MoC(Gr_tads_control_50000,Gr_tads_RNaseTreat_50000)
# MoC_NH4OAc_FL_50000 <- calculate_MoC(Gr_tads_NH4OAc_50000,Gr_tads_FL_50000)
# MoC_NH4OAc_RNaseTreat_50000 <- calculate_MoC(Gr_tads_NH4OAc_50000,Gr_tads_RNaseTreat_50000)
# MoC_FL_RNaseTreat_50000 <- calculate_MoC(Gr_tads_FL_50000,Gr_tads_RNaseTreat_50000)

MoC_control_NH4OAc_10000 <- calculate_MoC(Gr_tads_control_10000,Gr_tads_NH4OAc_10000)
MoC_control_FL_10000 <- calculate_MoC(Gr_tads_control_10000,Gr_tads_FL_10000)
MoC_control_RNase_10000 <- calculate_MoC(Gr_tads_control_10000,Gr_tads_RNase_10000)
MoC_NH4OAc_FL_10000 <- calculate_MoC(Gr_tads_NH4OAc_10000,Gr_tads_FL_10000)
MoC_NH4OAc_RNase_10000 <- calculate_MoC(Gr_tads_NH4OAc_10000,Gr_tads_RNase_10000)
MoC_FL_RNase_10000 <- calculate_MoC(Gr_tads_FL_10000,Gr_tads_RNase_10000)

##### Mean interaction frequencies over TAD boundaries
calculate_mean_interaction_frequencies_TAD <- function(directory, # with trailing slash in the end
                                                       sample,
                                                       bin_size,
                                                       template_dim){
  output_matrix = matrix(0, nrow=template_dim/bin_size, ncol=template_dim/bin_size)
  counter = 0
  for (i in hg38_chromosomes){
    # Load sparse matrix (check if file size is non-zero first!!!)
    if (file.size(paste0(directory, "HiC_contact_matrices/KR/", sample, "/", bin_size, "/",i,"_",i,"_",bin_size,".txt")) != 0){
      print(paste0("Calculating interaction frequencies for TADs in ",i))
      input_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/", sample, "/", bin_size, "/",i,"_",i,"_",bin_size,".txt"), stringsAsFactors = F)
      input_matrix_sparse[which(is.na(input_matrix_sparse[,3])),3] = 0 # removing NaN values
      input_matrix_sparse_bin = input_matrix_sparse
      input_matrix_sparse_bin[,1] = input_matrix_sparse_bin[,1] / bin_size + 1
      input_matrix_sparse_bin[,2] = input_matrix_sparse_bin[,2] / bin_size + 1
      # Load TADs (check if file exists first!!!)
      if (file.exists(paste0(directory, "tads/", sample, "/TADs_",i,"/10000_blocks.bedpe"))){
        input_tads = read.table(paste0(directory, "tads/", sample, "/TADs_",i,"/10000_blocks.bedpe"), stringsAsFactors = F)
        
        tad_boundaries = unique(c(input_tads[,2], input_tads[,3])) # to consider each boundary only once
        tad_boundaries_bin = tad_boundaries / bin_size + 1
        
        out_tad_boundaries = matrix(nrow = length(tad_boundaries), ncol = 3)
        out_tad_boundaries[,1] = i
        out_tad_boundaries[,2] = tad_boundaries
        
        for (j in 1:length(tad_boundaries_bin)){
          boundary_bin = tad_boundaries_bin[j]
          shift = template_dim/2
          shift_bin = shift / bin_size
          if (boundary_bin - shift_bin > 0 & boundary_bin + shift_bin - 1 <= hg38_lengths[i] %/% bin_size){
            counter = counter + 1
            
            upstream_bin = as.numeric(boundary_bin - shift_bin) 
            downstream_bin = as.numeric(boundary_bin + shift_bin - 1) 
            input_matrix_sparse_bin_template = input_matrix_sparse_bin[which(input_matrix_sparse_bin[,1] >= upstream_bin &
                                                                               input_matrix_sparse_bin[,1] <= downstream_bin &
                                                                               input_matrix_sparse_bin[,2] >= upstream_bin &
                                                                               input_matrix_sparse_bin[,2] <= downstream_bin),]
            
            input_matrix_sparse_bin_template_scaled <- input_matrix_sparse_bin_template
            input_matrix_sparse_bin_template_scaled[,1] = input_matrix_sparse_bin_template_scaled[,1] - upstream_bin + 1
            input_matrix_sparse_bin_template_scaled[,2] = input_matrix_sparse_bin_template_scaled[,2] - upstream_bin + 1
            input_matrix_sparse_bin_template_scaled = as.matrix(input_matrix_sparse_bin_template_scaled)
            
            output_matrix[input_matrix_sparse_bin_template_scaled[,1:2]] = output_matrix[input_matrix_sparse_bin_template_scaled[,1:2]] + input_matrix_sparse_bin_template_scaled[,3]
            
            # Calculate TAD boundary insulation score
            mean_L = sum(input_matrix_sparse_bin_template_scaled[which(input_matrix_sparse_bin_template_scaled[,1] <= shift_bin &
                                                                         input_matrix_sparse_bin_template_scaled[,2] <= shift_bin),3]) / (shift_bin * (shift_bin + 1) / 2)
            mean_R = sum(input_matrix_sparse_bin_template_scaled[which(input_matrix_sparse_bin_template_scaled[,1] > shift_bin &
                                                                         input_matrix_sparse_bin_template_scaled[,2] > shift_bin),3]) / (shift_bin * (shift_bin + 1) / 2)
            mean_X = sum(input_matrix_sparse_bin_template_scaled[which(input_matrix_sparse_bin_template_scaled[,1] <= shift_bin &
                                                                         input_matrix_sparse_bin_template_scaled[,2] > shift_bin),3]) / (shift_bin^2)
            insulation_score = max(mean_L, mean_R) / mean_X
            out_tad_boundaries[j,3] = insulation_score
            
            # if (j <= length(tad_boundaries) / 2){
            #   input_tads[j, 17] = insulation_score
            # } else {
            #   input_tads[j - length(tad_boundaries) / 2, 18] = insulation_score
            # }
            
          }
        }
        write.table(out_tad_boundaries, paste0(directory, "tads/", sample, "/TADs_",i,"/10000_boundary_insulation_score.txt"), row.names = F, col.names = F, sep = '\t', quote = F)
      }
    }
  }
  output_matrix_average = output_matrix / counter
  write.table(output_matrix_average, paste0(directory, "tads/", sample,"/mean_interaction_frequencies_TAD_heatmap_", bin_size, ".txt"), sep="\t", row.names = F, col.names = F, quote = F)
}

# output_matrix = matrix(0,4,4)
# index_matrix = matrix(c(2,2,4,2,4,3),2,3)
# output_matrix[index_matrix[,1:2]-1] <- output_matrix[index_matrix[,1:2]-1] + index_matrix[,3]


calculate_mean_interaction_frequencies_TAD("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                           "H1_control_merged",
                                           10000,
                                           2000000)
calculate_mean_interaction_frequencies_TAD("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                           "H1_NH4OAc_merged",
                                           10000,
                                           2000000)
calculate_mean_interaction_frequencies_TAD("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                           "H1_FL_merged",
                                           10000,
                                           2000000)
calculate_mean_interaction_frequencies_TAD("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                           "H1_RNase_merged",
                                           10000,
                                           2000000)

### Plotting MIF matrix
parsing_MIF_matrix <- function(directory,
                               sample,
                               bin_size,
                               template_dim){
  temp = read.table(paste0(directory, "tads/", sample,"/mean_interaction_frequencies_TAD_heatmap_", bin_size, ".txt"), stringsAsFactors = F)
  for (i in 1:nrow(temp)){
    for (j in 1:ncol(temp)){
      if (i > j){
        temp[i,j] = NA
      }
    }
  }
  
  mean_interaction_frequencies_TAD = log(as.matrix(temp) + 1)
  #mean_interaction_frequencies_TAD = as.matrix(temp)
  
  k = template_dim/bin_size
  x <- paste0("V",seq(1:k))
  y <- as.character(seq(1:k))
  data <- expand.grid(X=x, Y=y)
  temp_z = c()
  for (i in 1:nrow(mean_interaction_frequencies_TAD)){
    temp_z = c(temp_z,mean_interaction_frequencies_TAD[i,])
  }
  data$Z <- temp_z
  return(data)
}

calculate_insulation_score <- function(directory,
                                       sample,
                                       bin_size,
                                       template_dim){
  
  sum_L = 0
  sum_R = 0
  sum_X = 0
  shift_bin = template_dim/bin_size/2
  
  temp = read.table(paste0(directory, "tads/", sample,"/mean_interaction_frequencies_TAD_heatmap_", bin_size, ".txt"), stringsAsFactors = F)
  for (i in 1:nrow(temp)){
    for (j in 1:ncol(temp)){
      if (i <= shift_bin & j <= shift_bin & i <= j){
        sum_L = sum_L + temp[i,j]
      } else if (i > shift_bin & j > shift_bin & i <= j){
        sum_R = sum_R + temp[i,j]
      } else if (i <= shift_bin & j > shift_bin){
        sum_X = sum_X + temp[i,j]
      }
    }
  }
  mean_L = sum_L / (shift_bin * (shift_bin+1) / 2)
  mean_R = sum_R / (shift_bin * (shift_bin+1) / 2)
  mean_X = sum_X / (shift_bin^2)
  insulation_score = max(mean_L, mean_R) / mean_X
  return(insulation_score)
}


MIF_tad_control_heatmap = parsing_MIF_matrix("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                     "H1_control_merged",
                                     10000,
                                     2000000)
max(MIF_tad_control_heatmap$Z[!is.na(MIF_tad_control_heatmap$Z)])
min(MIF_tad_control_heatmap$Z[!is.na(MIF_tad_control_heatmap$Z)])
calculate_insulation_score("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                           "H1_control_merged",
                           10000,
                           2000000)

MIF_tad_NH4OAc_heatmap = parsing_MIF_matrix("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                     "H1_NH4OAc_merged",
                                     10000,
                                     2000000)
max(MIF_tad_NH4OAc_heatmap$Z[!is.na(MIF_tad_NH4OAc_heatmap$Z)])
min(MIF_tad_NH4OAc_heatmap$Z[!is.na(MIF_tad_NH4OAc_heatmap$Z)])
calculate_insulation_score("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                            "H1_NH4OAc_merged",
                                            10000,
                                            2000000)

MIF_tad_FL_heatmap = parsing_MIF_matrix("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                     "H1_FL_merged",
                                     10000,
                                     2000000)
max(MIF_tad_FL_heatmap$Z[!is.na(MIF_tad_FL_heatmap$Z)])
min(MIF_tad_FL_heatmap$Z[!is.na(MIF_tad_FL_heatmap$Z)])
calculate_insulation_score("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                        "H1_FL_merged",
                                        10000,
                                        2000000)

MIF_tad_RNaseTreat_heatmap = parsing_MIF_matrix("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                     "H1_RNase_merged",
                                     10000,
                                     2000000)
max(MIF_tad_RNaseTreat_heatmap$Z[!is.na(MIF_tad_RNaseTreat_heatmap$Z)])
min(MIF_tad_RNaseTreat_heatmap$Z[!is.na(MIF_tad_RNaseTreat_heatmap$Z)])
calculate_insulation_score("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", 
                                                "H1_RNase_merged",
                                                10000,
                                                2000000)


plot_MIF_tad_heatmap <- function(data, output_file){
  colorbar_limits = c(0.2,5.3)
  
  p <- ggplot(data, aes(X, rev(Y), fill= Z)) + 
    geom_tile() +
    scale_fill_gradientn(
      colours = c("blue", "yellow", "red", "black"),
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = colorbar_limits
    ) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.title = element_blank())
  
  pdf(output_file, height = 5, width = 5)
  print(p)
  dev.off()
}

plot_MIF_tad_heatmap(data=MIF_tad_control_heatmap, "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/MIF_tad_control_heatmap.pdf")
plot_MIF_tad_heatmap(data=MIF_tad_NH4OAc_heatmap, "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/MIF_tad_NH4OAc_heatmap.pdf")
plot_MIF_tad_heatmap(data=MIF_tad_FL_heatmap, "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/MIF_tad_FL_heatmap.pdf")
plot_MIF_tad_heatmap(data=MIF_tad_RNaseTreat_heatmap, "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/MIF_tad_RNase_heatmap.pdf")


##### Plot distribution of insulation scores per sample
generate_insulation_score_data_per_sample <- function(directory, # with trailing slash
                                                      sample){
  boundaries_with_scores = data.frame(V1 = character(),
                                      V2 = integer(),
                                      V3 = double())

  for (i in 1:length(hg38_chromosomes)){
    if(file.exists(paste0(directory,"tads/",sample,"/TADs_",hg38_chromosomes[i],"/10000_boundary_insulation_score.txt"))){
      input_file = read.table(paste0(directory,"tads/",sample,"/TADs_",hg38_chromosomes[i],"/10000_boundary_insulation_score.txt"), stringsAsFactors = F)
      input_file = input_file[which(!is.na(input_file$V3)),]
      boundaries_with_scores = rbind(boundaries_with_scores, input_file)
    }
  }
  return(boundaries_with_scores)
}

boundaries_with_scores_control = generate_insulation_score_data_per_sample("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", "H1_control_merged")
boundaries_with_scores_NH4OAc = generate_insulation_score_data_per_sample("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", "H1_NH4OAc_merged")
boundaries_with_scores_FL = generate_insulation_score_data_per_sample("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", "H1_FL_merged")
boundaries_with_scores_RNase = generate_insulation_score_data_per_sample("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/", "H1_RNase_merged")

boundaries_with_scores_control = cbind(boundaries_with_scores_control, "Control")
boundaries_with_scores_NH4OAc = cbind(boundaries_with_scores_NH4OAc, "NH4OAc")
boundaries_with_scores_FL = cbind(boundaries_with_scores_FL, "FL")
boundaries_with_scores_RNase = cbind(boundaries_with_scores_RNase, "RNase")

colnames(boundaries_with_scores_control) = c("chr","coord","insulation_score","sample")
colnames(boundaries_with_scores_NH4OAc) = c("chr","coord","insulation_score","sample")
colnames(boundaries_with_scores_FL) = c("chr","coord","insulation_score","sample")
colnames(boundaries_with_scores_RNase) = c("chr","coord","insulation_score","sample")

all_boundary_file = rbind(boundaries_with_scores_control,
                          boundaries_with_scores_FL,
           boundaries_with_scores_NH4OAc,
           boundaries_with_scores_RNase)

#all_boundary_file$sample = factor(all_boundary_file$sample, levels = c("Control","NH4OAc","FL","RNase"))

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/TAD_boundary_insulation_score_boxplot.png", width = 15, height = 8, units = "in", res = 200)
p1<-ggplot(all_boundary_file, aes(x=sample, y=insulation_score)) + 
  geom_boxplot() +
  labs(x="", y="TAD boundary insulation score") +
  #ylim(0,50) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm"))
p2<-ggplot(all_boundary_file, aes(x=sample, y=insulation_score)) + 
  geom_boxplot() +
  labs(x="", y="TAD boundary insulation score") +
  ylim(0,50) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm"))
plot_grid(p1,p2)
dev.off()

wilcox.test(boundaries_with_scores_control$insulation_score,
            boundaries_with_scores_NH4OAc$insulation_score,
            paired = F)

wilcox.test(boundaries_with_scores_control$insulation_score,
            boundaries_with_scores_FL$insulation_score,
            paired = F)

wilcox.test(boundaries_with_scores_control$insulation_score,
            boundaries_with_scores_RNase$insulation_score,
            paired = F)

############ Generate datasets of union of all the boundaries across the 4 samples. If a boundary is not shared among all
# the 4 samples, the insulation score at that coordinate will be recalculated for the sample where it was not a boundary.

all_boundary_file = rbind(boundaries_with_scores_control,
                          boundaries_with_scores_NH4OAc,
                          boundaries_with_scores_FL,
                          boundaries_with_scores_RNase)

write.table(all_boundary_file,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/all_boundary_file.txt", row.names = F, col.names = F, sep = "\t", quote = F)

all_boundary_file[which(all_boundary_file[,1]==union_boundary_file[1,1] &
                          all_boundary_file[,2]==union_boundary_file[1,2]),]

union_boundary_file[which(union_boundary_file[,1]==union_boundary_file[1,1] &
                            union_boundary_file[,2]==union_boundary_file[1,2]),]


# This function was run separately in union_tad_boundaries_script.r
generate_full_boundary_dataset <- function(all_boundary_file_path,
                                           directory,
                                           bin_size,
                                           template_dim){
  
  all_boundary_file = read.table(all_boundary_file_path, stringsAsFactors = F)
  colnames(all_boundary_file) = c("chr","coord","insulation_score","sample")
  union_boundary_file = all_boundary_file[!duplicated(all_boundary_file[,1:2]),]
  union_boundary_file$chr = factor(union_boundary_file$chr, levels = hg38_chromosomes)
  
  union_boundary_file = union_boundary_file[order(union_boundary_file$chr),] # sorting is useful so that we can use the chr_table_cumulative to load contact matrices only once per chromosome
  
  chr_table = table(as.factor(union_boundary_file$chr))
  chr_table_cumulative = chr_table
  for (i in 1:length(chr_table)){
    chr_table_cumulative[i] = sum(chr_table[1:i])
  }
  
  union_boundary_file[,1] = as.character(union_boundary_file[,1])
  union_boundary_file = union_boundary_file[1:100,] # Line for testing purposes
  
  output_file = matrix(NA, nrow = nrow(union_boundary_file), ncol = 10)
  colnames(output_file) = c("chr", "coord", "Ctrl", "Ctrl_insulation", "NH4OAc", "NH4OAc_insulation", "FL", "FL_insulation", "RNase", "RNase_insulation")
  output_file[,1] = as.character(union_boundary_file[,1])
  output_file[,2] = union_boundary_file[,2]
  sample_names = c("Ctrl", "NH4OAc", "FL", "RNase")
  
  ### Initialize contact matrices for 4 samples for chromosome 1
  if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_control_merged/", bin_size, "/chr1_chr1_",bin_size,".txt")) != 0){
    control_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_control_merged/", bin_size, "/chr1_chr1_",bin_size,".txt"), stringsAsFactors = F)
  } else {
    control_matrix_sparse = F
  }
  
  if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_NH4OAc_merged/", bin_size, "/chr1_chr1_",bin_size,".txt")) != 0){
    NH4OAc_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_NH4OAc_merged/", bin_size, "/chr1_chr1_",bin_size,".txt"), stringsAsFactors = F)
  } else {
    NH4OAc_matrix_sparse = F
  }
  
  if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_FL_merged/", bin_size, "/chr1_chr1_",bin_size,".txt")) != 0){
    FL_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_FL_merged/", bin_size, "/chr1_chr1_",bin_size,".txt"), stringsAsFactors = F)
  } else {
    FL_matrix_sparse = F
  }
  
  if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_RNaseTreat_merged/", bin_size, "/chr1_chr1_",bin_size,".txt")) != 0){
    RNase_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_RNaseTreat_merged/", bin_size, "/chr1_chr1_",bin_size,".txt"), stringsAsFactors = F)
  } else {
    RNase_matrix_sparse = F
  }
  current_chr = "chr1" # just for printing to console for testing purposes
  
  #### Loop over the union set of TAD boundaries
  for(i in 1:nrow(union_boundary_file)){
    temp_boundary = union_boundary_file[i,]
    temp_boundary_all = all_boundary_file[which(all_boundary_file[,1]==as.character(temp_boundary[1]) &
                                                  all_boundary_file[,2]==as.numeric(temp_boundary[2])),]
    
    output_file[i,as.character(temp_boundary_all$sample)] = 1
    output_file[i, setdiff(sample_names, as.character(temp_boundary_all$sample))] = 0
    output_file[i,paste0(as.character(temp_boundary_all$sample),"_insulation")] = as.numeric(temp_boundary_all$insulation_score)
    
    ### Check if the boundary was not present in a sample, in that case calculate the insulation score for that location
    sample_with_no_boundary = names(which(output_file[i,sample_names]=="0"))
    if (length(sample_with_no_boundary) > 0){
      for (sample in sample_with_no_boundary){
        # Assign to input_matrix_sparse the specific pre-loaded matrix depending on the sample
        if (sample == "Ctrl"){
          input_matrix_sparse = control_matrix_sparse
        } else if (sample == "NH4OAc") {
          input_matrix_sparse = NH4OAc_matrix_sparse
        } else if (sample == "FL") {
          input_matrix_sparse = FL_matrix_sparse
        } else if (sample == "RNase") {
          input_matrix_sparse = RNase_matrix_sparse
        }
        
        if (!is.logical(input_matrix_sparse)){
          input_matrix_sparse[which(is.na(input_matrix_sparse[,3])),3] = 0 # removing NaN values
          input_matrix_sparse_bin = input_matrix_sparse
          input_matrix_sparse_bin[,1] = input_matrix_sparse_bin[,1] / bin_size + 1
          input_matrix_sparse_bin[,2] = input_matrix_sparse_bin[,2] / bin_size + 1
          
          boundary_bin = temp_boundary[1,2] / bin_size + 1 # bin location of the "boundary" (not boundary for this sample)
          shift = template_dim/2
          shift_bin = shift / bin_size
          
          upstream_bin = as.numeric(boundary_bin - shift_bin) 
          downstream_bin = as.numeric(boundary_bin + shift_bin - 1) 
          input_matrix_sparse_bin_template = input_matrix_sparse_bin[which(input_matrix_sparse_bin[,1] >= upstream_bin &
                                                                             input_matrix_sparse_bin[,1] <= downstream_bin &
                                                                             input_matrix_sparse_bin[,2] >= upstream_bin &
                                                                             input_matrix_sparse_bin[,2] <= downstream_bin),]
          
          input_matrix_sparse_bin_template_scaled <- input_matrix_sparse_bin_template
          input_matrix_sparse_bin_template_scaled[,1] = input_matrix_sparse_bin_template_scaled[,1] - upstream_bin + 1
          input_matrix_sparse_bin_template_scaled[,2] = input_matrix_sparse_bin_template_scaled[,2] - upstream_bin + 1
          input_matrix_sparse_bin_template_scaled = as.matrix(input_matrix_sparse_bin_template_scaled)
          
          # Calculate insulation score
          mean_L = sum(input_matrix_sparse_bin_template_scaled[which(input_matrix_sparse_bin_template_scaled[,1] <= shift_bin &
                                                                       input_matrix_sparse_bin_template_scaled[,2] <= shift_bin),3]) / (shift_bin * (shift_bin + 1) / 2)
          mean_R = sum(input_matrix_sparse_bin_template_scaled[which(input_matrix_sparse_bin_template_scaled[,1] > shift_bin &
                                                                       input_matrix_sparse_bin_template_scaled[,2] > shift_bin),3]) / (shift_bin * (shift_bin + 1) / 2)
          mean_X = sum(input_matrix_sparse_bin_template_scaled[which(input_matrix_sparse_bin_template_scaled[,1] <= shift_bin &
                                                                       input_matrix_sparse_bin_template_scaled[,2] > shift_bin),3]) / (shift_bin^2)
          insulation_score = max(mean_L, mean_R) / mean_X
          
          output_file[i,paste0(sample,"_insulation")] = insulation_score
          }
      }
    }
    
    print(paste0(i,"_",current_chr))
    ### Load new contact matrices for the 4 samples if we reach the final index of the current chromosome in chr_table_cumulative
    for (my_chr in 1:length(chr_table_cumulative)){
      if (i == as.numeric(chr_table_cumulative[my_chr]) & i <= as.numeric(chr_table_cumulative["chrX"])){
        current_chr = names(chr_table_cumulative[my_chr+1])
        ### Initialize contact matrices for 4 samples for chromosome 1
        if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_control_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt")) != 0){
          control_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_control_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt"), stringsAsFactors = F)
        } else {
          control_matrix_sparse = F
        }
        
        if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_NH4OAc_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt")) != 0){
          NH4OAc_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_NH4OAc_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt"), stringsAsFactors = F)
        } else {
          NH4OAc_matrix_sparse = F
        }
        
        if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_FL_merged/",bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt")) != 0){
          FL_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_FL_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt"), stringsAsFactors = F)
        } else {
          FL_matrix_sparse = F
        }
        
        if (file.size(paste0(directory, "HiC_contact_matrices/KR/H1_RNaseTreat_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt")) != 0){
          RNase_matrix_sparse = read.table(paste0(directory, "HiC_contact_matrices/KR/H1_RNaseTreat_merged/", bin_size, "/",current_chr,"_",current_chr,"_",bin_size,".txt"), stringsAsFactors = F)
        } else {
          RNase_matrix_sparse = F
        }
      }
    }
  }
  write.table(output_file,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/result/dataset_union_boundaries.txt", row.names = F, col.names = T, sep = "\t", quote = F)
  return(output_file)
}

# temp = generate_full_boundary_dataset("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/TADs/result/all_boundary_file.txt",
#                                       "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/",
#                                       10000,
#                                       2000000)

dataset_union_boundaries = read.table("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/dataset_union_boundaries.txt", stringsAsFactors = F, header = T)
      
table(as.factor(dataset_union_boundaries$Control))
table(as.factor(dataset_union_boundaries$FL))

dataset_union_boundaries.melt = melt(dataset_union_boundaries, id.vars = c("chr","coord"), measure.vars = c("Control_insulation","NH4OAc_insulation","FL_insulation","RNase_insulation"))
dataset_union_boundaries.melt$variable = gsub("_insulation","",dataset_union_boundaries.melt$variable)
#dataset_union_boundaries.melt$variable = factor(dataset_union_boundaries.melt$variable, c("Control","NH4OAc","FL","RNase"))

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/union_TAD_boundary_insulation_score_boxplot.png", width = 15, height = 8, units = "in", res = 200)
p1<-ggplot(dataset_union_boundaries.melt, aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(x="", y="Union TAD boundary insulation score") +
  #ylim(0,50) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm"))
p2<-ggplot(dataset_union_boundaries.melt, aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(x="", y="Union TAD boundary insulation score") +
  ylim(0,50) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm"))
plot_grid(p1,p2)
dev.off()

wilcox.test(dataset_union_boundaries$Ctrl_insulation,
            dataset_union_boundaries$NH4OAc_insulation,
            paired = F)

wilcox.test(dataset_union_boundaries$Ctrl_insulation,
            dataset_union_boundaries$FL_insulation,
            paired = F)

wilcox.test(dataset_union_boundaries$Ctrl_insulation,
            dataset_union_boundaries$RNase_insulation,
            paired = F)


######### Upset plot
library(UpSetR)
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/result/union_TAD_boundary_upset_plot.png", width = 9, height = 5, units = "in", res = 200)
upset(dataset_union_boundaries, 
      sets = c("RNase","NH4OAc","FL","Control"),
      sets.x.label = "Number of TAD boundaries",
      order.by = "freq",
      keep.order = T,
      text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
      point.size = 3, line.size = 1)
dev.off()






