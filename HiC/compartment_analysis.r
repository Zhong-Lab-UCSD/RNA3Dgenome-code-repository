############## RNA-genome interaction project compartment analysis
options(scipen=999)

chrSize = read.table("/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes")

hg38_chromosomes = as.character(chrSize$V1)
hg38_lengths = as.numeric(chrSize$V2)
names(hg38_lengths) = hg38_chromosomes

### Make gene density files
annotation <- read.table('/dataOS/rcalandrelli/MARGI/Homo_sapiens.GRCh38.84.chr.gtf_to_geneTable.tsv', stringsAsFactors = F)
colnames(annotation)=annotation[1,]
annotation=annotation[-1,]
annotation$start = as.numeric(annotation$start)
annotation$end = as.numeric(annotation$end)

Gr_annotation <- GRanges(
  seqnames = Rle(annotation[,1]),
  ranges = IRanges(as.numeric(annotation[,2]), end = as.numeric(annotation[,3]), names = c(1:nrow(annotation))),
  strand = Rle(strand(annotation[,4])),
  gene_id = annotation[,5],
  gene_name = annotation[,6],
  gene_biotype = annotation[,7],
  transcript_id = annotation[,8])

Gr_hg38 <- GRanges(
  seqnames = Rle(names(hg38_lengths)),
  ranges = IRanges(rep(1,24), end = as.numeric(hg38_lengths), names = c(1:length(hg38_lengths))),
  strand = Rle(strand('*')))
seqlengths(Gr_hg38) <- hg38_lengths[names(seqlengths(Gr_hg38))]

make_gene_density_file <- function(annotation_gr, genome_gr, window_size) {
  genome_window <- tileGenome(seqinfo(genome_gr), tilewidth = as.numeric(window_size), cut.last.tile.in.chrom = T)
  overlaps = countOverlaps(genome_window,Gr_annotation)
  mcols(genome_window)['overlap'] = overlaps
  df = data.frame(genome_window)
  df$start = df$start - 1
  df = df[,c(1,2,3,6)]
  write.table(df,paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/gene_density_annotation/gene_density_",window_size,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
}

make_gene_density_file(Gr_annotation, Gr_hg38, "5000000")
make_gene_density_file(Gr_annotation, Gr_hg38, "2500000")
make_gene_density_file(Gr_annotation, Gr_hg38, "1000000")
make_gene_density_file(Gr_annotation, Gr_hg38, "500000")


### Make single file for all eigen value across all genome
make_eigen_all_chromosomes <- function(directory,
                                       sample,
                                       resolution # string
                                    ){
  gene_density = read.table(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/gene_density_annotation/gene_density_",resolution,".txt"), stringsAsFactors = F)
  eigen_all = data.frame()
  for (i in hg38_chromosomes){
    gene_density_chr = gene_density[which(gene_density$V1 == i),]
    gene_density_norm = gene_density_chr$V4 - mean(gene_density_chr$V4)
    temp_eigen = read.table(paste0(directory,"compartments/",sample,"/",resolution,"/eigen_",i,".txt"), stringsAsFactors = F)
    temp_eigen[which(is.na(temp_eigen)),1] = 0
    #print(cor(temp_eigen,gene_density_norm))
    if (cor(temp_eigen,gene_density_norm)<0){
      temp_eigen = temp_eigen*(-1)
    }
    chr = i
    coord = seq(0, nrow(temp_eigen)*as.numeric(resolution) - as.numeric(resolution), by = as.numeric(resolution))
    eigen = temp_eigen
    temp_eigen_mat = data.frame(chr, coord, eigen)
    
    #temp_eigen = read.table(paste0(directory,"compartments/H1_",sample,"_merged/",resolution,"/eigen_",i,".txt"), stringsAsFactors = F)
    #temp_eigen1 = read.table(paste0("/mnt/extraids/OceanStor-SysCmn-5/wenxingzhao/project/8_phase_separation/3_compartment/1_juicer_eigenvec/HiC_H1_",sample,"/",resolution,"/",i,"_KR_eigen_",resolution,".txt"), stringsAsFactors = F)

    eigen_all = rbind(eigen_all, temp_eigen_mat)
    #eigen_all1 = c(eigen_all1, t(temp_eigen1))
  }
  colnames(eigen_all)[3] = "eigen"
  write.table(eigen_all, paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/",gsub("_merged","",sample),"_",resolution,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  return(eigen_all)
}

####### 2.5mb
eigen_control_2500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                   "H1_control_merged",
                                                   "2500000")
eigen_NH4OAc_2500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                   "H1_NH4OAc_merged",
                                                   "2500000")
eigen_FL_2500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                   "H1_FL_merged",
                                                   "2500000")
eigen_RNase_2500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                   "H1_RNase_merged",
                                                   "2500000")

####### 1mb
eigen_control_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                   "H1_control_merged",
                                                   "1000000")
eigen_NH4OAc_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                  "H1_NH4OAc_merged",
                                                  "1000000")
eigen_FL_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                              "H1_FL_merged",
                                              "1000000")
eigen_RNase_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                      "H1_RNase_merged",
                                                      "1000000")

eigen_RNaseCtrl_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                 "H1_RNaseCtrl_merged",
                                                 "1000000")

eigen_HFF_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                             "HFF",
                             "1000000")

eigen_K562_1000000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                               "K562",
                                               "1000000")


####### 500kb
eigen_control_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                   "H1_control_merged",
                                                   "500000")
eigen_NH4OAc_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                  "H1_NH4OAc_merged",
                                                  "500000")
eigen_FL_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                              "H1_FL_merged",
                                              "500000")
eigen_RNase_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                      "H1_RNase_merged",
                                                      "500000")

eigen_RNaseCtrl_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                "H1_RNaseCtrl_merged",
                                                "500000")

eigen_HFF_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                               "HFF",
                                               "500000")

eigen_K562_500000 = make_eigen_all_chromosomes("/dataOS/rcalandrelli/phase_separation/HiC/",
                                                "K562",
                                                "500000")

###################### Plot OE/PCC heatmaps and eigen tracks to make sure that there was not any wrong flipping due to a correlation with gene density track around 0 (but in a sample positve, in another negative).
# These false sign flippings are clearly noticeable if the eigen track is completely flipped over across different samples
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
require(dplyr)
library(igraph)

edge_to_matrix <- function(df,attrname="V3"){
  # df is source target wrights 
  g=graph.data.frame(df,directed = FALSE)
  matrix <- as_adjacency_matrix(g,type = "both", names = TRUE, sparse = FALSE, attr = attrname)
  return(matrix)
}

make_heatmap_eigen_plot <- function(directory="/dataOS/rcalandrelli/phase_separation/HiC/",
                                    plot_type, # "pcc" or "oe"
                                    sample,
                                    resolution,
                                    chr,
                                    fig.title="Title",
                                    force_manual_switch=F){
  # draw a combine plot 
  # hicmap = matrix(c(1,NA,NA,NA,NA,1,5,6,3,5,1,7,4,6,7,1),nrow = 4, ncol=4, byrow = T)
  # hicmap[is.na(hicmap)] = 0
  # hicmap_pcc <- cor(hicmap, method = "pearson")
  
  ### Heatmap
  hicmap <- read.table(paste0(directory,"HiC_contact_matrices/oe_KR/",sample,"/",resolution,"/",chr,"_",chr,"_",resolution,".txt")) %>% edge_to_matrix()
  # Even though hicmap is a dense matrix there are still row and col indexes brought from the sparse matrix. This means that there is a gap in 
  # rownames and colnames over the centromere that is just plotted as "empty space" or "no data space", which is not NA.
  
  if (plot_type == "pcc"){
    hicmap[is.na(hicmap)] = 0
    hicmap_pcc <- cor(hicmap, method = "pearson") # pearson correlation matrix of the OE
    #entropy = entropySpecificity(hicmap_pcc)
    m <- melt(hicmap_pcc)
    temp1 = m[which(is.na(m$value)),]
    m$Var2 = as.numeric(rownames(hicmap)[nrow(hicmap)]) - m$Var2 # to make the (0,0) vertex on the upper-left corner and have a classic heatmap visualization
    
    hm <- ggplot(data = m, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", limits=c(-1,1), na.value = "#969696") + #na.value = "#969696"
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.title.x = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position="none")
  } 
  else if (plot_type == "oe"){
    #entropy = entropySpecificity(hicmap)
    m <- melt(hicmap)
    #m[which(m$Var1==11000000 & m$Var2==41000000),]
    m$Var2 = as.numeric(rownames(hicmap)[nrow(hicmap)]) - m$Var2 # to make the (0,0) vertex on the upper-left corner and have a classic heatmap visualization
    #m[which(m$value==1),"value"] = 0
    m$value = log(m$value)
    #m[which(is.infinite(m$value)),"value"] = 0
    
    hm <- ggplot(data = m, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() +
      #scale_fill_distiller(name = "Legend title", palette = "Reds", direction = 1, na.value = "transparent") +
      scale_fill_gradient2(low = "blue", high = "red", na.value = "#969696") +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.title.x = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position="none")
  }
  
  ### Eigenvector
  gene_density = read.table(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/gene_density_annotation/gene_density_",resolution,".txt"), stringsAsFactors = F)
  gene_density_chr = gene_density[which(gene_density$V1 == chr),]
  gene_density_norm = gene_density_chr$V4 - mean(gene_density_chr$V4)
  
  ev <- read.table(paste0(directory,"compartments/",sample,"/",resolution,"/eigen_",chr,".txt"), header = F)
  
  temp_ev = ev
  temp_ev[which(is.na(temp_ev)),1] = 0
  if (cor(temp_ev,gene_density_norm)<0){
    ev = ev*(-1)
  }
  
  if (force_manual_switch){
    ev = -1*ev
  }
  
  ev %>% as.data.frame()
  colnames(ev) <- "V1"
  
  ev$x <- seq(1,length(ev$V1))
  bar_ev <- ggplot()+
    geom_bar(data=ev,aes(x=x,y=V1,color=V1>0,fill=V1>0),stat = "identity")+ 
    scale_fill_manual(name="compartment",labels=c("A","B"),values = c("#6b6665","#a31500"))+
    scale_color_manual(name="compartment",labels=c("A","B"),values = c("#6b6665","#a31500"))+
    xlab("Bins") +
    theme_bw() + 
    ylab("EV1") +
    #ylim(-1,1) +
    theme(axis.text.x = element_text(size = 12), 
          axis.title.x = element_text(size = 12, margin = margin(5,0,0,0)),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position="none")
  
  ### Entropy track
  # ev$V2 = NaN
  # ev[ev[which(!is.na(ev$V1)),"x"],"V2"] = entropy
  # bar_entropy <- ggplot()+
  #   geom_bar(data=ev,aes(x=x,y=V2), color="blue", fill = "blue", stat="identity") + 
  #   # scale_fill_manual(name="compartment",labels=c("A","B"),values = c("#6b6665","#a31500"))+
  #   # scale_color_manual(name="compartment",labels=c("A","B"),values = c("#6b6665","#a31500"))+
  #   xlab("Bins") +
  #   theme_bw() + 
  #   ylab("Entropy") +
  #   #ylim(-1,1) +
  #   theme(axis.text.x = element_text(size = 12), 
  #         axis.title.x = element_text(size = 12, margin = margin(5,0,0,0)),
  #         panel.border = element_blank(), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),panel.background = element_blank(),
  #         axis.title.y = element_text(size = 12, margin = margin(5,0,0,0)), 
  #         axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank(), 
  #         legend.position="none")
  
  # grob
  grob.title <- textGrob(fig.title, hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 20))
  p<-grid.arrange(hm, bar_ev,
                  nrow = 2, ncol =1,
                  heights = c(30, 10), top = grob.title)
  #return(list(p,entropy))
  return(p)
}

### Make plots
for (my_chr in hg38_chromosomes){
  my_plot_type = "pcc"
  #my_resolution = "2500000"
  #my_resolution = "1000000"
  my_resolution = "500000"
  g1 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                                plot_type = my_plot_type,
                                sample = "H1_control_merged",
                                resolution = my_resolution,
                                chr = my_chr,
                                fig.title = paste0("Control_",my_chr))
  
  g2 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                                plot_type = my_plot_type,
                                sample = "H1_NH4OAc_merged",
                                resolution = my_resolution,
                                chr = my_chr,
                                fig.title = paste0("NH4OAc_",my_chr))
  
  g3 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                                plot_type = my_plot_type,
                                sample = "H1_FL_merged",
                                resolution = my_resolution,
                                chr = my_chr,
                                fig.title = paste0("FL_",my_chr))
  
  g4 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                                plot_type = my_plot_type,
                                sample = "H1_RNase_merged",
                                resolution = my_resolution,
                                chr = my_chr,
                                fig.title = paste0("Rnase_",my_chr))
  
  g5 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                                plot_type = my_plot_type,
                                sample = "H1_RNaseCtrl_merged",
                                resolution = my_resolution,
                                chr = my_chr,
                                fig.title = paste0("RnaseCtrl_",my_chr))
  
  png(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/",my_chr,"_",my_plot_type,"_heatmap_eigen_",my_resolution,".png"), width = 17.5, height = 5, units = "in", res = 200)
  grid.arrange(g1,g2,g3,g4,g5,nrow=1)
  dev.off()
}

### Calculate entropy data
for (my_plot_type in c("oe","pcc")){
  for (my_resolution in c("500000","1000000")){
    entropy_list = list()
    for (my_chr in hg38_chromosomes){
      g1 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/",
                                    plot_type = my_plot_type,
                                    sample = "H1_control_merged",
                                    resolution = my_resolution,
                                    chr = my_chr,
                                    fig.title = paste0("Control_",my_chr))
      
      g2 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/",
                                    plot_type = my_plot_type,
                                    sample = "H1_NH4OAc_merged",
                                    resolution = my_resolution,
                                    chr = my_chr,
                                    fig.title = paste0("NH4OAc_",my_chr))
      
      g3 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/",
                                    plot_type = my_plot_type,
                                    sample = "H1_FL_merged",
                                    resolution = my_resolution,
                                    chr = my_chr,
                                    fig.title = paste0("FL_",my_chr))
      
      g4 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/",
                                    plot_type = my_plot_type,
                                    sample = "H1_RNaseTreat_merged",
                                    resolution = my_resolution,
                                    chr = my_chr,
                                    fig.title = paste0("Rnase_",my_chr))
      
      # E1 = mean(g1[[2]][which(!is.na(g1[[2]]))])
      # E2 = mean(g2[[2]][which(!is.na(g2[[2]]))])
      # E3 = mean(g3[[2]][which(!is.na(g3[[2]]))])
      # E4 = mean(g4[[2]][which(!is.na(g4[[2]]))])
      # entropy_list[[my_chr]] = c(E1,E2,E3,E4)
      
      temp_list = list(g1[[2]], g2[[2]], g3[[2]], g4[[2]])
      temp_out = do.call("cbind", lapply(temp_list, function(x) x[match(names(temp_list[[1]]), names(x))]))
      temp_out = cbind(my_chr, rownames(temp_out), temp_out)
      colnames(temp_out) = c("chr","coord","Control","NH4OAc","FL", "RNase")
      entropy_list[[my_chr]] = temp_out
    }
    out = do.call("rbind",entropy_list)
    #colnames(out) = c("chr","coord","Control","NH4OAc","FL", "RNase")
    write.table(out, paste0("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy/entropy_",my_plot_type,"_",my_resolution,"_all_rows.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  }
}

### Entropy data analysis

entropy_oe_500000 = read.table("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy/entropy_oe_500000_all_rows.txt", header = T)
entropy_oe_1000000 = read.table("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy/entropy_oe_1000000_all_rows.txt", header = T)
# entropy_pcc_500000 = read.table("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy/entropy_pcc_500000_all_rows.txt", header = T)
# entropy_pcc_1000000 = read.table("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy/entropy_pcc_1000000_all_rows.txt", header = T)

entropy_oe_500000 = entropy_oe_500000[which(!is.na(rowSums(entropy_oe_500000[,3:6]))),]
entropy_oe_1000000 = entropy_oe_1000000[which(!is.na(rowSums(entropy_oe_1000000[,3:6]))),]
# entropy_pcc_500000 = entropy_pcc_500000[which(!is.na(rowSums(entropy_pcc_500000[,3:6]))),]
# entropy_pcc_1000000 = entropy_pcc_1000000[which(!is.na(rowSums(entropy_pcc_1000000[,3:6]))),]

entropy_oe = rbind(cbind(melt(entropy_oe_500000, id.vars = c("chr","coord")),resolution="500000"),
                   cbind(melt(entropy_oe_1000000, id.vars = c("chr","coord")),resolution="1000000"))
# entropy_pcc = rbind(cbind(melt(entropy_pcc_500000, id.vars = c("chr","coord")),resolution="500000"),
#                    cbind(melt(entropy_pcc_1000000, id.vars = c("chr","coord")),resolution="1000000"))

png("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy_boxplot.png", width = 6, height = 6, units = "in", res = 200)
ggplot(entropy_oe, aes(x=variable, y=value, fill=resolution)) + 
  geom_boxplot() +
  labs(x="", y="Entropy") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("OE contact matrix")
dev.off()


# p2<-ggplot(entropy_pcc, aes(x=variable, y=value, fill=resolution)) + 
#   geom_boxplot() +
#   labs(x="", y="Entropy") +
#   #ylim(0,500) +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=20),
#         axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 20),
#         axis.ticks.length = unit(0.2, "cm")) +
#   ggtitle("Pearson correlation matrix")
# png("/dataOS/rcalandrelli/phase_separation/compartments/result/entropy_boxplots.png", width = 12, height = 6, units = "in", res = 200)
# plot_grid(p1,p2)
# dev.off()

wilcox.test(entropy_oe_500000$Control, entropy_oe_500000$NH4OAc)
wilcox.test(entropy_oe_500000$Control, entropy_oe_500000$FL)
wilcox.test(entropy_oe_500000$Control, entropy_oe_500000$RNase)

wilcox.test(entropy_oe_1000000$Control, entropy_oe_1000000$NH4OAc)
wilcox.test(entropy_oe_1000000$Control, entropy_oe_1000000$FL)
wilcox.test(entropy_oe_1000000$Control, entropy_oe_1000000$RNase)

wilcox.test(entropy_pcc_500000$Control, entropy_pcc_500000$NH4OAc)
wilcox.test(entropy_pcc_500000$Control, entropy_pcc_500000$FL)
wilcox.test(entropy_pcc_500000$Control, entropy_pcc_500000$RNase)

wilcox.test(entropy_pcc_1000000$Control, entropy_pcc_1000000$NH4OAc)
wilcox.test(entropy_pcc_1000000$Control, entropy_pcc_1000000$FL)
wilcox.test(entropy_pcc_1000000$Control, entropy_pcc_1000000$RNase)



### NB!!!!! Change sign again for eigen of chromosome/s and sample/s where the sign has been wrongly flipped

# 1Mb
eigen_RNase_1000000[which(eigen_RNase_1000000$chr == "chr19"),"eigen"] = -1*eigen_RNase_1000000[which(eigen_RNase_1000000$chr == "chr19"),"eigen"]
eigen_RNase_1000000[which(eigen_RNase_1000000$chr == "chr21"),"eigen"] = -1*eigen_RNase_1000000[which(eigen_RNase_1000000$chr == "chr21"),"eigen"]

eigen_RNaseCtrl_1000000[which(eigen_RNaseCtrl_1000000$chr == "chr19"),"eigen"] = -1*eigen_RNaseCtrl_1000000[which(eigen_RNaseCtrl_1000000$chr == "chr19"),"eigen"]


write.table(eigen_control_1000000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_1000000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_NH4OAc_1000000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_NH4OAc_1000000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_FL_1000000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_FL_1000000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_RNase_1000000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_1000000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_RNaseCtrl_1000000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNaseCtrl_1000000.txt",sep="\t",row.names = F, col.names = T, quote = F)


# 500kb
eigen_RNase_500000[which(eigen_RNase_500000$chr == "chr21"),"eigen"] = -1*eigen_RNase_500000[which(eigen_RNase_500000$chr == "chr21"),"eigen"]

write.table(eigen_control_500000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_NH4OAc_500000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_NH4OAc_500000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_FL_500000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_FL_500000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_RNase_500000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_500000.txt",sep="\t",row.names = F, col.names = T, quote = F)
write.table(eigen_RNaseCtrl_500000,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNaseCtrl_500000.txt",sep="\t",row.names = F, col.names = T, quote = F)


################ Summary of A/B compartments
resolution_vector = c(rep("1000000",5), rep("500000",5))#, rep("500000",4))
sample_vector = rep(c("Control","NH4OAc","FL","RNase","RNaseCtrl"),2)
A_bins = c(#length(which(eigen_control_5000000$eigen>0)), length(which(eigen_NH4OAc_5000000$eigen>0)), length(which(eigen_FL_5000000$eigen>0)), length(which(eigen_RNaseTreat_5000000$eigen>0)),
             length(which(eigen_control_1000000$eigen>0)), length(which(eigen_NH4OAc_1000000$eigen>0)), length(which(eigen_FL_1000000$eigen>0)), length(which(eigen_RNase_1000000$eigen>0)),length(which(eigen_RNaseCtrl_1000000$eigen>0)),
             length(which(eigen_control_500000$eigen>0)), length(which(eigen_NH4OAc_500000$eigen>0)), length(which(eigen_FL_500000$eigen>0)), length(which(eigen_RNase_500000$eigen>0)),length(which(eigen_RNaseCtrl_500000$eigen>0)))            
B_bins = c(#length(which(eigen_control_5000000$eigen<0)), length(which(eigen_NH4OAc_5000000$eigen<0)), length(which(eigen_FL_5000000$eigen<0)), length(which(eigen_RNaseTreat_5000000$eigen<0)),
             length(which(eigen_control_1000000$eigen<0)), length(which(eigen_NH4OAc_1000000$eigen<0)), length(which(eigen_FL_1000000$eigen<0)), length(which(eigen_RNase_1000000$eigen<0)),length(which(eigen_RNaseCtrl_1000000$eigen<0)),
             length(which(eigen_control_500000$eigen<0)), length(which(eigen_NH4OAc_500000$eigen<0)), length(which(eigen_FL_500000$eigen<0)), length(which(eigen_RNase_500000$eigen<0)),length(which(eigen_RNaseCtrl_500000$eigen<0)))   

df_summary_eigen = data.frame(resolution_vector,sample_vector,A_bins,B_bins)
colnames(df_summary_eigen)[1:2] = c("resolution","sample")
write.table(df_summary_eigen,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/summary_eigen.txt", row.names = F, col.names = T, sep = '\t', quote = F)

# Summary numbers
sum(df_summary_eigen[which(df_summary_eigen$resolution=="1000000"),"A_bins"])
sum(df_summary_eigen[which(df_summary_eigen$resolution=="1000000"),"B_bins"])

sum(df_summary_eigen[which(df_summary_eigen$resolution=="500000"),"A_bins"])
sum(df_summary_eigen[which(df_summary_eigen$resolution=="500000"),"B_bins"])

##### Summary of A/B compartment switch between control and treated
switch_stats <- function(x){
  if (x[1]==0 | x[2]==0){
    out = NA
  } else {
    if (x[1]>0 & x[2]>0){
      out = "AA"
    } else if (x[1]<0 & x[2]<0){
      out = "BB"
    } else if (x[1]>0 & x[2]<0){
      out = "AB"
    } else if (x[1]<0 & x[2]>0){
      out = "BA"
    }
  }
  return(out)
}

# 5mb
# full_eigen_data_5000000 = cbind(eigen_control_5000000,eigen_NH4OAc_5000000,eigen_FL_5000000,eigen_RNaseTreat_5000000)
# full_eigen_data_5000000 = full_eigen_data_5000000[,c(-4,-5,-7,-8,-10,-11)]
# colnames(full_eigen_data_5000000)[3:6] = c("eigen_control","eigen_NH4OAc","eigen_FL","eigen_RNase")
# 
# full_eigen_data_5000000$switch_NH4OAc <- apply(full_eigen_data_5000000[,c(3,4)],1,switch_stats)
# full_eigen_data_5000000$switch_FL <- apply(full_eigen_data_5000000[,c(3,5)],1,switch_stats)
# full_eigen_data_5000000$switch_RNaseTreat <- apply(full_eigen_data_5000000[,c(3,6)],1,switch_stats)
# 
# label_summary_5000000 = c(nrow(full_eigen_data_5000000[which(full_eigen_data_5000000$eigen_control>0),]),nrow(full_eigen_data_5000000[which(full_eigen_data_5000000$eigen_control<0),]),
#                           table(as.factor(full_eigen_data_5000000$switch_NH4OAc))["AB"],table(as.factor(full_eigen_data_5000000$switch_NH4OAc))["BA"],table(as.factor(full_eigen_data_5000000$switch_NH4OAc))["AA"],table(as.factor(full_eigen_data_5000000$switch_NH4OAc))["BB"],
#                           table(as.factor(full_eigen_data_5000000$switch_FL))["AB"],table(as.factor(full_eigen_data_5000000$switch_FL))["BA"],table(as.factor(full_eigen_data_5000000$switch_FL))["AA"],table(as.factor(full_eigen_data_5000000$switch_FL))["BB"],
#                           table(as.factor(full_eigen_data_5000000$switch_RNaseTreat))["AB"],table(as.factor(full_eigen_data_5000000$switch_RNaseTreat))["BA"],table(as.factor(full_eigen_data_5000000$switch_RNaseTreat))["AA"],table(as.factor(full_eigen_data_5000000$switch_RNaseTreat))["BB"])

# 1mb
full_eigen_data_1000000 = cbind(eigen_control_1000000,eigen_NH4OAc_1000000,eigen_FL_1000000,eigen_RNase_1000000,eigen_RNaseCtrl_1000000)
full_eigen_data_1000000 = full_eigen_data_1000000[,c(-4,-5,-7,-8,-10,-11,-13,-14)]
colnames(full_eigen_data_1000000)[3:7] = c("eigen_control","eigen_NH4OAc","eigen_FL","eigen_RNase","eigen_RNaseCtrl")

full_eigen_data_1000000$switch_NH4OAc <- apply(full_eigen_data_1000000[,c(3,4)],1,switch_stats)
full_eigen_data_1000000$switch_FL <- apply(full_eigen_data_1000000[,c(3,5)],1,switch_stats)
full_eigen_data_1000000$switch_RNase <- apply(full_eigen_data_1000000[,c(3,6)],1,switch_stats)
full_eigen_data_1000000$switch_RNaseCtrl <- apply(full_eigen_data_1000000[,c(3,7)],1,switch_stats)

label_summary_1000000 = c(nrow(full_eigen_data_1000000[which(full_eigen_data_1000000$eigen_control>0),]),nrow(full_eigen_data_1000000[which(full_eigen_data_1000000$eigen_control<0),]),
                          table(as.factor(full_eigen_data_1000000$switch_NH4OAc))["AB"],table(as.factor(full_eigen_data_1000000$switch_NH4OAc))["BA"],table(as.factor(full_eigen_data_1000000$switch_NH4OAc))["AA"],table(as.factor(full_eigen_data_1000000$switch_NH4OAc))["BB"],
                          table(as.factor(full_eigen_data_1000000$switch_FL))["AB"],table(as.factor(full_eigen_data_1000000$switch_FL))["BA"],table(as.factor(full_eigen_data_1000000$switch_FL))["AA"],table(as.factor(full_eigen_data_1000000$switch_FL))["BB"],
                          table(as.factor(full_eigen_data_1000000$switch_RNase))["AB"],table(as.factor(full_eigen_data_1000000$switch_RNase))["BA"],table(as.factor(full_eigen_data_1000000$switch_RNase))["AA"],table(as.factor(full_eigen_data_1000000$switch_RNase))["BB"],
                          table(as.factor(full_eigen_data_1000000$switch_RNaseCtrl))["AB"],table(as.factor(full_eigen_data_1000000$switch_RNaseCtrl))["BA"],table(as.factor(full_eigen_data_1000000$switch_RNaseCtrl))["AA"],table(as.factor(full_eigen_data_1000000$switch_RNaseCtrl))["BB"])

# 500kb
full_eigen_data_500000 = cbind(eigen_control_500000,eigen_NH4OAc_500000,eigen_FL_500000,eigen_RNase_500000,eigen_RNaseCtrl_500000)
full_eigen_data_500000 = full_eigen_data_500000[,c(-4,-5,-7,-8,-10,-11,-13,-14)]
colnames(full_eigen_data_500000)[3:7] = c("eigen_control","eigen_NH4OAc","eigen_FL","eigen_RNase","eigen_RNaseCtrl")

full_eigen_data_500000$switch_NH4OAc <- apply(full_eigen_data_500000[,c(3,4)],1,switch_stats)
full_eigen_data_500000$switch_FL <- apply(full_eigen_data_500000[,c(3,5)],1,switch_stats)
full_eigen_data_500000$switch_RNase <- apply(full_eigen_data_500000[,c(3,6)],1,switch_stats)
full_eigen_data_500000$switch_RNaseCtrl <- apply(full_eigen_data_500000[,c(3,7)],1,switch_stats)

label_summary_500000 = c(nrow(full_eigen_data_500000[which(full_eigen_data_500000$eigen_control>0),]),nrow(full_eigen_data_500000[which(full_eigen_data_500000$eigen_control<0),]),
                         table(as.factor(full_eigen_data_500000$switch_NH4OAc))["AB"],table(as.factor(full_eigen_data_500000$switch_NH4OAc))["BA"],table(as.factor(full_eigen_data_500000$switch_NH4OAc))["AA"],table(as.factor(full_eigen_data_500000$switch_NH4OAc))["BB"],
                         table(as.factor(full_eigen_data_500000$switch_FL))["AB"],table(as.factor(full_eigen_data_500000$switch_FL))["BA"],table(as.factor(full_eigen_data_500000$switch_FL))["AA"],table(as.factor(full_eigen_data_500000$switch_FL))["BB"],
                         table(as.factor(full_eigen_data_500000$switch_RNase))["AB"],table(as.factor(full_eigen_data_500000$switch_RNase))["BA"],table(as.factor(full_eigen_data_500000$switch_RNase))["AA"],table(as.factor(full_eigen_data_500000$switch_RNase))["BB"],
                         table(as.factor(full_eigen_data_500000$switch_RNaseCtrl))["AB"],table(as.factor(full_eigen_data_500000$switch_RNaseCtrl))["BA"],table(as.factor(full_eigen_data_500000$switch_RNaseCtrl))["AA"],table(as.factor(full_eigen_data_500000$switch_RNaseCtrl))["BB"])

### Make summary dataframe
resolution_vector = c(rep("1000000",18), rep("500000",18))
sample_vector = rep(c(rep("Control",2), rep("NH4OAc",4),rep("FL",4),rep("RNase",4),rep("RNaseCtrl",4)),2)
switch_label = rep(c("AA","BB",rep(c("AB","BA","AA","BB"),4)),2)
switch_vector = as.numeric(c(label_summary_1000000,label_summary_500000))

df_switch_summary = data.frame(resolution_vector, sample_vector, switch_label, switch_vector)
df_switch_summary$switch_label = factor(df_switch_summary$switch_label, levels = c("AB","BA","AA","BB"))
df_switch_summary$sample_vector = factor(df_switch_summary$sample_vector, levels = c("Control","NH4OAc","FL","RNase","RNaseCtrl"))

### Save data in a table
table_switch_summary = matrix(nrow=0,ncol=6)
for (i in c("1000000","500000")){
  temp_mat = df_switch_summary[which(df_switch_summary$resolution_vector == i),]
  temp_mat_1 = cbind(temp_mat[which(temp_mat$sample_vector=="NH4OAc"),"switch_vector"],
                     temp_mat[which(temp_mat$sample_vector=="FL"),"switch_vector"],
                     temp_mat[which(temp_mat$sample_vector=="RNase"),"switch_vector"],
                     temp_mat[which(temp_mat$sample_vector=="RNaseCtrl"),"switch_vector"])
  temp_mat_1 = cbind(i,c("AB","BA","AA","BB"),temp_mat_1)
  table_switch_summary = rbind(table_switch_summary,temp_mat_1)
}
colnames(table_switch_summary) = c("resolution","compartment_switch","NH4OAc","FL","RNase","RNaseCtrl")
write.table(table_switch_summary,"/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/summary_eigen_switch.txt", row.names = F, col.names = T, sep = '\t', quote = F)

######## Distribution of length of consecutive switches: we want to see how many are actually "meaningful" switches that comprise more than one bin

calculate_length_consecutive_switches <- function(eigen_data,
                                                  sample,
                                                  switch_type){
  i = 1
  out_list = list()
  while (i <= nrow(eigen_data)){
    
    if(!is.na(eigen_data[i,paste0("switch_", sample)]) & eigen_data[i,paste0("switch_", sample)] == switch_type){
      counter = 1
      k = i + 1
      while(k <= nrow(eigen_data)){
        if (!is.na(eigen_data[k,paste0("switch_", sample)]) & eigen_data[k,paste0("switch_", sample)] == switch_type){
          counter = counter + 1
          k = k + 1
        } else {
          out = c(as.character(eigen_data[i,1]), eigen_data[i,2], counter, switch_type)
          out_list[[i]] = out
          i = k
          break
        }
      }
    } else {
      i = i + 1
    }
  }
  out_df = do.call("rbind", out_list)
  colnames(out_df) = c("chr", "coord", "consecutive_switches", "switch_label")
  return(out_df)
}

##### 1 Mb
switch_length_distribution_1000000_list = list()
switch_length_distribution_1000000_list[["NH4OAc"]] = list()
switch_length_distribution_1000000_list[["FL"]] = list()
switch_length_distribution_1000000_list[["RNase"]] = list()
switch_length_distribution_1000000_list[["RNaseCtrl"]] = list()

for (x in c("NH4OAc","FL","RNase","RNaseCtrl")){
  for (y in c("AB", "BA")){
    switch_length_distribution_1000000_list[[x]][[y]] = calculate_length_consecutive_switches(eigen_data = full_eigen_data_1000000, sample = x, switch_type = y)
  }
}


summary(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"]) >= 3)
switch_length_distribution_1000000_list[["RNase"]][["AB"]][which(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"] >= 3),]

summary(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"]) >= 3)
switch_length_distribution_1000000_list[["RNase"]][["BA"]][which(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"] >= 3),]


### 500 kb
switch_length_distribution_500000_list = list()
switch_length_distribution_500000_list[["NH4OAc"]] = list()
switch_length_distribution_500000_list[["FL"]] = list()
switch_length_distribution_500000_list[["RNase"]] = list()
switch_length_distribution_500000_list[["RNaseCtrl"]] = list()

for (x in c("NH4OAc","FL","RNase","RNaseCtrl")){
  for (y in c("AB", "BA")){
    switch_length_distribution_500000_list[[x]][[y]] = calculate_length_consecutive_switches(eigen_data = full_eigen_data_500000, sample = x, switch_type = y)
  }
}

summary(as.numeric(switch_length_distribution_500000_list[["RNase"]][["AB"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_500000_list[["RNase"]][["AB"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]) > 4)

summary(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]) > 2)


##################################### 1 MB
switch_length_distribution_1000000_summary = matrix(data = 0, nrow=7, ncol=9)
switch_length_distribution_1000000_summary[,1] = seq(1,7)
colnames(switch_length_distribution_1000000_summary) = c("x",paste0(c("NH4OAc","FL","RNase","RNaseCtrl"),"_AB"),paste0(c("NH4OAc","FL","RNase","RNaseCtrl"),"_BA"))

for (i in c("NH4OAc","FL","RNase","RNaseCtrl")){
  for (j in c("AB", "BA")){
    temp = table(as.factor(as.numeric(switch_length_distribution_1000000_list[[i]][[j]][,"consecutive_switches"])))
    switch_length_distribution_1000000_summary[as.numeric(names(temp)),paste0(i,"_",j)] = as.numeric(temp)
  }
}

switch_length_distribution_1000000_summary = data.frame(switch_length_distribution_1000000_summary)

switch_length_distribution_1000000_summary[1,] / colSums(switch_length_distribution_1000000_summary)

df_AB_1000000 = melt(switch_length_distribution_1000000_summary[,1:5], id.vars = "x")
df_AB_1000000 = df_AB_1000000[which(df_AB_1000000$x %in% seq(1,max(df_AB_1000000[which(df_AB_1000000$value != 0),"x"]))),]

df_BA_1000000 = melt(switch_length_distribution_1000000_summary[,c(1,6:9)], id.vars = "x")
df_BA_1000000 = df_BA_1000000[which(df_BA_1000000$x %in% seq(1,max(df_BA_1000000[which(df_BA_1000000$value != 0),"x"]))),]

p1<-ggplot(data = df_AB_1000000, aes(x=x, y=value, color=variable)) +
  geom_point(size = 0.8) +
  geom_line() +
  xlab("Number of consecutive switches") +
  ylab("Frequency") +
  scale_x_continuous(labels = as.character(seq(1,max(df_AB_1000000$x))), breaks = seq(1,max(df_AB_1000000$x))) +
  ylim(c(0,80)) +
  scale_color_manual(labels = c("NH4OAc", "FL", "RNase", "RNaseCtrl"), values = c("red", "green","blue", "orange")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(size=8),
        legend.position = "none") +
  geom_vline(xintercept = 3, linetype = "dashed") +
        ggtitle("AB switches - 1 Mb")

p2<-ggplot(data = df_BA_1000000, aes(x=x, y=value, color=variable)) +
  geom_point(size = 0.8) +
  geom_line() +
  xlab("Number of consecutive switches") +
  ylab("Frequency") +
  scale_x_continuous(labels = as.character(seq(1,max(df_BA_1000000$x))), breaks = seq(1,max(df_BA_1000000$x))) +
  ylim(c(0,80)) +
  scale_color_manual(labels = c("NH4OAc", "FL", "RNase","RNaseCtrl"), values = c("red", "green","blue","orange")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(size=8),
        legend.position = "none") +
  geom_vline(xintercept = 3, linetype = "dashed") +
  ggtitle("BA switches - 1 Mb")

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/consecutive_switches_histogram_1000000.png", width = 5, height = 2.5, res = 300, units = "in")
plot_grid(p1,p2)
dev.off()


##################################### 500 kb
switch_length_distribution_500000_summary = matrix(data = 0, nrow=16, ncol=9)
switch_length_distribution_500000_summary[,1] = seq(1,16)
colnames(switch_length_distribution_500000_summary) = c("x",paste0(c("NH4OAc","FL","RNase","RNaseCtrl"),"_AB"),paste0(c("NH4OAc","FL","RNase","RNaseCtrl"),"_BA"))

for (i in c("NH4OAc","FL","RNase","RNaseCtrl")){
  for (j in c("AB", "BA")){
    temp = table(as.factor(as.numeric(switch_length_distribution_500000_list[[i]][[j]][,"consecutive_switches"])))
    switch_length_distribution_500000_summary[as.numeric(names(temp)),paste0(i,"_",j)] = as.numeric(temp)
  }
}

switch_length_distribution_500000_summary = data.frame(switch_length_distribution_500000_summary)

switch_length_distribution_500000_summary[1,] / colSums(switch_length_distribution_500000_summary)

df_AB_500000 = melt(switch_length_distribution_500000_summary[,1:5], id.vars = "x")
df_AB_500000 = df_AB_500000[which(df_AB_500000$x %in% seq(1,max(df_AB_500000[which(df_AB_500000$value != 0),"x"]))),]

df_BA_500000 = melt(switch_length_distribution_500000_summary[,c(1,6:9)], id.vars = "x")
df_BA_500000 = df_BA_500000[which(df_BA_500000$x %in% seq(1,max(df_BA_500000[which(df_BA_500000$value != 0),"x"]))),]

p1<-ggplot(data = df_AB_500000, aes(x=x, y=value, color=variable)) +
  geom_point(size = 0.8) +
  geom_line() +
  xlab("Number of consecutive switches") +
  ylab("Frequency") +
  scale_x_continuous(labels = as.character(seq(1,max(df_AB_500000$x), by=2)), breaks = seq(1,max(df_AB_500000$x),by=2)) +
  ylim(c(0,135)) +
  scale_color_manual(labels = c("NH4OAc", "FL", "RNase","RNaseCtrl"), values = c("red", "green","blue","orange")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(size = 8),
        legend.position = "none") +
  geom_vline(xintercept = 4, linetype = "dashed") +
  ggtitle("AB switches - 500 kb")

p2<-ggplot(data = df_BA_500000, aes(x=x, y=value, color=variable)) +
  geom_point(size = 0.8) +
  geom_line() +
  xlab("Number of consecutive switches") +
  ylab("Frequency") +
  scale_x_continuous(labels = as.character(seq(1,max(df_BA_500000$x),by = 2)), breaks = seq(1,max(df_BA_500000$x),by=2)) +
  ylim(c(0,135)) +
  scale_color_manual(labels = c("NH4OAc", "FL", "RNase","RNaseCtrl"), values = c("red", "green","blue","orange")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(size = 8),
        legend.position = "none") +
  geom_vline(xintercept = 3, linetype = "dashed") +
  ggtitle("BA switches - 500 kb")


png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/consecutive_switches_histogram_500000.png", width = 5, height = 2.5, res = 300, units = "in")
plot_grid(p1,p2)
dev.off()





####### Switch stats
df.table_switch_summary = data.frame(table_switch_summary[which(table_switch_summary[,"resolution"]=="500000"),])
df.table_switch_summary$NH4OAc = as.numeric(as.character(df.table_switch_summary$NH4OAc))
df.table_switch_summary$FL = as.numeric(as.character(df.table_switch_summary$FL))
df.table_switch_summary$RNase = as.numeric(as.character(df.table_switch_summary$RNase))

average_perc_switch_sample <- function(sample){
  return(sum(df.table_switch_summary[which(df.table_switch_summary$compartment_switch %in% c("AB","BA")),sample])/sum(df.table_switch_summary[,sample]))
}

average_perc_switch_sample("NH4OAc")
average_perc_switch_sample("FL")
average_perc_switch_sample("RNase")

switch_ratio_all <- function(switch_type){
  out = c()
  for (i in c("500000")){
    for (j in c("NH4OAc","FL","RNase")){
      switch_ratio = df.table_switch_summary[which(df.table_switch_summary$resolution == i &
                                                     df.table_switch_summary$compartment_switch == switch_type),j]/sum(df.table_switch_summary[which(df.table_switch_summary$resolution == i),j])
      out = c(out,switch_ratio)
    }
  }
  return(out)
}

AB_all_ratios = switch_ratio_all("AB")
BA_all_ratios = switch_ratio_all("BA")

t.test(BA_all_ratios-AB_all_ratios)

######### Chi-square tests

# to test the association between AB switches and treatment
temp = rbind(colSums(df.table_switch_summary[which(df.table_switch_summary$compartment_switch=="AB"),c("RNase","FL","NH4OAc")]),
      colSums(df.table_switch_summary[which(df.table_switch_summary$compartment_switch!="AB"),c("RNase","FL","NH4OAc")]))
chisq.test(temp)

# to test the association between BA switches and treatment
temp = rbind(colSums(df.table_switch_summary[which(df.table_switch_summary$compartment_switch=="BA"),c("RNase","FL","NH4OAc")]),
             colSums(df.table_switch_summary[which(df.table_switch_summary$compartment_switch!="BA"),c("RNase","FL","NH4OAc")]))
chisq.test(temp)

# to test the association between switches and treatment
temp = rbind(colSums(df.table_switch_summary[which(df.table_switch_summary$compartment_switch %in% c("AB","BA")),c("RNase","FL","NH4OAc")]),
             colSums(df.table_switch_summary[which(df.table_switch_summary$compartment_switch %in% c("AA","BB")),c("RNase","FL","NH4OAc")]))
chisq.test(temp)

### Table update to test association between RNase and NonRNase
chisquare_sample <- function(sample,
                             switch_type){
  
  df.table_switch_summary_mod = data.frame(cbind(as.character(df.table_switch_summary[,"compartment_switch"]), 
                                                 rowSums(df.table_switch_summary[,setdiff(c("FL","NH4OAc","RNase"), sample)]), 
                                                 df.table_switch_summary[,sample]))
  
  colnames(df.table_switch_summary_mod) = c("compartment_switch",paste0("Non",sample), sample)
  df.table_switch_summary_mod[,2] = as.numeric(as.character(df.table_switch_summary_mod[,2]))
  df.table_switch_summary_mod[,3] = as.numeric(as.character(df.table_switch_summary_mod[,3]))
  
  # to test the association between AB switches and RNase treatment
  if (switch_type == "AB"){
    temp = rbind(colSums(df.table_switch_summary_mod[which(df.table_switch_summary_mod$compartment_switch=="AB"),3:2]),
                 colSums(df.table_switch_summary_mod[which(df.table_switch_summary_mod$compartment_switch!="AB"),3:2]))
    out = chisq.test(temp)
    or = temp[1,1]*temp[2,2]/(temp[1,2]*temp[2,1])
  }
  else if (switch_type == "BA"){
    # to test the association between BA switches and RNase treatment
    temp = rbind(colSums(df.table_switch_summary_mod[which(df.table_switch_summary_mod$compartment_switch=="BA"),3:2]),
                 colSums(df.table_switch_summary_mod[which(df.table_switch_summary_mod$compartment_switch!="BA"),3:2]))
    out = chisq.test(temp)
    or = temp[1,1]*temp[2,2]/(temp[1,2]*temp[2,1])
  } 
  else if (switch_type == "both") {
    # to test the association between switches and RNase treatment
    temp = rbind(colSums(df.table_switch_summary_mod[which(df.table_switch_summary_mod$compartment_switch %in% c("AB","BA")),3:2]),
                 colSums(df.table_switch_summary_mod[which(df.table_switch_summary_mod$compartment_switch %in% c("AA","BB")),3:2]))
    out = chisq.test(temp)
    or = temp[1,1]*temp[2,2]/(temp[1,2]*temp[2,1])
  }
  return(list(out,or))
}

chisquare_sample("NH4OAc","both")
chisquare_sample("FL","both")
chisquare_sample("RNase","both")

chisquare_sample("RNase","AB")
chisquare_sample("RNase","BA")


### Plot data
make_barplot_switch <- function(df_switch_summary,resolution,label_size=10){
  df_switch_summary_res = df_switch_summary[which(df_switch_summary$resolution_vector == resolution),]
  p <- ggplot(df_switch_summary_res, aes(y=switch_vector, x=sample_vector, fill=switch_label)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c('#b2bc9b','#6ca7d2','#a81900','#716c6b')) +
    labs(x = "", y = "Number of cases") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size = label_size, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = label_size, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.title = element_blank()) +
    ggtitle(paste0("Resolution: ",resolution))
  return(p)
}

make_scatter_plot_switch <- function(full_eigen_data,
                                     sample,
                                     resolution,
                                     manual_switch = T, # set to F to plot data without manual switch of false eigen signs to show where these bins would be
                                     label_size=8
                                     ){
  if (manual_switch == T){
    full_eigen_data_temp = full_eigen_data[,c("eigen_control",paste0("eigen_",sample),paste0("switch_",sample))]
    colnames(full_eigen_data_temp) = c("x","y","switch")
    full_eigen_data_temp = full_eigen_data_temp[which(!is.na(full_eigen_data_temp$switch)),] # to do not consider NA compartments
    
    full_eigen_data_temp$switch = factor(full_eigen_data_temp$switch, levels = c("AB","BA","AA","BB"))
    max_axis_lim = max(full_eigen_data_temp$x,full_eigen_data_temp$y)
    min_axis_lim = min(full_eigen_data_temp$x,full_eigen_data_temp$y)
    
    AB_num = as.character(table(as.factor(full_eigen_data_temp$switch))["AB"])
    BA_num = as.character(table(as.factor(full_eigen_data_temp$switch))["BA"])
    AB_perc = paste0(round(table(as.factor(full_eigen_data_temp$switch))["AB"] / nrow(full_eigen_data_temp) * 100,2),"%")
    BA_perc = paste0(round(table(as.factor(full_eigen_data_temp$switch))["BA"] / nrow(full_eigen_data_temp) * 100,2),"%")
    scc = round(cor(full_eigen_data_temp$x,full_eigen_data_temp$y,method = "spearman"), digits = 2)
    
    p <- ggplot(full_eigen_data_temp) +
      geom_point(aes(y=y, x=x), size = 0.2) +
      xlim(min_axis_lim, max_axis_lim) +
      ylim(min_axis_lim, max_axis_lim) +
      #geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("SCC = ",scc), size = 2) +
      # geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("A to B:\n",AB_num,"(",AB_perc,")")) +
      # geom_text(aes(x=min_axis_lim*0.6,y=max_axis_lim*0.8),label=paste0("B to A:\n",BA_num,"(",BA_perc,")")) +
      #scale_color_manual(values=c('#b2bc9b','#6ca7d2','#a81900','#716c6b')) +
      scale_color_manual(values=c('#808080','#808080','black','black')) +
      labs(x = "Control", y = sample) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            text = element_text(size=label_size),
            axis.text.x = element_text(size = label_size, color = "black"),
            axis.text.y = element_text(size = label_size, color = "black"),
            #axis.ticks.length = unit(0.2, "cm"),
            legend.title = element_blank(),
            legend.position = "none")
      #ggtitle(sample)
  } else {
    full_eigen_data_temp = full_eigen_data[,c("chr","eigen_control",paste0("eigen_",sample),paste0("switch_",sample))]
    full_eigen_data_temp$shape = "round"
    if (resolution == "1000000"){
      full_eigen_data_temp[which(full_eigen_data_temp$chr %in% c("chr19","chr21")),"shape"] = "cross"
    } else if (resolution == "500000"){
      full_eigen_data_temp[which(full_eigen_data_temp$chr == "chr21"),"shape"] = "cross"
    }
    full_eigen_data_temp = full_eigen_data_temp[,c(2,3,4,5)]
    colnames(full_eigen_data_temp) = c("x","y","switch","shape")
    full_eigen_data_temp = full_eigen_data_temp[which(!is.na(full_eigen_data_temp$switch)),] # to do not consider NA compartments
    
    full_eigen_data_temp$switch = factor(full_eigen_data_temp$switch, levels = c("AB","BA","AA","BB"))
    max_axis_lim = max(full_eigen_data_temp$x,full_eigen_data_temp$y)
    min_axis_lim = min(full_eigen_data_temp$x,full_eigen_data_temp$y)
    
    AB_num = as.character(table(as.factor(full_eigen_data_temp$switch))["AB"])
    BA_num = as.character(table(as.factor(full_eigen_data_temp$switch))["BA"])
    AB_perc = paste0(round(table(as.factor(full_eigen_data_temp$switch))["AB"] / nrow(full_eigen_data_temp) * 100,2),"%")
    BA_perc = paste0(round(table(as.factor(full_eigen_data_temp$switch))["BA"] / nrow(full_eigen_data_temp) * 100,2),"%")
    
    p <- ggplot(full_eigen_data_temp) +
      geom_point(aes(y=y, x=x, color=switch, shape=shape)) +
      xlim(min_axis_lim, max_axis_lim) +
      ylim(min_axis_lim, max_axis_lim) +
      geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("A to B:\n",AB_num,"(",AB_perc,")")) +
      geom_text(aes(x=min_axis_lim*0.6,y=max_axis_lim*0.8),label=paste0("B to A:\n",BA_num,"(",BA_perc,")")) +
      scale_color_manual(values=c('#b2bc9b','#6ca7d2','#a81900','#716c6b')) +
      scale_shape_manual(values=c(4,19), guide = F) +
      labs(x = "eigen_control", y = paste0("eigen_",sample)) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            text = element_text(size=label_size),
            axis.text.x = element_text(size = label_size, color = "black"),
            axis.text.y = element_text(size = label_size, color = "black"),
            axis.ticks.length = unit(0.2, "cm"),
            legend.title = element_blank()) +
      ggtitle(sample)
  }
  return(p)
}

p5 <- make_barplot_switch(df_switch_summary,"1000000")
p6 <- make_scatter_plot_switch(full_eigen_data_1000000,"NH4OAc","1000000")
p7 <- make_scatter_plot_switch(full_eigen_data_1000000,"FL","1000000")
p8 <- make_scatter_plot_switch(full_eigen_data_1000000,"RNase","1000000")

p9 <- make_barplot_switch(df_switch_summary,"500000")
p10 <- make_scatter_plot_switch(full_eigen_data_500000,"NH4OAc","500000")
p11 <- make_scatter_plot_switch(full_eigen_data_500000,"FL","500000")
p12 <- make_scatter_plot_switch(full_eigen_data_500000,"RNase","500000")

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/eigen_switch_summary_plot_1000000.png", width = 13, height = 3, units = "in", res = 200)
plot_grid(p5,p7,p6,p8, nrow = 1)
dev.off()

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/eigen_switch_summary_plot_500000.png", width = 13, height = 3, units = "in", res = 200)
plot_grid(p9,p11,p10,p12, nrow = 1)
dev.off()

### Make plots for RNase data without manual switch
p5 <- make_barplot_switch(df_switch_summary,"1000000")
p6 <- make_scatter_plot_switch(full_eigen_data_1000000,"NH4OAc","1000000",T)
p7 <- make_scatter_plot_switch(full_eigen_data_1000000,"FL","1000000",T)
p8 <- make_scatter_plot_switch(full_eigen_data_1000000,"RNase","1000000",F)

p9 <- make_barplot_switch(df_switch_summary,"500000")
p10 <- make_scatter_plot_switch(full_eigen_data_500000,"NH4OAc","500000",T)
p11 <- make_scatter_plot_switch(full_eigen_data_500000,"FL","500000",T)
p12 <- make_scatter_plot_switch(full_eigen_data_500000,"RNase","500000",F)

png("/dataOS/rcalandrelli/phase_separation/compartments/result/eigen_switch_summary_plot_NO_MANUAL_SWITCH.png", width = 13, height = 6, units = "in", res = 200)
plot_grid(p5,p6,p7,p8,p9,p10,p11,p12, nrow = 2)
dev.off()


### Make scatter plots with linear regression line and mark dots outside the region -sd_factor*sd:+sd_factor*sd
# make_scatter_plot_regression <- function(full_eigen_data,
#                                      sample,
#                                      resolution,
#                                      sd_factor = 3, 
#                                      label_size=12
# ){
#   
#   na_switches_compartments = union(union(rownames(full_eigen_data[which(is.na(full_eigen_data$switch_NH4OAc)),]),
#                                    rownames(full_eigen_data[which(is.na(full_eigen_data$switch_FL)),])),
#                                    rownames(full_eigen_data[which(is.na(full_eigen_data$switch_RNase)),]))
#   full_eigen_data = full_eigen_data[-as.numeric(na_switches_compartments),]
#   
#   if (sample == "NH4OAc"){
#     full_eigen_data_temp = full_eigen_data[,c("eigen_control",paste0("eigen_",sample),paste0("switch_",sample))]
#     colnames(full_eigen_data_temp) = c("x","y","switch")
# 
#     ### Linear regression model
#     linearMod <- lm(y ~ x, data = full_eigen_data_temp)
#     x.data = data.frame(x = full_eigen_data_temp$x)
#     y.pred = predict(linearMod, newdata = x.data)
#     my_sd = sd(abs(full_eigen_data_temp$y - y.pred))
# 
#     full_eigen_data_temp$y.pred = y.pred
#     full_eigen_data_temp$position = "within"
#     full_eigen_data_temp$diff = full_eigen_data_temp$y - full_eigen_data_temp$y.pred
#     full_eigen_data_temp[which(abs(full_eigen_data_temp$diff) > sd_factor*my_sd &
#                                  full_eigen_data_temp$diff > 0),"position"] = "above"
#     full_eigen_data_temp[which(abs(full_eigen_data_temp$diff) > sd_factor*my_sd &
#                                  full_eigen_data_temp$diff < 0),"position"] = "below"
# 
#     max_axis_lim = max(full_eigen_data_temp$x, full_eigen_data_temp$y, full_eigen_data_temp$y.pred+sd_factor*my_sd)
#     min_axis_lim = min(full_eigen_data_temp$x,full_eigen_data_temp$y, full_eigen_data_temp$y.pred-sd_factor*my_sd)
#     
#     above_num = as.character(table(as.factor(full_eigen_data_temp$position))["above"])
#     below_num = as.character(table(as.factor(full_eigen_data_temp$position))["below"])
#     above_perc = paste0(round(table(as.factor(full_eigen_data_temp$position))["above"] / nrow(full_eigen_data_temp) * 100,2),"%")
#     below_perc = paste0(round(table(as.factor(full_eigen_data_temp$position))["below"] / nrow(full_eigen_data_temp) * 100,2),"%")
#     
#     p <- ggplot(full_eigen_data_temp) +
#       geom_point(aes(y=y, x=x, color=position), size = 2, alpha = 0.5) +
#       geom_line(aes(x=x, y=y.pred)) +
#       geom_line(aes(x=x, y=y.pred+sd_factor*my_sd), linetype = "dashed") +
#       geom_line(aes(x=x, y=y.pred-sd_factor*my_sd), linetype = "dashed") +
#       xlim(min_axis_lim, max_axis_lim) +
#       ylim(min_axis_lim, max_axis_lim) +
#       geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("<-", sd_factor, "*sd:\n",below_num,"(",below_perc,")")) +
#       geom_text(aes(x=min_axis_lim*0.6,y=max_axis_lim*0.8),label=paste0(">", sd_factor, "*sd:\n",above_num,"(",above_perc,")")) +
#       scale_color_manual(values=c('blue','red','gray')) +
#       labs(x = "eigen_control", y = paste0("eigen_",sample)) + 
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background = element_blank(), axis.line = element_line(colour = "black"),
#             text = element_text(size=label_size),
#             axis.text.x = element_text(size = label_size, color = "black"),
#             axis.text.y = element_text(size = label_size, color = "black"),
#             axis.ticks.length = unit(0.2, "cm"),
#             legend.position = "none") +
#       ggtitle(sample)
#   } else {
#     
#     full_eigen_data_temp = full_eigen_data[,c("eigen_control",paste0("eigen_",sample),paste0("switch_",sample))]
#     colnames(full_eigen_data_temp) = c("x","y","switch")
# 
#     ### Linear regression model
#     linearMod <- lm(y ~ x, data = full_eigen_data_temp)
#     x.data = data.frame(x = full_eigen_data_temp$x)
#     y.pred = predict(linearMod, newdata = x.data)
#     my_sd = sd(abs(full_eigen_data_temp$y - y.pred))
#     
#     full_eigen_data_temp$y.pred = y.pred
#     full_eigen_data_temp$position = "within"
#     full_eigen_data_temp$diff = full_eigen_data_temp$y - full_eigen_data_temp$y.pred
#     full_eigen_data_temp[which(abs(full_eigen_data_temp$diff) > sd_factor*my_sd &
#                                  full_eigen_data_temp$diff > 0),"position"] = "above"
#     full_eigen_data_temp[which(abs(full_eigen_data_temp$diff) > sd_factor*my_sd &
#                                  full_eigen_data_temp$diff < 0),"position"] = "below"
#     
#     ### NH4OAc
#     # full_eigen_data_temp = full_eigen_data[,c("eigen_control",paste0("eigen_",sample),paste0("switch_",sample))]
#     # colnames(full_eigen_data_temp) = c("x","y","switch")
#     # full_eigen_data_temp = full_eigen_data_temp[which(!is.na(full_eigen_data_temp$switch)),] # to do not consider NA compartments
#     
#     ### Linear regression model
#     linearMod <- lm(eigen_NH4OAc ~ eigen_control, data = full_eigen_data)
#     x.data = data.frame(eigen_control = full_eigen_data$eigen_control)
#     y.pred = predict(linearMod, newdata = x.data)
#     my_sd_NH4OAc = sd(abs(full_eigen_data$eigen_NH4OAc - y.pred))
#     
#     full_eigen_data_temp$position_NH4OAc <- "within"
#   
#     full_eigen_data_temp$y.pred_NH4OAc = y.pred
#     full_eigen_data_temp$diff_NH4OAc = full_eigen_data$eigen_NH4OAc - full_eigen_data_temp$y.pred_NH4OA
#     full_eigen_data_temp[which(abs(full_eigen_data_temp$diff_NH4OAc) > sd_factor*my_sd_NH4OAc &
#                                  full_eigen_data_temp$diff_NH4OAc > 0),"position_NH4OAc"] = "above"
#     full_eigen_data_temp[which(abs(full_eigen_data_temp$diff_NH4OAc) > sd_factor*my_sd_NH4OAc &
#                                  full_eigen_data_temp$diff_NH4OAc < 0),"position_NH4OAc"] = "below"
#     
#     max_axis_lim = max(full_eigen_data_temp$x, full_eigen_data_temp$y, full_eigen_data_temp$y.pred+sd_factor*my_sd)
#     min_axis_lim = min(full_eigen_data_temp$x,full_eigen_data_temp$y, full_eigen_data_temp$y.pred-sd_factor*my_sd)
#     
#     above_num = as.character(table(as.factor(full_eigen_data_temp$position))["above"])
#     below_num = as.character(table(as.factor(full_eigen_data_temp$position))["below"])
#     above_perc = paste0(round(table(as.factor(full_eigen_data_temp$position))["above"] / nrow(full_eigen_data_temp) * 100,2),"%")
#     below_perc = paste0(round(table(as.factor(full_eigen_data_temp$position))["below"] / nrow(full_eigen_data_temp) * 100,2),"%")
#     
#     full_eigen_data_temp$position = factor(full_eigen_data_temp$position, c("within", "below", "above"))
#     full_eigen_data_temp$position_NH4OAc = factor(full_eigen_data_temp$position_NH4OAc, c("above", "below", "within"))
#     
#     p <- ggplot(full_eigen_data_temp) +
#       geom_point(aes(y=y, x=x, color=position_NH4OAc), size = 2, alpha = 0.5) +
#       geom_line(aes(x=x, y=y.pred)) +
#       geom_line(aes(x=x, y=y.pred+sd_factor*my_sd), linetype = "dashed") +
#       geom_line(aes(x=x, y=y.pred-sd_factor*my_sd), linetype = "dashed") +
#       xlim(min_axis_lim, max_axis_lim) +
#       ylim(min_axis_lim, max_axis_lim) +
#       geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("<-", sd_factor, "*sd:\n",below_num,"(",below_perc,")")) +
#       geom_text(aes(x=min_axis_lim*0.6,y=max_axis_lim*0.8),label=paste0(">", sd_factor, "*sd:\n",above_num,"(",above_perc,")")) +
#       scale_color_manual(values=c('blue', 'red', 'gray'), breaks = levels(full_eigen_data_temp$position_NH4OAc)) +
#       labs(x = "eigen_control", y = paste0("eigen_",sample)) + 
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background = element_blank(), axis.line = element_line(colour = "black"),
#             text = element_text(size=label_size),
#             axis.text.x = element_text(size = label_size, color = "black"),
#             axis.text.y = element_text(size = label_size, color = "black"),
#             axis.ticks.length = unit(0.2, "cm"),
#             legend.position = "none") +
#       ggtitle(sample)
#   }
#   return(p)
# }


# make_scatter_plot_outliers <- function(full_eigen_data,
#                                sample,
#                                resolution,
#                                sd_factor = 2, 
#                                label_size=12
# ){
#   
#   na_switches_compartments = union(union(rownames(full_eigen_data[which(is.na(full_eigen_data$switch_NH4OAc)),]),
#                                          rownames(full_eigen_data[which(is.na(full_eigen_data$switch_FL)),])),
#                                    rownames(full_eigen_data[which(is.na(full_eigen_data$switch_RNase)),]))
#   full_eigen_data = full_eigen_data[-as.numeric(na_switches_compartments),]
#   
#   # Process NH4OAc to color dots differently
#   full_eigen_data$diff_NH4OAc = full_eigen_data$eigen_NH4OAc - full_eigen_data$eigen_control
#   sd_NH4OAc = sd(abs(full_eigen_data$diff_NH4OAc))
#   
#   full_eigen_data$position_NH4OAc = "within"
#   full_eigen_data[which(abs(full_eigen_data$diff_NH4OAc) > sd_factor*sd_NH4OAc &
#                           full_eigen_data$diff_NH4OAc > 0),"position_NH4OAc"] = "above"
#   full_eigen_data[which(abs(full_eigen_data$diff_NH4OAc) > sd_factor*sd_NH4OAc &
#                                full_eigen_data$diff_NH4OAc < 0),"position_NH4OAc"] = "below"
#   
#   # Process selected sample
#   colnames(full_eigen_data)[which(colnames(full_eigen_data) == "eigen_control")] = "x"
#   colnames(full_eigen_data)[which(colnames(full_eigen_data) == paste0("eigen_",sample))] = "y"
#   
#   full_eigen_data$diff = full_eigen_data$y - full_eigen_data$x
#   my_sd = sd(abs(full_eigen_data$diff))
#   
#   max_axis_lim = max(full_eigen_data$x, full_eigen_data$y, full_eigen_data$x+sd_factor*my_sd)
#   min_axis_lim = min(full_eigen_data$x, full_eigen_data$y, full_eigen_data$x-sd_factor*my_sd)
#   
#   full_eigen_data$position = "within"
#   full_eigen_data[which(abs(full_eigen_data$diff) > sd_factor*my_sd &
#                           full_eigen_data$diff > 0),"position"] = "above"
#   full_eigen_data[which(abs(full_eigen_data$diff) > sd_factor*my_sd &
#                           full_eigen_data$diff < 0),"position"] = "below"
#   
#   above_num = as.character(table(as.factor(full_eigen_data$position))["above"])
#   below_num = as.character(table(as.factor(full_eigen_data$position))["below"])
#   above_perc = paste0(round(table(as.factor(full_eigen_data$position))["above"] / nrow(full_eigen_data) * 100,2),"%")
#   below_perc = paste0(round(table(as.factor(full_eigen_data$position))["below"] / nrow(full_eigen_data) * 100,2),"%")
#   
#   p <- ggplot(full_eigen_data) +
#     geom_point(aes(y=y, x=x, color=position_NH4OAc), size = 2, alpha = 0.5) +
#     geom_line(aes(x=x, y=x)) +
#     geom_line(aes(x=x, y=x+sd_factor*my_sd), linetype = "dashed") +
#     geom_line(aes(x=x, y=x-sd_factor*my_sd), linetype = "dashed") +
#     # xlim(min_axis_lim, max_axis_lim) +
#     # ylim(min_axis_lim, max_axis_lim) +
#     geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("<-", sd_factor, "*sd:\n",below_num,"(",below_perc,")")) +
#     geom_text(aes(x=min_axis_lim*0.6,y=max_axis_lim*0.8),label=paste0(">", sd_factor, "*sd:\n",above_num,"(",above_perc,")")) +
#     scale_color_manual(values=c('blue','red','gray')) +
#     labs(x = "eigen_control", y = paste0("eigen_",sample)) + 
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"),
#           text = element_text(size=label_size),
#           axis.text.x = element_text(size = label_size, color = "black"),
#           axis.text.y = element_text(size = label_size, color = "black"),
#           axis.ticks.length = unit(0.2, "cm"),
#           legend.position = "none") +
#     ggtitle(sample)
#   return(p)
# }


make_scatter_plot_outliers <- function(full_eigen_data,
                                       sample,
                                       resolution,
                                       multiple_factor = 2, 
                                       label_size=16){
  
  na_switches_compartments = union(union(rownames(full_eigen_data[which(is.na(full_eigen_data$switch_NH4OAc)),]),
                                         rownames(full_eigen_data[which(is.na(full_eigen_data$switch_FL)),])),
                                   rownames(full_eigen_data[which(is.na(full_eigen_data$switch_RNase)),]))
  full_eigen_data = full_eigen_data[-as.numeric(na_switches_compartments),]
  
  # Process NH4OAc to color dots differently
  full_eigen_data$position_NH4OAc = "within"
  full_eigen_data[which((full_eigen_data$eigen_NH4OAc > multiple_factor * full_eigen_data$eigen_control &
                           !(full_eigen_data$eigen_NH4OAc < 0 &
                           full_eigen_data$eigen_control < 0)) |
                           (full_eigen_data$eigen_NH4OAc > 1 / multiple_factor * full_eigen_data$eigen_control &
                             full_eigen_data$eigen_NH4OAc < 0 &
                             full_eigen_data$eigen_control < 0)),"position_NH4OAc"] = "above"
  
  full_eigen_data[which((full_eigen_data$eigen_NH4OAc < 1 / multiple_factor * full_eigen_data$eigen_control &
                           !(full_eigen_data$eigen_NH4OAc < 0 &
                               full_eigen_data$eigen_control < 0)) |
                          (full_eigen_data$eigen_NH4OAc < multiple_factor * full_eigen_data$eigen_control &
                             full_eigen_data$eigen_NH4OAc < 0 &
                             full_eigen_data$eigen_control < 0)),"position_NH4OAc"] = "below"  
  
  # Process selected sample
  colnames(full_eigen_data)[which(colnames(full_eigen_data) == "eigen_control")] = "x"
  colnames(full_eigen_data)[which(colnames(full_eigen_data) == paste0("eigen_",sample))] = "y"
  
  # max_axis_lim = max(full_eigen_data$x, full_eigen_data$y, full_eigen_data$x*multiple_factor)
  # min_axis_lim = min(full_eigen_data$x, full_eigen_data$y, -full_eigen_data$x*multiple_factor)
  
  full_eigen_data$position = "within"
  full_eigen_data[which((full_eigen_data$y > multiple_factor * full_eigen_data$x &
                           !(full_eigen_data$y < 0 &
                               full_eigen_data$x < 0)) |
                          (full_eigen_data$y > 1 / multiple_factor * full_eigen_data$x &
                             full_eigen_data$y < 0 &
                             full_eigen_data$x < 0)),"position"] = "above"
  
  full_eigen_data[which((full_eigen_data$y < 1 / multiple_factor * full_eigen_data$x &
                           !(full_eigen_data$y < 0 &
                               full_eigen_data$x < 0)) |
                          (full_eigen_data$y < multiple_factor * full_eigen_data$x &
                             full_eigen_data$y < 0 &
                             full_eigen_data$x < 0)),"position"] = "below"  
  
  above_num = as.character(table(as.factor(full_eigen_data$position))["above"])
  below_num = as.character(table(as.factor(full_eigen_data$position))["below"])
  above_perc = paste0(round(table(as.factor(full_eigen_data$position))["above"] / nrow(full_eigen_data) * 100,2),"%")
  below_perc = paste0(round(table(as.factor(full_eigen_data$position))["below"] / nrow(full_eigen_data) * 100,2),"%")
  
  annotations <- data.frame(
    xpos = c(Inf,-Inf),
    ypos =  c(-Inf,Inf),
    annotateText = c(paste0(below_num," (",below_perc,")"),
                     paste0(above_num," (",above_perc,")")),
    hjustvar = c(1.5,-0.5),
    vjustvar = c(-2,2))
  
  p <- ggplot(full_eigen_data) +
    geom_point(aes(y=y, x=x, color=position_NH4OAc), size = 2, alpha = 0.5) +
    geom_line(aes(x=x, y=x)) +
    geom_line(aes(x=x, y=x*multiple_factor), linetype = "dashed") +
    geom_line(aes(x=x, y=x/multiple_factor), linetype = "dashed") +
    # xlim(min_axis_lim, max_axis_lim) +
    # ylim(min_axis_lim, max_axis_lim) +
    # geom_text(aes(x=max_axis_lim*0.6,y=min_axis_lim*0.8),label=paste0("<-", sd_factor, "*sd:\n",below_num,"(",below_perc,")")) +
    # geom_text(aes(x=min_axis_lim*0.6,y=max_axis_lim*0.8),label=paste0(">", sd_factor, "*sd:\n",above_num,"(",above_perc,")")) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=5) +
    scale_color_manual(values=c('blue','red','gray')) +
    labs(x = "eigen_control", y = paste0("eigen_",sample)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size = label_size, color = "black"),
          axis.text.y = element_text(size = label_size, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none") +
    ggtitle(sample)
  return(p)
}




p1 <- make_scatter_plot_outliers(full_eigen_data_1000000,"NH4OAc","1000000")
p2 <- make_scatter_plot_outliers(full_eigen_data_1000000,"FL","1000000")
p3 <- make_scatter_plot_outliers(full_eigen_data_1000000,"RNase","1000000")
p4 <- make_scatter_plot_outliers(full_eigen_data_500000,"NH4OAc","500000")
p5 <- make_scatter_plot_outliers(full_eigen_data_500000,"FL","500000")
p6 <- make_scatter_plot_outliers(full_eigen_data_500000,"RNase","500000")

png("/dataOS/rcalandrelli/phase_separation/compartments/result/scatter_plots_outliers.png", width = 12, height = 8, units = "in", res = 200)
plot_grid(p1,p2,p3,p4,p5,p6, ncol = 3)
dev.off()


### Venn diagram to check if switched bins are shared between samples
library(VennDiagram)
make_venn_diagram <- function(full_eigen_data,
                              resolution,
                              switch_type){
  
  NH4OAc <- which(full_eigen_data$switch_NH4OAc == switch_type)
  FL <- which(full_eigen_data$switch_FL == switch_type)
  RNase <- which(full_eigen_data$switch_RNase == switch_type)
  
  library(RColorBrewer)
  myCol <- brewer.pal(3, "Pastel2")
  
  # Chart
  venn.diagram(
    x = list(NH4OAc, FL, RNase),
    category.names = c(paste0("NH4OAc (",length(NH4OAc),")"),
                       paste0("FL (",length(FL),")"),
                       paste0("RNase (",length(RNase),")")),
    filename = paste0('/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/venn_diagrams/venn_diagram_',switch_type,"_",resolution,'.png'),
    output=T,
    
    # Main title
    main = switch_type,
    main.fontfamily = "sans",
    main.fontface = "bold",
    main.cex = 1,
    
    # Output features
    imagetype="png" ,
    height = 650 , 
    width = 650 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
}

### 1 Mb

make_venn_diagram(full_eigen_data_1000000,"1000000","AB")
make_venn_diagram(full_eigen_data_1000000,"1000000","BA")
make_venn_diagram(full_eigen_data_1000000,"1000000","AA")
make_venn_diagram(full_eigen_data_1000000,"1000000","BB")

make_venn_diagram(full_eigen_data_500000,"500000","AB")
make_venn_diagram(full_eigen_data_500000,"500000","BA")
make_venn_diagram(full_eigen_data_500000,"500000","AA")
make_venn_diagram(full_eigen_data_500000,"500000","BB")


switches = c("AB","BA","AA","BB")
resolutions = c("1000000","500000")
figures_list = list()
k = 1
for (i in resolutions){
  for (j in switches){
    figures_list[[k]] = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/venn_diagrams/venn_diagram_",j,"_",i,".png"))
    k = k + 1
  }
}

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/venn_diagrams_1000000.png", width = 8, height = 2.2, units = "in", res = 200)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]), ncol=4, left = "Res: 1000000")
dev.off()

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/venn_diagrams_500000.png", width = 8, height = 2.2, units = "in", res = 200)
grid.arrange(rasterGrob(figures_list[[5]]),
             rasterGrob(figures_list[[6]]),
             rasterGrob(figures_list[[7]]),
             rasterGrob(figures_list[[8]]), ncol=4, left = "Res: 500000")
dev.off()

######### Upset plot
library(UpSetR)

make_upset_switch <- function(full_eigen_data,
                              switch_type,
                              resolution){
  temp_upset = full_eigen_data
  temp_upset = temp_upset[which(temp_upset$switch_NH4OAc %in% c("AA","BB","AB","BA")),]
  temp_upset = temp_upset[which(temp_upset$switch_FL %in% c("AA","BB","AB","BA")),]
  temp_upset = temp_upset[which(temp_upset$switch_RNase %in% c("AA","BB","AB","BA")),]
  
  temp_upset$switch_NH4OAc <- ifelse(temp_upset$switch_NH4OAc == switch_type, 1, 0)
  temp_upset$switch_FL <- ifelse(temp_upset$switch_FL == switch_type, 1, 0)  
  temp_upset$switch_RNase <- ifelse(temp_upset$switch_RNase == switch_type, 1, 0)  
  
  colnames(temp_upset)[7:9] =  c("NH4OAc","FL","RNase")
  p<-upset(temp_upset, 
           sets = c("RNase","FL","NH4OAc"),
           sets.x.label = paste0(switch_type, " switches"),
           mainbar.y.label = paste0(switch_type, " intersection size"),
           order.by = "freq",
           keep.order = T,
           text.scale = c(3, 2.5, 2, 1.5, 2, 2),
           point.size = 3, line.size = 1)
  png(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/upset_plots/",resolution,"_",switch_type,".png"), height = 5.5, width = 6, res = 200, units = "in")
  print(p)
  dev.off()
}

make_upset_switch(full_eigen_data_1000000,"AB","1000000")
make_upset_switch(full_eigen_data_1000000,"BA","1000000")
make_upset_switch(full_eigen_data_1000000,"AA","1000000")
make_upset_switch(full_eigen_data_1000000,"BB","1000000")

make_upset_switch(full_eigen_data_500000,"AB","500000")
make_upset_switch(full_eigen_data_500000,"BA","500000")
make_upset_switch(full_eigen_data_500000,"AA","500000")
make_upset_switch(full_eigen_data_500000,"BB","500000")

switches = c("AB","BA","AA","BB")
resolutions = c("1000000","500000")
figures_list = list()
k = 1
for (i in resolutions){
  for (j in switches){
    figures_list[[k]] = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/upset_plots/",i,"_",j,".png"))
    k = k + 1
  }
}

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/upset_plots_1000000.png", width = 5, height = 4, units = "in", res = 400)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]),
             ncol=2, left = "Res: 1000000")
dev.off()

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/upset_plots_500000.png", width = 5, height = 4, units = "in", res = 400)
grid.arrange(rasterGrob(figures_list[[5]]),
             rasterGrob(figures_list[[6]]),
             rasterGrob(figures_list[[7]]),
             rasterGrob(figures_list[[8]]), ncol=2, left = "Res: 500000")
dev.off()


# Plot again chr21 with manual switches

my_chr = "chr21"
my_resolution = "500000"
my_plot_type = "pcc"

g1 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                              plot_type = my_plot_type,
                              sample = "H1_control_merged",
                              resolution = my_resolution,
                              chr = my_chr,
                              fig.title = paste0("Control_",my_chr))

g2 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                              plot_type = my_plot_type,
                              sample = "H1_NH4OAc_merged",
                              resolution = my_resolution,
                              chr = my_chr,
                              fig.title = paste0("NH4OAc_",my_chr))

g3 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                              plot_type = my_plot_type,
                              sample = "H1_FL_merged",
                              resolution = my_resolution,
                              chr = my_chr,
                              fig.title = paste0("FL_",my_chr))

g4 <- make_heatmap_eigen_plot(directory = "/dataOS/rcalandrelli/phase_separation/HiC/",
                              plot_type = my_plot_type,
                              sample = "H1_RNase_merged",
                              resolution = my_resolution,
                              chr = my_chr,
                              fig.title = paste0("Rnase_",my_chr),
                              force_manual_switch = T)

png(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/",my_chr,"_",my_plot_type,"_heatmap_eigen_",my_resolution,"_MANUAL_SWITCH.png"), width = 14, height = 5, units = "in", res = 200)
grid.arrange(g1,g2,g3,g4,nrow=1)
dev.off()



##### Compaction analysis
chrSize = read.table("/dataOS/chenweizhong/Project____InterChrom/Data/4DN_files_reference/4DNFI823LSII.hg38.mainonly.chrom.sizes")

calculate_compaction_track <- function(directory, # with trailing slash in the end
                                       sample,
                                       bin_size,
                                       local_limit){
  if (!dir.exists(paste0(directory, "compartments/DLR_analysis/", bin_size))){
    dir.create(paste0(directory, "compartments/DLR_analysis/", bin_size))
  }
  
  if (!dir.exists(paste0(directory, "compartments/DLR_analysis/", bin_size, "/", sample))){
    dir.create(paste0(directory, "compartments/DLR_analysis/", bin_size, "/", sample))
  }
  
  output_list = list()
  # counter = 0
  local_limit_bin = local_limit/bin_size
  for (i in hg38_chromosomes){
    chr_size_bin = round((chrSize[which(chrSize$V1==i),2])/bin_size)
    # Load sparse matrix (check if file size is non-zero first!!!)
    if (file.size(paste0(directory, "HiC_contact_matrices/KR/", sample, "/", bin_size, "/",i,"_",i,"_",bin_size,".txt")) != 0){
      print(paste0("Performing compaction analysis for ",i))
      input_matrix_sparse = data.frame(fread(paste0(directory, "HiC_contact_matrices/KR/", sample, "/", bin_size, "/",i,"_",i,"_",bin_size,".txt"), stringsAsFactors = F))
      input_matrix_sparse[which(is.na(input_matrix_sparse[,3])),3] = 0 # removing NaN values
      input_matrix_sparse_bin = input_matrix_sparse
      input_matrix_sparse_bin[,1] = input_matrix_sparse_bin[,1] / bin_size + 1
      input_matrix_sparse_bin[,2] = input_matrix_sparse_bin[,2] / bin_size + 1
      out_track_chr = matrix(NA,nrow=chr_size_bin,ncol=7)
      colnames(out_track_chr) = c("chr","bin","local_up","local_down","distal_up","distal_down","DLR")
      out_track_chr = data.frame(out_track_chr)
      out_track_chr$chr = i
      out_track_chr$bin = 1:nrow(out_track_chr)
      for (j in 1:nrow(out_track_chr)){
        temp_bin_mat = input_matrix_sparse_bin[which(input_matrix_sparse_bin$V1 == j |
                                                       input_matrix_sparse_bin$V2 == j),]
        if (nrow(temp_bin_mat) > 0){
          upstream_local =  sum(temp_bin_mat[which(temp_bin_mat$V2 == j &
                                                 temp_bin_mat$V1 %in% seq(j-local_limit_bin,j-1)),3])
          downstream_local =  sum(temp_bin_mat[which(temp_bin_mat$V1 == j &
                                                 temp_bin_mat$V2 %in% seq(j+1,j+local_limit_bin)),3])
          upstream_distal =  sum(temp_bin_mat[which(temp_bin_mat$V2 == j &
                                                 temp_bin_mat$V1 %in% seq(1,j-local_limit_bin-1)),3])
          downstream_distal =  sum(temp_bin_mat[which(temp_bin_mat$V1 == j &
                                                   temp_bin_mat$V2 %in% seq(j+local_limit_bin+1,nrow(out_track_chr))),3])
          
          out_track_chr[j,"local_up"] = upstream_local
          out_track_chr[j,"local_down"] = downstream_local
          out_track_chr[j,"distal_up"] = upstream_distal
          out_track_chr[j,"distal_down"] = downstream_distal
          out_track_chr[j,"DLR"] = log2((upstream_distal+downstream_distal)/(upstream_local+downstream_local))
        }
      }

      write.table(out_track_chr, paste0(directory, "compartments/DLR_analysis/", bin_size, "/", sample, "/",i,".txt"), row.names = F, col.names = T, sep = '\t', quote = F)
      output_list[[i]] = out_track_chr
    }
  }
  return(output_list)
}

# 500 kb resolution
DLR_control_500000 = calculate_compaction_track(directory="/dataOS/rcalandrelli/phase_separation/HiC/", 
                                           sample="H1_control_merged",
                                           bin_size=500000,
                                           local_limit=3000000)
DLR_NH4OAc_500000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                           "H1_NH4OAc_merged",
                                           500000,
                                           3000000)
DLR_FL_500000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                           "H1_FL_merged",
                                           500000,
                                           3000000)
DLR_RNase_500000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                           "H1_RNase_merged",
                                           500000,
                                           3000000)


# DLR analysis
df_DRL_500000 = data.frame("Control" = do.call("rbind",DLR_control_500000)[,"DLR"],
                     "FL" = do.call("rbind",DLR_FL_500000)[,"DLR"],
                     "NH4OAc" = do.call("rbind",DLR_NH4OAc_500000)[,"DLR"],
                     "RNase" = do.call("rbind",DLR_RNase_500000)[,"DLR"])

df1 = melt(df_DRL_500000)

mean(df1[which(df1$variable=="Control" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])
mean(df1[which(df1$variable=="NH4OAc" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])
mean(df1[which(df1$variable=="FL" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])
mean(df1[which(df1$variable=="RNase" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])

p1<-ggplot(df1, aes(y=value, x=variable)) +
  geom_boxplot() +
  labs(x = "", y = "DLR") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Res: 500 kb")

wilcox.test(df1[which(df1$variable=="Control"),"value"], df1[which(df1$variable=="NH4OAc"),"value"])
wilcox.test(df1[which(df1$variable=="Control"),"value"], df1[which(df1$variable=="FL"),"value"])
wilcox.test(df1[which(df1$variable=="Control"),"value"], df1[which(df1$variable=="RNase"),"value"])

# Local and distal interactions
df_local = melt(data.frame("Control" = rowSums(do.call("rbind",DLR_control_500000)[,c("local_up","local_down")]),
                     "FL" = rowSums(do.call("rbind",DLR_FL_500000)[,c("local_up","local_down")]),
                     "NH4OAc" = rowSums(do.call("rbind",DLR_NH4OAc_500000)[,c("local_up","local_down")]),
                     "RNase" = rowSums(do.call("rbind",DLR_RNase_500000)[,c("local_up","local_down")])))
df_local$label = "Local"

wilcox.test(df_local[which(df_local$variable=="Control"),"value"], df_local[which(df_local$variable=="NH4OAc"),"value"])
wilcox.test(df_local[which(df_local$variable=="Control"),"value"], df_local[which(df_local$variable=="FL"),"value"])
wilcox.test(df_local[which(df_local$variable=="Control"),"value"], df_local[which(df_local$variable=="RNase"),"value"])

df_distal = melt(data.frame("Control" = rowSums(do.call("rbind",DLR_control_500000)[,c("distal_up","distal_down")]),
                           "FL" = rowSums(do.call("rbind",DLR_FL_500000)[,c("distal_up","distal_down")]),
                           "NH4OAc" = rowSums(do.call("rbind",DLR_NH4OAc_500000)[,c("distal_up","distal_down")]),
                           "RNase" = rowSums(do.call("rbind",DLR_RNase_500000)[,c("distal_up","distal_down")])))
df_distal$label = "Distal"

wilcox.test(df_distal[which(df_distal$variable=="Control"),"value"], df_distal[which(df_distal$variable=="NH4OAc"),"value"])
wilcox.test(df_distal[which(df_distal$variable=="Control"),"value"], df_distal[which(df_distal$variable=="FL"),"value"])
wilcox.test(df_distal[which(df_distal$variable=="Control"),"value"], df_distal[which(df_distal$variable=="RNase"),"value"])

df2 = rbind(df_local,df_distal)

p2<-ggplot(df2, aes(y=value, x=variable, fill=label)) +
  geom_boxplot() +
  labs(x = "", y = "Number of interactions") + 
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Res: 500 kb")

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/DLR_500000_boxplots.png", width = 12, height = 6, units = "in", res = 200)
plot_grid(p1,p2)
dev.off()


# 50 kb resolution
DLR_control_50000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                                "H1_control_merged",
                                                50000,
                                                3000000)
DLR_NH4OAc_50000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                               "H1_NH4OAc_merged",
                                               50000,
                                               3000000)
DLR_FL_50000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                           "H1_FL_merged",
                                           50000,
                                           3000000)
DLR_RNase_50000 = calculate_compaction_track("/dataOS/rcalandrelli/phase_separation/HiC/", 
                                              "H1_RNase_merged",
                                              50000,
                                              3000000)

DLR_control_50000 = list()
DLR_FL_50000= list()
DLR_NH4OAc_50000 = list()
DLR_RNase_50000 = list()

for (i in hg38_chromosomes){
  temp = read.table(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/DLR_analysis/50000/H1_control_merged/",i,".txt"), header = T)
  DLR_control_50000[[i]] = temp
  
  temp = read.table(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/DLR_analysis/50000/H1_FL_merged/",i,".txt"), header = T)
  DLR_FL_50000[[i]] = temp
  
  temp = read.table(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/DLR_analysis/50000/H1_NH4OAc_merged/",i,".txt"), header = T)
  DLR_NH4OAc_50000[[i]] = temp
  
  temp = read.table(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/DLR_analysis/50000/H1_RNase_merged/",i,".txt"), header = T)
  DLR_RNase_50000[[i]] = temp
}


df_DRL_50000 = data.frame("Control" = do.call("rbind",DLR_control_50000)[,"DLR"],
                          "FL" = do.call("rbind",DLR_FL_50000)[,"DLR"],
                           "NH4OAc" = do.call("rbind",DLR_NH4OAc_50000)[,"DLR"],
                           "RNase" = do.call("rbind",DLR_RNase_50000)[,"DLR"])

df1 = melt(df_DRL_50000)

mean(df1[which(df1$variable=="Control" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])
mean(df1[which(df1$variable=="NH4OAc" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])
mean(df1[which(df1$variable=="FL" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])
mean(df1[which(df1$variable=="RNase" & !is.na(df1$value) & !is.infinite(df1$value)), "value"])

p1<-ggplot(df1, aes(y=value, x=variable)) +
  geom_boxplot() +
  labs(x = "", y = "DLR") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Res: 50 kb")

wilcox.test(df[which(df$variable=="Control"),"value"], df[which(df$variable=="NH4OAc"),"value"])
wilcox.test(df[which(df$variable=="Control"),"value"], df[which(df$variable=="FL"),"value"])
wilcox.test(df[which(df$variable=="Control"),"value"], df[which(df$variable=="RNase"),"value"])

# Local and distal interactions
df_local = melt(data.frame("Control" = rowSums(do.call("rbind",DLR_control_50000)[,c("local_up","local_down")]),
                           "NH4OAc" = rowSums(do.call("rbind",DLR_NH4OAc_50000)[,c("local_up","local_down")]),
                           "FL" = rowSums(do.call("rbind",DLR_FL_50000)[,c("local_up","local_down")]),
                           "RNase" = rowSums(do.call("rbind",DLR_RNase_50000)[,c("local_up","local_down")])))
df_local$label = "Local"

wilcox.test(df_local[which(df_local$variable=="Control"),"value"], df_local[which(df_local$variable=="NH4OAc"),"value"])
wilcox.test(df_local[which(df_local$variable=="Control"),"value"], df_local[which(df_local$variable=="FL"),"value"])
wilcox.test(df_local[which(df_local$variable=="Control"),"value"], df_local[which(df_local$variable=="RNase"),"value"])

df_distal = melt(data.frame("Control" = rowSums(do.call("rbind",DLR_control_50000)[,c("distal_up","distal_down")]),
                            "NH4OAc" = rowSums(do.call("rbind",DLR_NH4OAc_50000)[,c("distal_up","distal_down")]),
                            "FL" = rowSums(do.call("rbind",DLR_FL_50000)[,c("distal_up","distal_down")]),
                            "RNase" = rowSums(do.call("rbind",DLR_RNase_50000)[,c("distal_up","distal_down")])))
df_distal$label = "Distal"

wilcox.test(df_distal[which(df_distal$variable=="Control"),"value"], df_distal[which(df_distal$variable=="NH4OAc"),"value"])
wilcox.test(df_distal[which(df_distal$variable=="Control"),"value"], df_distal[which(df_distal$variable=="FL"),"value"])
wilcox.test(df_distal[which(df_distal$variable=="Control"),"value"], df_distal[which(df_distal$variable=="RNase"),"value"])

df2 = rbind(df_local,df_distal)

p2<-ggplot(df2, aes(y=value, x=variable, fill=label)) +
  geom_boxplot() +
  labs(x = "", y = "Number of interactions") + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Res: 50 kb")

png("/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/DLR_50000_boxplots.png", width = 12, height = 6, units = "in", res = 200)
plot_grid(p1,p2)
dev.off()




################################################ Repetitive element analysis

hg38_repeatMasker = data.frame(fread("/dataOS/wenxingzhao/database/Repeat_masker_human/hg38.fa.out.bed"))

temp = names(table(as.factor(hg38_repeatMasker$V4)))

grep(pattern = "B1", temp)

"LINE/L1"
"LINE/L2"
"SINE/Alu"
"LTR/ERV1"

compartment_repeat_enrichment <- function(eigen_data,
                                          repeat_label,
                                          sample_name,
                                          resolution_flanking_region=10000,
                                          flanking_region_size=1000000){
  
  ### Repeat data
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == repeat_label),]
  Gr_repeat_data = GRanges(
    seqnames = Rle(repeat_data[,1]),
    ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
    strand = Rle(strand(repeat_data[,6])))
  
  ### Eigen data
  eigen_switch_index = c()
  for (i in 1:(nrow(eigen_data)-1)){
    if (eigen_data[i,3] * eigen_data[i+1,3] < 0){
      eigen_switch_index = c(eigen_switch_index, i)
    }
  }
  
  eigen_switch_regions = eigen_data[eigen_switch_index,1:2]
  eigen_switch_regions = eigen_switch_regions[which(eigen_switch_regions$coord - flanking_region_size > 0),]
  eigen_switch_regions$coord = eigen_switch_regions$coord - flanking_region_size + 500000
  colnames(eigen_switch_regions)[2] = "start"
  eigen_switch_regions$end = eigen_switch_regions$start + 2*flanking_region_size
  
  
  repeat_density_track_switch_region <- function(x){
    Gr_switch_region = GRanges(
      seqnames = Rle(eigen_switch_regions[x,1]),
      ranges = IRanges(as.numeric(eigen_switch_regions[x,2]), end = as.numeric(eigen_switch_regions[x,3]), names = rownames(eigen_switch_regions[x,])),
      strand = Rle(strand("*")))
    
    Gr_switch_region_tile = tile(Gr_switch_region, 2*flanking_region_size/resolution_flanking_region)[[1]]
    temp = as.numeric(countOverlaps(Gr_switch_region_tile, Gr_repeat_data, ignore.strand = T))
    return(temp)
  }
  
  out_density_list = lapply(1:nrow(eigen_switch_regions), repeat_density_track_switch_region)
  df_heatmap = do.call("rbind", out_density_list)
  
  df_heatmap = df_heatmap[sort_order_heatmap[,"index"],]
  
  png(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/heatmap_",sample_name,"_",gsub(".*/","",repeat_label),".png"), width = 3, height = 5, res = 200, units = "in")
  pheatmap(df_heatmap,
           color = colorRampPalette(c("blue", "yellow"))(100),
           cluster_rows=FALSE,
           cluster_cols=FALSE)
  dev.off()
  
}

# To sort by decreasing B1 values in the first half of the heatmap
eigen_data=eigen_control_500000
repeat_label="SINE/Alu"
sample_name="H1_control"
resolution_flanking_region = 10000
flanking_region_size = 1000000

sort_order_heatmap = rowSums(df_heatmap[,1:(ncol(df_heatmap)/2)])
sort_order_heatmap = cbind(sort_order_heatmap, seq(1,length(sort_order_heatmap)))
colnames(sort_order_heatmap) = c("value","index")
sort_order_heatmap = sort_order_heatmap[order(sort_order_heatmap[,"value"], decreasing = T),]



dev.off()
compartment_repeat_enrichment(eigen_data = eigen_control_500000,
                              repeat_label = "LINE/L1",
                              sample_name = "H1_control",
                              resolution_flanking_region = 10000,
                              flanking_region_size = 1000000)

compartment_repeat_enrichment(eigen_data=eigen_control_500000,
                              repeat_label="SINE/Alu",
                              sample_name="H1_control")

compartment_repeat_enrichment(eigen_data=eigen_control_500000,
                              repeat_label="LINE/L2",
                              sample_name="H1_control")

compartment_repeat_enrichment(eigen_data=eigen_control_500000,
                              repeat_label="LTR/ERV1",
                              sample_name="H1_control")



######## Plotting tracks

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


plot_compartment_repetitive_tracks <- function(eigen_file,
                                               sample_name,
                                               chromosomes,
                                               resolution=500000,
                                               tick_dist=50000000,
                                               minor_tick_dist=10000000
)
{
  
  ### Load compartment data
  compartment = read.table(eigen_file, header = T)
  compartment$eigen_pos = compartment$eigen
  compartment$eigen_neg = compartment$eigen
  
  compartment[,4][which(compartment[,4]<0)] = NA
  compartment[,5][which(compartment[,5]>=0)] = NA

  colnames(compartment) = c("chr", "coord", "eigen", "eigen_pos", "eigen_neg")
  compartment = compartment[which(compartment$chr != "chrY"),]
  
  
  ##### Repetitive elements
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == "SINE/Alu"),]
  Gr_Alu = GRanges(
    seqnames = Rle(repeat_data[,1]),
    ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
    strand = Rle(strand(repeat_data[,6])))
  
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == "LINE/L1"),]
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
  cov_data_repetitive$log2_Alu_L1 = log2(cov_data_repetitive$Alu / cov_data_repetitive$L1)
  
  cov_data_repetitive$log2_Alu_L1_pos = cov_data_repetitive$log2_Alu_L1
  cov_data_repetitive[which(cov_data_repetitive$log2_Alu_L1_pos < 0),"log2_Alu_L1_pos"] = NA
  
  cov_data_repetitive$log2_Alu_L1_neg = cov_data_repetitive$log2_Alu_L1
  cov_data_repetitive[which(cov_data_repetitive$log2_Alu_L1_neg > 0),"log2_Alu_L1_neg"] = NA
  
  cor_index = which(!is.na(compartment$eigen) & !is.na(cov_data_repetitive$log2_Alu_L1) & !is.infinite(cov_data_repetitive$log2_Alu_L1), arr.ind = T)
  out_cor = cor(compartment$eigen[cor_index], cov_data_repetitive$log2_Alu_L1[cor_index], method = "spearman")
  
  
  for (my_chr in chromosomes){
    
      print(paste0(sample_name, " - ", my_chr))
      start_coord = 0
      end_coord = hg38_lengths[my_chr]
      detail.region <- toGRanges(data.frame(my_chr, start_coord, end_coord))
      
      y_text = 0.1 # track label height
      
      line_width = 1 # width of the line track
      label_cex = 0.6
      
      # Alu area
      r1_Alu = 0.95
      r0_Alu = r1_Alu - 0.15
      
      # L1 area
      r1_L1 = r0_Alu - 0.1
      r0_L1 = r1_L1 - 0.15
      
      # log2 repetitive element area
      r1_repetitive = r0_L1 - 0.1
      r0_repetitive = r1_repetitive - 0.15
      temp = cov_data_repetitive[which(cov_data_repetitive$seqnames==my_chr),]
      ymin_repetitive = round(min(temp[which(!is.na(temp$log2_Alu_L1) & !is.infinite(temp$log2_Alu_L1)),"log2_Alu_L1"]),1)
      ymax_repetitive = round(max(temp[which(!is.na(temp$log2_Alu_L1) & !is.infinite(temp$log2_Alu_L1)),"log2_Alu_L1"]),1)
      
      y_repetitive = min(ymax_repetitive,-ymin_repetitive)
      
      # Compartment area
      r1_comp = r0_repetitive - 0.1
      r0_comp = r1_comp - 0.15
      ymin_comp = round(min(compartment[which(compartment$chr==my_chr),"eigen_neg"][which(!is.na(compartment[which(compartment$chr==my_chr),"eigen_neg"]))]),1)
      ymax_comp = round(max(compartment[which(compartment$chr==my_chr),"eigen_pos"][which(!is.na(compartment[which(compartment$chr==my_chr),"eigen_pos"]))]),1)
      
      
      pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/",sample_name,"_",my_chr,"_tracks.pdf"), width = 4, height = 2.5)
      plot_params <- getDefaultPlotParams(plot.type=1)
      plot_params$data1inmargin <- 5
      plot_params$ideogramheight <- 10
      plot_params$bottommargin <- 30
      plot_params$topmargin <- 0
      plot_params$leftmargin <- 0.12
      plot_params$rightmargin <- 0.02
      kp <- plotKaryotype(genome=Gr_hg38, plot.type=1, plot.params = plot_params, zoom = detail.region, cex=0.8)
      kpAddBaseNumbers(kp, tick.dist = tick_dist, tick.len = 6, tick.col="#4d4d4d", cex=0.7,
                       minor.tick.dist = minor_tick_dist, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
      #kpAddMainTitle(kp, main="test", col="red")
      
      # Alu barplot
      kpDataBackground(kp, data.panel = 1, r0=r0_Alu, r1=r1_Alu, color = "#ffffff")
      kpAxis(kp, ymin=0, ymax=max(cov_data_repetitive$Alu), r0=r0_Alu, r1=r1_Alu, col="gray50", cex=label_cex, numticks = 3)
      kpBars(kp, chr=cov_data_repetitive$seqnames, x0=cov_data_repetitive$start, x1=cov_data_repetitive$end, y1=cov_data_repetitive$Alu,
             col="black", ymin=0, ymax=max(cov_data_repetitive$Alu), r0=r0_Alu, r1=r1_Alu, border = "black")
      
      # L1 barplot
      kpDataBackground(kp, data.panel = 1, r0=r0_L1, r1=r1_L1, color = "#ffffff")
      kpAxis(kp, ymin=0, ymax=max(cov_data_repetitive$L1), r0=r0_L1, r1=r1_L1, col="gray50", cex=label_cex, numticks = 3)
      kpBars(kp, chr=cov_data_repetitive$seqnames, x0=cov_data_repetitive$start, x1=cov_data_repetitive$end, y1=cov_data_repetitive$L1,
             col="gray40", ymin=0, ymax=max(cov_data_repetitive$L1), r0=r0_L1, r1=r1_L1, border = "gray40")
      
      # Repetitive log2ratio barplot
      kpDataBackground(kp, data.panel = 1, r0=r0_repetitive, r1=r1_repetitive, color = "#ffffff")
      kpAxis(kp, ymin=-y_repetitive, ymax=y_repetitive, r0=r0_repetitive, r1=r1_repetitive, col="gray50", cex=label_cex, numticks = 3)
      kpBars(kp, chr=cov_data_repetitive$seqnames, x0=cov_data_repetitive$start, x1=cov_data_repetitive$end, y1=cov_data_repetitive$log2_Alu_L1_pos,
             col="black", ymin=-y_repetitive, ymax=y_repetitive, r0=r0_repetitive, r1=r1_repetitive, border = "black")
      kpBars(kp, chr=cov_data_repetitive$seqnames, x0=cov_data_repetitive$start, x1=cov_data_repetitive$end, y1=cov_data_repetitive$log2_Alu_L1_neg,
             col="gray40", ymin=-y_repetitive, ymax=y_repetitive, r0=r0_repetitive, r1=r1_repetitive, border = "gray40")
      
      # Compartment barplot
      kpDataBackground(kp, data.panel = 1, r0=r0_comp, r1=r1_comp, color = "#ffffff")
      kpAxis(kp, ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, col="gray50", cex=label_cex, numticks = 3)
      kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_pos,
             col="#E41A1C", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, border = "#E41A1C")
      kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_neg,
             col="#0cad01", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, border = "#0cad01")
      
      
      dev.off()
      
  }
  
  return(out_cor)
}


plot_compartment_repetitive_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                               sample_name="H1_control",
                                               chromosomes=c("chr1","chr11"))

plot_compartment_repetitive_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_NH4OAc_500000.txt",
                                   sample_name="H1_NH4OAc",
                                   chromosomes=c("chr11"))

plot_compartment_repetitive_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_FL_500000.txt",
                                   sample_name="H1_FL",
                                   chromosomes=c("chr11"))

plot_compartment_repetitive_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_500000.txt",
                                   sample_name="H1_RNase",
                                   chromosomes=c("chr11"))



################ Make dataframe with RAL of repetitive elements
make_RAL_cov_data <- function(input_margi,
                              sample,
                              resolution,
                              repeat_label = NA,
                              minimum_overlap_vector=0){
  options(scipen=999)
  
  ### Load iMARGI data
  input_bedpe = data.frame(fread(input_margi))
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
  
  if (is.na(repeat_label)){
    
    for (min_overlap in minimum_overlap_vector){
      ### Alu
      repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == "SINE/Alu"),]
      repeat_data = repeat_data[which(repeat_data$V1 %in% hg38_chromosomes),]
      
      Gr_Alu = GRanges(
        seqnames = Rle(repeat_data[,1]),
        ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
        strand = Rle(strand(repeat_data[,6])))
      
      overlap = countOverlaps(Gr_RNA_end, Gr_Alu, ignore.strand = T, minoverlap = min_overlap)
      mcols(Gr_RNA_end)["overlap_Alu"] = overlap
      Gr_RNA_end_Alu = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap_Alu"] >= 1]
      Gr_DNA_end_Alu = Gr_DNA_end[names(Gr_RNA_end_Alu)]
      
      ### L1
      repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == "LINE/L1"),]
      repeat_data = repeat_data[which(repeat_data$V1 %in% hg38_chromosomes),]
      
      Gr_L1 = GRanges(
        seqnames = Rle(repeat_data[,1]),
        ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
        strand = Rle(strand(repeat_data[,6])))
      
      overlap = countOverlaps(Gr_RNA_end, Gr_L1, ignore.strand = T, minoverlap = min_overlap)
      mcols(Gr_RNA_end)["overlap_L1"] = overlap
      Gr_RNA_end_L1 = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap_L1"] >= 1]
      Gr_DNA_end_L1 = Gr_DNA_end[names(Gr_RNA_end_L1)]
      
      for (res in resolution){
        cov_data = generate_karyoploteR_data(tags = list(Gr_DNA_end_Alu,Gr_DNA_end_L1),
                                             genome_gr = Gr_hg38, 
                                             window_size = res,
                                             amplifier = c(1,1), 
                                             threshold = NA, 
                                             names = c("RAL_Alu","RAL_L1"))
        cov_data$log2_Alu_L1 = log2(cov_data$RAL_Alu / cov_data$RAL_L1)
        
        if (min_overlap == 0){
          write.table(cov_data, paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/",res,"/",sample,"_Alu_L1_RAL.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
        } else {
          write.table(cov_data, paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/",res,"/",sample,"_Alu_L1_RAL_minoverlap_",min_overlap,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
        }
      }
    }
    

    } else {
    repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == repeat_label),]
    repeat_data = repeat_data[which(repeat_data$V1 %in% hg38_chromosomes),]
    
    Gr_repeat = GRanges(
      seqnames = Rle(repeat_data[,1]),
      ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
      strand = Rle(strand(repeat_data[,6])))
    
    for (min_overlap in minimum_overlap_vector){
      overlap = countOverlaps(Gr_RNA_end, Gr_repeat, ignore.strand = T, minoverlap = min_overlap)
      mcols(Gr_RNA_end)["overlap"] = overlap
      Gr_RNA_end_repeat = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap"] >= 1]
      Gr_DNA_end_repeat = Gr_DNA_end[names(Gr_RNA_end_repeat)]
      
      for (res in resolution){
        cov_data = generate_karyoploteR_data(tags = list(Gr_DNA_end_repeat),
                                             genome_gr = Gr_hg38, 
                                             window_size = res,
                                             amplifier = c(1), 
                                             threshold = NA, 
                                             names = c(paste0("RAL_",gsub(".*/","",repeat_label))))
        if (minimum_overlap == 0){
          write.table(cov_data, paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/",res,"/",sample,"_",gsub(".*/","",repeat_label),"_RAL.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
        } else {
          write.table(cov_data, paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/",res,"/",sample,"_",gsub(".*/","",repeat_label),"_RAL_minoverlap_",min_overlap,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
        }
      }
    }
  }
  return(0)
}

for (i in c(all_imargi_samples, "iMARGI_HFF_control","iMARGI_K562_control")){
  for (j in c("1k","200k")){
    sample = gsub("iMARGI_","",i)
    make_RAL_cov_data(input_margi=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",i,"/",i,".mapq30.",j,".bedpe.gz"),
                      sample=paste0(sample,"_",j),
                      resolution=c(1000000,500000,200000,100000),
                      repeat_label = NA,
                      minimum_overlap_vector = c(0,30))
  }
}


i = "iMARGI_H1_RNaseCtrl"
for (j in c("1k","200k")){
  sample = gsub("iMARGI_","",i)
  make_RAL_cov_data(input_margi=paste0("/dataOS/frankyan/phase_separation/std_results/merged_iMARGI/",i,"/",i,".mapq30.",j,".bedpe.gz"),
                    sample=paste0(sample,"_",j),
                    resolution=c(1000000,500000,200000,100000),
                    repeat_label = NA,
                    minimum_overlap_vector = c(0,30))
}



################ Plot tracks with RAL
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
                                                   add_on_label = "" # additional custom label to be added to the figure filenames
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
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == "SINE/Alu"),]
  repeat_data = repeat_data[which(repeat_data$V1 %in% hg38_chromosomes),]
  Alu_norm_factor = sum(repeat_data$V3-repeat_data$V2) / sum(hg38_lengths)
  
  Gr_Alu = GRanges(
    seqnames = Rle(repeat_data[,1]),
    ranges = IRanges(as.numeric(repeat_data[,2]) + 1, end = as.numeric(repeat_data[,3]), names = c(1:nrow(repeat_data))),
    strand = Rle(strand(repeat_data[,6])))
  
  repeat_data = hg38_repeatMasker[which(hg38_repeatMasker$V4 == "LINE/L1"),]
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
  ma_H1 <- fread("/dataOS/wenxingzhao/database/4DN/SPIN/H1_new.SPIN.JAWG.25kb.9_state.bed")
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
      r1_spin = r0_comp - 0.06
      r0_spin = r1_spin - 0.03
      
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
  
    
    if (normalized == T){
      pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks_norm/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_norm",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    } else if (normalized_by_RAL == T){
      pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution), "/", sample_name,"_",my_chr,"_tracks_RAL_norm_by_RAL",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
    } else {
      pdf(paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_tracks/",as.character(resolution),"/",sample_name,"_",my_chr,"_tracks_RAL",add_on_label,"_minoverlap_", minoverlap,".pdf"), width = 4, height = 2.5)
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
    kpPlotRegions(kp, data=H1_spin, col=H1_spin$color, r0=r0_spin , r1=r1_spin, data.panel=1)

    
    dev.off()
    
  }
  
  return(out_cor)
}

### To pick high enough number of consecutive switches
summary(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"]) >= 3)
temp = switch_length_distribution_1000000_list[["RNase"]][["AB"]][which(switch_length_distribution_1000000_list[["RNase"]][["AB"]][,"consecutive_switches"] >= 3),]
write.table(temp, "/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/switch_AB_1000000_thres.txt", row.names = F, col.names = T, sep = "\t", quote = F)

summary(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"]) >= 3)
temp = switch_length_distribution_1000000_list[["RNase"]][["BA"]][which(switch_length_distribution_1000000_list[["RNase"]][["BA"]][,"consecutive_switches"] >= 3),]
write.table(temp, "/dataOS/rcalandrelli/phase_separation/HiC/compartments/result/switch_BA_1000000_thres.txt", row.names = F, col.names = T, sep = "\t", quote = F)




summary(as.numeric(switch_length_distribution_500000_list[["RNase"]][["AB"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_500000_list[["RNase"]][["AB"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_500000_list[["RNase"]][["AB"]][,"consecutive_switches"]) >= 4)
switch_length_distribution_500000_list[["RNase"]][["AB"]][which(switch_length_distribution_500000_list[["RNase"]][["AB"]][,"consecutive_switches"] >= 4),]


summary(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]))
quantile(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]),0.95)
sum(as.numeric(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"]) >= 3)
switch_length_distribution_500000_list[["RNase"]][["BA"]][which(switch_length_distribution_500000_list[["RNase"]][["BA"]][,"consecutive_switches"] >= 3),]





### Not normalized

for (i in c("_with_separate_RAL")){
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_control_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_control_1k",
                                         chromosomes=c("chr4","chr6","chr7","chr9","chrX"),
                                         resolution=1000000,
                                         tick_dist=50000000,
                                         minor_tick_dist=10000000,
                                         add_on_label = i)
  
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_control_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_control_1k",
                                         chromosomes=c("chr15"),
                                         resolution=1000000,
                                         tick_dist=20000000,
                                         minor_tick_dist=10000000,
                                         add_on_label = i)
  
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_control_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_control_1k",
                                         chromosomes=c("chr21"),
                                         resolution=1000000,
                                         tick_dist=10000000,
                                         minor_tick_dist=5000000,
                                         add_on_label = i)
  
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_RNase_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_RNase_1k",
                                         chromosomes=c("chr4","chr6","chr7","chr9","chrX"),
                                         resolution=1000000,
                                         tick_dist=50000000,
                                         minor_tick_dist=10000000,
                                         add_on_label = i)
  
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_RNase_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_RNase_1k",
                                         chromosomes=c("chr15"),
                                         resolution=1000000,
                                         tick_dist=20000000,
                                         minor_tick_dist=10000000,
                                         add_on_label = i)
  
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_RNase_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_RNase_1k",
                                         chromosomes=c("chr21"),
                                         resolution=1000000,
                                         tick_dist=10000000,
                                         minor_tick_dist=5000000,
                                         add_on_label = i)
  
  plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNaseCtrl_1000000.txt",
                                         repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/1000000/H1_RNaseCtrl_1k_Alu_L1_RAL.txt",
                                         sample_name="H1_RNaseCtrl_1k",
                                         chromosomes=c("chr21"),
                                         resolution=1000000,
                                         tick_dist=10000000,
                                         minor_tick_dist=5000000,
                                         add_on_label = i)
}


####### Paper plots
plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/500000/H1_control_200k_Alu_L1_RAL_minoverlap_30.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr7","chr11"),
                                       resolution=500000,
                                       tick_dist=20000000,
                                       minor_tick_dist=10000000,
                                       add_on_label = "")

plot_compartment_repetitive_RAL_tracks(input_bedpe_file = "/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                       eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/500000/H1_control_200k_Alu_L1_RAL_minoverlap_30.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr7","chr11"),
                                       resolution=500000,
                                       tick_dist=20000000,
                                       minor_tick_dist=10000000,
                                       add_on_label = "_with_global_RAL")

plot_compartment_repetitive_RAL_tracks(input_bedpe_file = "/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.200k.final.bedpe.gz",
                                       eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/500000/H1_control_200k_Alu_L1_RAL_minoverlap_30.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr7","chr11"),
                                       resolution=500000,
                                       tick_dist=20000000,
                                       minor_tick_dist=10000000,
                                       normalized_by_RAL = T,
                                       add_on_label = "_with_global_RAL")



plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/H1_control_200k_Alu_L1_RAL.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr1","chr4","chr11","chr13","chr14","chr21"))

plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/H1_RNase_200k_Alu_L1_RAL.txt",
                                       sample_name="H1_RNase_200k",
                                       chromosomes=c("chr1","chr4","chr11","chr13","chr14","chr21"))




### Normalized
plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/H1_control_1k_Alu_L1_RAL.txt",
                                       sample_name="H1_control_1k",
                                       chromosomes=c("chr1","chr4","chr11","chr13","chr14","chr21"),
                                       normalized = T)

plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/H1_RNase_1k_Alu_L1_RAL.txt",
                                       sample_name="H1_RNase_1k",
                                       chromosomes=c("chr1","chr4","chr11","chr13","chr14","chr21"),
                                       normalized = T)



plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/H1_control_200k_Alu_L1_RAL.txt",
                                       sample_name="H1_control_200k",
                                       chromosomes=c("chr1","chr4","chr11","chr13","chr14","chr21"),
                                       normalized = T)

plot_compartment_repetitive_RAL_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_RNase_500000.txt",
                                       repetitive_RAL_file = "/dataOS/rcalandrelli/phase_separation/HiC/compartments/repetitive_elements/RAL_repetitive/H1_RNase_200k_Alu_L1_RAL.txt",
                                       sample_name="H1_RNase_200k",
                                       chromosomes=c("chr1","chr4","chr11","chr13","chr14","chr21"),
                                       normalized = T)




















