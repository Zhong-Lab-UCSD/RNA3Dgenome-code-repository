###################### iMARGI analysis at the compartment level
library(grid)
library(ggplotify)
library(cowplot)
library(png)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(scales)
library(data.table)
library(regioneR)
library("karyoploteR", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")

options(scipen=999)


chrSize = read.table("/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes")
hg38_chromosomes = as.character(chrSize$V1)
hg38_lengths = as.numeric(chrSize$V2)
names(hg38_lengths) = hg38_chromosomes

Gr_hg38 <- GRanges(
  seqnames = Rle(names(hg38_lengths[1:23])),
  ranges = IRanges(rep(1,23), end = as.numeric(hg38_lengths[1:23]), names = c(1:length(hg38_lengths[1:23]))),
  strand = Rle(strand('*')))
seqlengths(Gr_hg38) <- hg38_lengths[names(seqlengths(Gr_hg38))]

annotation_hg38 <- read.table('/dataOS/rcalandrelli/MARGI/Homo_sapiens.GRCh38.84.chr.gtf_to_geneTable.tsv', stringsAsFactors = F, header = T)

# # Load compartment data
# eigen_data_list = list()
# eigen_control = read.table("/dataOS/rcalandrelli/phase_separation/compartments/H1_control_500000.txt", header = T)
# eigen_NH4OAc = read.table("/dataOS/rcalandrelli/phase_separation/compartments/H1_NH4OAc_500000.txt", header = T)
# eigen_FL = read.table("/dataOS/rcalandrelli/phase_separation/compartments/H1_FL_500000.txt", header = T)
# eigen_RNase = read.table("/dataOS/rcalandrelli/phase_separation/compartments/H1_RNase_500000.txt", header = T)
# 
# eigen_data_list[["control"]] = eigen_control
# eigen_data_list[["NH4OAc"]] = eigen_NH4OAc
# eigen_data_list[["FL"]] = eigen_FL
# eigen_data_list[["RNase"]] = eigen_RNase
# 
# # Initialize iMARGI sample names
# imargi_samples = c(paste0("iMARGI_H1_control_",seq(1,5)),
#                    paste0("iMARGI_H1_NH4OAc_",seq(1,2)),
#                    paste0("iMARGI_H1_FL_",seq(1,2)),
#                    paste0("iMARGI_H1_RNase_",seq(1,2)))
# 
# # Load iMARGI data
# load_iMARGI <- function(path, dist_filter=NA){
#   temp = read.table(path)
#   temp = temp[which(temp[,1] %in% hg38_chromosomes & temp[,4] %in% hg38_chromosomes),]
#   temp[,1] = as.character(temp[,1])
#   temp[,4] = as.character(temp[,4])
#   # if (!is.na(dist_filter)){
#   #   temp = temp[which((temp[,1] == temp[,4] &
#   #                       abs(temp[,2]-temp[,5]) > dist_filter) |
#   #                       temp[,1] != temp[,4]),]
#   # }
#   return(temp)
# }
# 
# # Create a list with iMARGI contact data in each element
# imargi_data_list_1k = list()
# for (i in imargi_samples[3:4]){
#   imargi_data_list_1k[[i]] = load_iMARGI(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",i,"/","filter1k_",i,".txt"))
# }

### Plot A/B compartment track with karyoploteR
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

plot_compartment_imargi_tracks <- function(eigen_file,
                                           input_bedpe_file,
                                           sample_name,
                                           chromosomes,
                                           resolution,
                                           start_coord=0,
                                           end_coord=0,
                                           tick_dist=50000000,
                                           minor_tick_dist=10000000,
                                           log_margi=F,
                                           MALAT1_cov=F,
                                           y_max_RNA=NA,
                                           y_max_DNA=NA){
  
  ### Load compartment data
  compartment = read.table(eigen_file, header = T)
  compartment = cbind(compartment,compartment[,3])
  compartment[,3][which(compartment[,3]<0)] = NA
  compartment[,4][which(compartment[,4]>=0)] = NA
  colnames(compartment) = c("chr", "coord", "eigen_pos", "eigen_neg")
  #compartment = compartment[which(compartment$chr != "chrY"),]
  
  ### Load iMARGI data
  input_bedpe = data.frame(fread(input_bedpe_file))
  input_data = input_bedpe[which(input_bedpe[,1] %in% chromosomes & 
                                   input_bedpe[,4] %in% chromosomes),1:6]
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
  
  if (MALAT1_cov == T){
    malat1 = annotation_hg38[which(annotation_hg38$gene_name=="MALAT1"),][1,]
    Gr_malat1 = GRanges(
      seqnames = Rle(as.character(malat1[1])),
      ranges = IRanges(as.numeric(malat1[2]), end = as.numeric(malat1[3]), names = 1),
      strand = Rle(strand('*')))
    
    overlap = countOverlaps(Gr_RNA_end, Gr_malat1, ignore.strand = T)
    mcols(Gr_RNA_end)["overlap_malat1"] = overlap
    Gr_RNA_end = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap_malat1"] >= 1]
    
    Gr_DNA_end = Gr_DNA_end[names(Gr_RNA_end)]
  }
  

  cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
                                       genome_gr = Gr_hg38, 
                                       window_size = resolution,
                                       amplifier = c(1,1), 
                                       threshold = NA, 
                                       names = c("RNA_end","DNA_end"))
  
  if (log_margi == T){
    cov_data$RNA_end = log(cov_data$RNA_end)
    cov_data$RNA_end[which(is.infinite(cov_data$RNA_end))] = 0

    # cov_data$DNA_end = log(cov_data$DNA_end)
    # cov_data$DNA_end[which(is.infinite(cov_data$DNA_end))] = 0
  }
  
  for (my_chr in chromosomes){
    print(paste0(sample_name, " - ", my_chr))
    start_coord = 0
    end_coord = hg38_lengths[my_chr]
    detail.region <- toGRanges(data.frame(my_chr, start_coord, end_coord))
    
    y_text = 0.1 # track label height
    
    # Compartment area
    r1_comp = 0.9
    r0_comp = r1_comp - 0.2
    ymin_comp = round(min(compartment[which(compartment$chr==my_chr),4][which(!is.na(compartment[which(compartment$chr==my_chr),4]))]),1)
    ymax_comp = round(max(compartment[which(compartment$chr==my_chr),3][which(!is.na(compartment[which(compartment$chr==my_chr),3]))]),1)
    
    # RNA_end area
    r1_RNA_end = r0_comp - 0.1
    r0_RNA_end = r1_RNA_end - 0.2
    if (is.na(y_max_RNA)) {
      y_max_RNA_end = round(max(cov_data[which(cov_data$seqnames==my_chr),"RNA_end"]))
    } else {
      y_max_RNA_end = y_max_RNA
    }
    
    # DNA_end area
    r1_DNA_end = r0_RNA_end - 0.1
    r0_DNA_end = r1_DNA_end - 0.2
    if (is.na(y_max_DNA)) {
      y_max_DNA_end = round(max(cov_data[which(cov_data$seqnames==my_chr),"DNA_end"]))
    } else {
      y_max_DNA_end = y_max_DNA
    }
    
    line_width = 1 # width of the line track
    label_cex = 0.5
    
    png(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/iMARGI_",sample_name,"/", my_chr,"_tracks.png"), width = 5.5, height = 2, units = "in", res = 400)
    plot_params <- getDefaultPlotParams(plot.type=1)
    plot_params$data1inmargin <- 5
    plot_params$ideogramheight <- 10
    plot_params$bottommargin <- 30
    plot_params$topmargin <- 0
    plot_params$leftmargin <- 0.11
    plot_params$rightmargin <- 0.02
    kp <- plotKaryotype(genome="hg38", plot.type=1, plot.params = plot_params, zoom = detail.region, cex=0.7)
    kpAddBaseNumbers(kp, tick.dist = tick_dist, tick.len = 6, tick.col="#4d4d4d", cex=0.6,
                     minor.tick.dist = minor_tick_dist, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
    #kpAddMainTitle(kp, main="test", col="red")
    
    # Compartment barplot
    kpDataBackground(kp, data.panel = 1, r0=r0_comp, r1=r1_comp, color = "#ffffff")
    kpAxis(kp, ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp, col="gray50", cex=label_cex, numticks = 3)
    kpText(kp, chr=compartment$chr, x=mean(c(start(detail.region),end(detail.region))), 
           y=y_text, col="#000000", r0=r1_comp + 0.04, r1=r1_comp + 0.06, labels="A/B compartments", cex=label_cex)
    kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_pos,
           col="#E41A1C", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp)
    kpBars(kp, chr=compartment$chr, x0=compartment$coord, x1=compartment$coord+resolution, y1=compartment$eigen_neg,
           col="#0cad01", ymin=ymin_comp, ymax=ymax_comp, r0=r0_comp, r1=r1_comp)
    kpAddLabels(kp, sample_name, label.margin=0.095, side="left", pos=NULL, offset=2, r0=0.5, r1=1, data.panel=1, srt = 90, cex = 0.7)
    
    # RNA end track
    kpDataBackground(kp, data.panel = 2, r0=r0_RNA_end, r1=r1_RNA_end, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=y_max_RNA_end, r0=r0_RNA_end, r1=r1_RNA_end, col="gray50", cex=label_cex, numticks = 3)
    kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
           y=y_text, col="#000000", r0=r1_RNA_end + 0.04, r1=r1_RNA_end + 0.06, labels="iMARGI RNA end", cex=label_cex)
    kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$RNA_end,
            col="blue", ymin=0, ymax=y_max_RNA_end, r0=r0_RNA_end, r1=r1_RNA_end, lwd=line_width)
    
    # DNA end track
    kpDataBackground(kp, data.panel = 2, r0=r0_DNA_end, r1=r1_DNA_end, color = "#ffffff")
    kpAxis(kp, ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, col="gray50", cex=label_cex, numticks = 3)
    kpText(kp, chr=cov_data$seqnames, x=mean(c(start(detail.region),end(detail.region))), 
           y=y_text, col="#000000", r0=r1_DNA_end + 0.04, r1=r1_DNA_end + 0.06, labels="iMARGI DNA end", cex=label_cex)
    kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$DNA_end,
            col="purple", ymin=0, ymax=y_max_DNA_end, r0=r0_DNA_end, r1=r1_DNA_end, lwd=line_width)
    dev.off()
  }
}

################################### Plotting

all_imargi_samples = list.dirs(path = "/dataOS/rcalandrelli/phase_separation/MARGI/data", full.names = F, recursive = F)

##### All read pairs (with the same y limits than control)

for (sample in all_imargi_samples[1:4]){
  sample_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/", sample)
  if (!dir.exists(sample_dir)){
    dir.create(sample_dir)
  }
  sample_name=gsub("iMARGI_","",sample)
  temp = strsplit(sample,"_")[[1]][3]
  plot_compartment_imargi_tracks(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", temp, "_500000.txt"),
                                 input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                 sample_name=sample_name,
                                 chromosomes=hg38_chromosomes[1:23],
                                 resolution=500000,
                                 #start_coord=0,
                                 #end_coord=0,
                                 tick_dist=20000000,
                                 minor_tick_dist=5000000,
                                 log_margi = F,
                                 MALAT1_cov=F,
                                 y_max_RNA=2821931, # max FL
                                 y_max_DNA=71326) # max control
}

figures_list = list()
k = 1
for (i in all_imargi_samples[1:4]){
  temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/all_reads_same_y/noLog/",i,"/chr11_tracks.png"))
  figures_list[[k]] = temp
  k = k + 1
}

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/all_reads_same_y/tracks_plots_noLog_chr11.png", width = 3, height = 1, units = "in", res = 1000)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]), ncol=2)
dev.off()


##### All read pairs

for (sample in all_imargi_samples[1:4]){
  sample_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/", sample)
  if (!dir.exists(sample_dir)){
    dir.create(sample_dir)
  }
  sample_name=gsub("iMARGI_","",sample)
  temp = strsplit(sample,"_")[[1]][3]
  plot_compartment_imargi_tracks(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", temp, "_500000.txt"),
                                 input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                 sample_name=sample_name,
                                 chromosomes=hg38_chromosomes[1:23],
                                 resolution=500000,
                                 #start_coord=0,
                                 #end_coord=0,
                                 tick_dist=20000000,
                                 minor_tick_dist=5000000,
                                 log_margi = T)
}

figures_list = list()
k = 1
for (i in all_imargi_samples[1:4]){
  # temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/logRNA/",i,"/chr11_tracks.png"))
  temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/loglog/",i,"/chr11_tracks.png"))
  # temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/noLog/",i,"/chr11_tracks.png"))
  figures_list[[k]] = temp
  k = k + 1
}

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/all_reads/tracks_plots_loglog_chr11.png", width = 3, height = 1, units = "in", res = 1000)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]), ncol=2)
dev.off()


##### Proximal read pairs
for (sample in all_imargi_samples[1:4]){
  sample_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/", sample)
  if (!dir.exists(sample_dir)){
    dir.create(sample_dir)
  }
  sample_name=gsub("iMARGI_","",sample)
  temp = strsplit(sample,"_")[[1]][3]
  plot_compartment_imargi_tracks(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", temp, "_500000.txt"),
                                 input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k_200k.final.bedpe.gz"),
                                 sample_name=sample_name,
                                 chromosomes=hg38_chromosomes[1:23],
                                 resolution=500000,
                                 #start_coord=0,
                                 #end_coord=0,
                                 tick_dist=20000000,
                                 minor_tick_dist=5000000,
                                 log_margi = F)
}

figures_list = list()
k = 1
for (i in all_imargi_samples[1:4]){
  #temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/proximal_reads/logRNA/",i,"/chr11_tracks.png"))
  temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/proximal_reads/noLog/",i,"/chr11_tracks.png"))
  figures_list[[k]] = temp
  k = k + 1
}

#png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/proximal_reads/tracks_plots_logRNA_chr11.png", width = 3, height = 1, units = "in", res = 1000)
png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/proximal_reads/tracks_plots_noLog_chr11.png", width = 3, height = 1, units = "in", res = 1000)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]), ncol=2)
dev.off()


##### Distal read pairs
for (sample in all_imargi_samples[1:4]){
  sample_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/", sample)
  if (!dir.exists(sample_dir)){
    dir.create(sample_dir)
  }
  sample_name=gsub("iMARGI_","",sample)
  temp = strsplit(sample,"_")[[1]][3]
  plot_compartment_imargi_tracks(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", temp, "_500000.txt"),
                                 input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz"),
                                 sample_name=sample_name,
                                 chromosomes=hg38_chromosomes[1:23],
                                 resolution=500000,
                                 #start_coord=0,
                                 #end_coord=0,
                                 tick_dist=20000000,
                                 minor_tick_dist=5000000,
                                 log_margi = T)
}


figures_list = list()
k = 1
for (i in all_imargi_samples[1:4]){
  temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/distal_reads/logRNA/",i,"/chr11_tracks.png"))
  #temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/distal_reads/noLog/",i,"/chr11_tracks.png"))
  figures_list[[k]] = temp
  k = k + 1
}

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/distal_reads/tracks_plots_logRNA_chr11.png", width = 3, height = 1, units = "in", res = 1000)
#png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/distal_reads/tracks_plots_noLog_chr11.png", width = 3, height = 1, units = "in", res = 1000)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]), ncol=2)
dev.off()


##### MALAT1 targets
for (sample in all_imargi_samples[1:4]){
  sample_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/", sample)
  if (!dir.exists(sample_dir)){
    dir.create(sample_dir)
  }
  sample_name=gsub("iMARGI_","",sample)
  temp = strsplit(sample,"_")[[1]][3]
  plot_compartment_imargi_tracks(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", temp, "_500000.txt"),
                                 input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz"),
                                 sample_name=sample_name,
                                 chromosomes=hg38_chromosomes[1:23],
                                 resolution=500000,
                                 #start_coord=0,
                                 #end_coord=0,
                                 tick_dist=20000000,
                                 minor_tick_dist=5000000,
                                 log_margi = F,
                                 MALAT1_cov = T)
}


figures_list = list()
k = 1
for (i in all_imargi_samples[1:4]){
  #temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/malat1/logRNA/",i,"/chr11_tracks.png"))
  temp = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/malat1_distal_reads/noLog/",i,"/chr11_tracks.png"))
  figures_list[[k]] = temp
  k = k + 1
}

#png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/malat1/tracks_plots_logRNA_chr11.png", width = 3, height = 1, units = "in", res = 1000)
png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/tracks_plots/malat1_distal_reads/tracks_plots_noLog_chr11.png", width = 3, height = 1, units = "in", res = 1000)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]), ncol=2)
dev.off()




# plot_compartment_imargi_tracks(eigen_file="/dataOS/rcalandrelli/phase_separation/compartments/H1_control_500000.txt",
#                                input_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control_4/filter200k_iMARGI_H1_control_4.final.bedpe",
#                                sample_name="test",
#                                my_chr="chr11",
#                                resolution=500000,
#                                start_coord=0,
#                                end_coord=hg38_lengths["chr11"],
#                                tick_dist=20000000,
#                                minor_tick_dist=5000000,
#                                log_margi = F,
#                                out_filename=paste0("sample",".png"))



######################## Calculate correlation between eigenvector and read coverage on the RNA and DNA ends
compartment_margi_correlation <- function(eigen_file,
                                          input_bedpe_file,
                                          resolution,
                                          log_margi=F,
                                          fig_title,
                                          force_compartment_control=F,
                                          force_imargi_control=F,
                                          MALAT1_cov=F){
  
  if (force_compartment_control == F){
    compartment = read.table(eigen_file, header = T)
    colnames(compartment) = c("chr","coord","eigen")
    compartment = compartment[which(compartment$chr != "chrY"),]
  } else {
    compartment = read.table("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt", header = T)
    colnames(compartment) = c("chr","coord","eigen")
    compartment = compartment[which(compartment$chr != "chrY"),]
  }

  if (force_imargi_control == F){
    input_bedpe = data.frame(fread(input_bedpe_file))
  } else {
    input_bedpe = data.frame(fread("/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe.gz"))
  }
  
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
  
  if (MALAT1_cov == T){
    malat1 = annotation_hg38[which(annotation_hg38$gene_name=="MALAT1"),][1,]
    Gr_malat1 = GRanges(
      seqnames = Rle(as.character(malat1[1])),
      ranges = IRanges(as.numeric(malat1[2]), end = as.numeric(malat1[3]), names = 1),
      strand = Rle(strand('*')))
    
    overlap = countOverlaps(Gr_RNA_end, Gr_malat1, ignore.strand = T)
    mcols(Gr_RNA_end)["overlap_malat1"] = overlap
    Gr_RNA_end = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap_malat1"] >= 1]
    
    Gr_DNA_end = Gr_DNA_end[names(Gr_RNA_end)]
  }
  
  cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
                                       genome_gr = Gr_hg38, 
                                       window_size = resolution,
                                       amplifier = c(1,1), 
                                       threshold = NA, 
                                       names = c("RNA_end","DNA_end"))
  

  # rm(input_data)
  # rm(Gr_RNA_end)
  # rm(Gr_DNA_end)
  
  # if (log_margi == T){
  #   cov_data$RNA_end = log(cov_data$RNA_end)
  #   cov_data$RNA_end[which(is.infinite(cov_data$RNA_end))] = 0
  # 
  #   cov_data$DNA_end = log(cov_data$DNA_end)
  #   cov_data$DNA_end[which(is.infinite(cov_data$DNA_end))] = 0
  # }
  
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
    vjustvar = c(-0.5))
  
  label_size = 12
  
  y_lab = ifelse(force_compartment_control == T, "Eigenvector control", "Eigenvector")
  
  p1<-ggplot(df, aes(x=cov_data.RNA_end, y=compartment.eigen)) +
    geom_point(color = "blue") +
    xlab("RNA end read coverage") +
    ylab(y_lab) +
    scale_x_continuous(labels = scientific) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size)) +
    ggtitle(fig_title) +
    geom_text(data=annotations_RNA,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  
  p2<-ggplot(df, aes(x=log(cov_data.RNA_end), y=compartment.eigen)) + 
    geom_point(color = "blue") +
    xlab("log(RNA end read coverage)") +
    ylab(y_lab) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size)) +
    ggtitle(fig_title) +
    geom_text(data=annotations_RNA,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  
  jj
  p3<-ggplot(df, aes(x=cov_data.DNA_end, y=compartment.eigen)) + 
    geom_point(color = "purple") +
    #xlab("DNA end read coverage") +
    xlab("RAL") +
    ylab(y_lab) +
    scale_x_continuous(labels = scientific) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size)) +
    ggtitle(fig_title) +
    geom_text(data=annotations_DNA,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  
  p4 <- plot_grid(p1,p3) # noLog
  p5 <- plot_grid(p2,p3) # logRNA
  
  png(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/scatter_plots/iMARGI_",fig_title,"_noLog.png"), width = 8, height = 4, res = 200, units = "in")
  print(p4)
  dev.off()
  
  png(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/scatter_plots/iMARGI_",fig_title,"_logRNA.png"), width = 8, height = 4, res = 200, units = "in")
  print(p5)
  dev.off()
  
  return(list(output_list,p4,p5))
}


# eigen_file="/dataOS/rcalandrelli/phase_separation/compartments/H1_control_500000.txt"
# input_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control_4/filter200k_iMARGI_H1_control_4.final.bedpe"
# resolution=500000
# log_margi=F
# fig_title="test"

all_imargi_samples = c("iMARGI_H1_control", "iMARGI_H1_FL", "iMARGI_H1_NH4OAc", "iMARGI_H1_RNase")
all_imargi_bio_replicates = c(paste0("iMARGI_H1_control_", c(1,3,4,5)),
                             paste0("iMARGI_H1_FL_", c(1,2)),
                             paste0("iMARGI_H1_NH4OAc_", c(1,2)),
                             paste0("iMARGI_H1_RNase_", c(1,2)))


############# Merged replicates
correlation_list = list()
correlation_list_compartment_control = list()
correlation_list_imargi_control = list()

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                        input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                        resolution=500000,
                                                        log_margi = F,
                                                        fig_title=sample_name,
                                                        force_compartment_control=F,
                                                        force_imargi_control=F)
}

temp_list_control = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                  input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                  resolution=500000,
                                                  log_margi = F,
                                                  fig_title=sample_name,
                                                  force_compartment_control=F,
                                                  force_imargi_control=F)

summary(temp_list_control[[1]][[2]])

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_compartment_control[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                             input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                             resolution=500000,
                                                             log_margi = F,
                                                             fig_title=sample_name,
                                                             force_compartment_control=T,
                                                             force_imargi_control=F)
}

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_imargi_control[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                            input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                                            resolution=500000,
                                                                            log_margi = F,
                                                                            fig_title=sample_name,
                                                                            force_compartment_control=F,
                                                                            force_imargi_control=T)
  
}


############# Biological replicates
correlation_list_bio_rep = list()
correlation_list_compartment_control_bio_rep = list()
correlation_list_imargi_control_bio_rep = list()

for (sample in all_imargi_bio_replicates){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_bio_rep[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                             input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                             resolution=500000,
                                                             log_margi = F,
                                                             fig_title=sample_name,
                                                             force_compartment_control=F,
                                                             force_imargi_control=F)
}

for (sample in all_imargi_bio_replicates){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_compartment_control_bio_rep[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                                 input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                                                 resolution=500000,
                                                                                 log_margi = F,
                                                                                 fig_title=sample_name,
                                                                                 force_compartment_control=T,
                                                                                 force_imargi_control=F)
}

for (sample in all_imargi_bio_replicates){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_imargi_control_bio_rep[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                            input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                                            resolution=500000,
                                                                            log_margi = F,
                                                                            fig_title=sample_name,
                                                                            force_compartment_control=F,
                                                                            force_imargi_control=T)
  
}


############# Merged replicates proximal and distal separately
correlation_list_proximal = list()
correlation_list_distal = list()

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_proximal[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                     input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k_200k.final.bedpe.gz"),
                                                                     resolution=500000,
                                                                     log_margi = F,
                                                                     fig_title=sample_name,
                                                                     force_compartment_control=F,
                                                                     force_imargi_control=F)
}

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_distal[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                                         input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz"),
                                                                                         resolution=500000,
                                                                                         log_margi = F,
                                                                                         fig_title=sample_name,
                                                                                         force_compartment_control=F,
                                                                                         force_imargi_control=F)
}


############# MALAT1 targets
correlation_list_malat1 = list()
correlation_list_malat1_distal = list()

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_malat1[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                      input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                                      resolution=500000,
                                                                      log_margi = F,
                                                                      fig_title=sample_name,
                                                                      force_compartment_control=F,
                                                                      force_imargi_control=F,
                                                                      MALAT1_cov = T)
}

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  correlation_list_malat1_distal[[sample]] = compartment_margi_correlation(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                                    input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.200k.final.bedpe.gz"),
                                                                    resolution=500000,
                                                                    log_margi = F,
                                                                    fig_title=sample_name,
                                                                    force_compartment_control=F,
                                                                    force_imargi_control=F,
                                                                    MALAT1_cov = T)
}


#################################### Parse dataframe with all the correlation values
corr_data = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list)){
  corr_data = rbind(corr_data, correlation_list[[i]][1][[1]])
}
rownames(corr_data) = names(correlation_list)
write.table(corr_data, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data.txt", row.names = T, col.names = T, sep="\t", quote = F)


corr_data_force_comp = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_compartment_control)){
  corr_data_force_comp = rbind(corr_data_force_comp, correlation_list_compartment_control[[i]][1][[1]])
}
rownames(corr_data_force_comp) = names(correlation_list_compartment_control)
write.table(corr_data_force_comp, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_force_compartment_control.txt", row.names = T, col.names = T, sep="\t", quote = F)


corr_data_force_imargi = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_imargi_control)){
  corr_data_force_imargi = rbind(corr_data_force_imargi, correlation_list_imargi_control[[i]][1][[1]])
}
rownames(corr_data_force_imargi) = names(correlation_list_imargi_control)
write.table(corr_data_force_imargi, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_force_imargi_control.txt", row.names = T, col.names = T, sep="\t", quote = F)


corr_data_proximal = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_proximal)){
  corr_data_proximal = rbind(corr_data_proximal, correlation_list_proximal[[i]][1][[1]])
}
rownames(corr_data_proximal) = names(correlation_list_proximal)
write.table(corr_data_proximal, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_proximal.txt", row.names = T, col.names = T, sep="\t", quote = F)


corr_data_distal = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_distal)){
  corr_data_distal = rbind(corr_data_distal, correlation_list_distal[[i]][1][[1]])
}
rownames(corr_data_distal) = names(correlation_list_distal)
write.table(corr_data_distal, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_distal.txt", row.names = T, col.names = T, sep="\t", quote = F)


corr_data_malat1 = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_malat1)){
  corr_data_malat1 = rbind(corr_data_malat1, correlation_list_malat1[[i]][1][[1]])
}
rownames(corr_data_malat1) = names(correlation_list_malat1)
write.table(corr_data_malat1, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_malat1.txt", row.names = T, col.names = T, sep="\t", quote = F)


corr_data_malat1_distal = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_malat1_distal)){
  corr_data_malat1_distal = rbind(corr_data_malat1_distal, correlation_list_malat1_distal[[i]][1][[1]])
}
rownames(corr_data_malat1_distal) = names(correlation_list_malat1_distal)
write.table(corr_data_malat1_distal, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_malat1_distal.txt", row.names = T, col.names = T, sep="\t", quote = F)


################################# Plot correlation data
df_corr = as.data.frame(corr_data)
df_corr$sample = gsub("iMARGI_","",rownames(df_corr))
df_corr$sample_type = sapply(df_corr$sample, function(x){paste0("H1_",strsplit(x,"_")[[1]][2])})
df = melt(df_corr)

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_barplot.png", width = 7, height = 5, res = 200, units = "in")
ggplot(df, aes(x=variable, y=value, fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.7)) +
  xlab("") +
  ylab("Correlation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_blank())
dev.off()

##### Plot correlation data to compare different correlation analysis
df_corr = as.data.frame(corr_data)
df_corr$sample = gsub("iMARGI_","",rownames(df_corr))
df_corr$type = "Within condition"
df_corr_melt = melt(df_corr)

df_corr_comp_control = as.data.frame(corr_data_force_comp)
df_corr_comp_control$sample = gsub("iMARGI_","",rownames(df_corr_comp_control))
df_corr_comp_control$type = "Compartment control"
df_corr_comp_control_melt = melt(df_corr_comp_control)

df_corr_imargi_control = as.data.frame(corr_data_force_imargi)
df_corr_imargi_control$sample = gsub("iMARGI_","",rownames(df_corr_imargi_control))
df_corr_imargi_control$type = "iMARGI control"
df_corr_imargi_control_melt = melt(df_corr_imargi_control)

df = rbind(df_corr_melt,df_corr_comp_control_melt,df_corr_imargi_control_melt)
df$width = 1

### Add a bar for the SCC of Hi-C's eigenvalues between each treatment and the control
SCC_FL_control_eigen = cor(eigen_control_500000$eigen, eigen_FL_500000$eigen, method = "spearman")
SCC_NH4OAc_control_eigen = cor(eigen_control_500000$eigen, eigen_NH4OAc_500000$eigen, method = "spearman")
SCC_RNase_control_eigen = cor(eigen_control_500000$eigen, eigen_RNase_500000$eigen, method = "spearman")

df1 = data.frame(sample=c("H1_FL","H1_NH4OAc","H1_RNase"),
                 type="SCC_eigen",
                 variable="SCC_eigen",
                 value=c(SCC_FL_control_eigen,SCC_NH4OAc_control_eigen,SCC_RNase_control_eigen),
                 width = 0.33)

df = rbind(df,df1)

### Remove PCC
df = df[which(df$variable != "PCC_RNA" & df$variable != "PCC_DNA"),]
df$type = factor(df$type, levels = c("Within condition", "Compartment control", "iMARGI control"))

p1<-ggplot(df[which(df$sample == "H1_FL"),], aes(x=variable, y=value, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("Correlation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_blank()) +
  ggtitle("H1_FL") +
  geom_hline(aes(yintercept = 0.5183008, color = "H1_control_SCC_RNA")) +
  geom_hline(aes(yintercept = 0.6339987, color = "H1_control_SCC_DNA")) +
  geom_hline(aes(yintercept = SCC_FL_control_eigen, color = "SCC_eigen_control")) +
  scale_color_manual(name = "", values = c(H1_control_SCC_RNA = "blue", H1_control_SCC_DNA = "purple", SCC_eigen_control = "red"))


p2<-ggplot(df[which(df$sample == "H1_NH4OAc"),], aes(x=variable, y=value, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("Correlation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_blank()) +
  ggtitle("H1_NH4OAc") +
  geom_hline(aes(yintercept = 0.5183008, color = "H1_control_SCC_RNA")) +
  geom_hline(aes(yintercept = 0.6339987, color = "H1_control_SCC_DNA")) +
  geom_hline(aes(yintercept = SCC_NH4OAc_control_eigen, color = "SCC_eigen_control")) +
  scale_color_manual(name = "", values = c(H1_control_SCC_RNA = "blue", H1_control_SCC_DNA = "purple", SCC_eigen_control = "red"))

p3<-ggplot(df[which(df$sample == "H1_RNase"),], aes(x=variable, y=value, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("Correlation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_blank()) +
  ggtitle("H1_RNase") +
  geom_hline(aes(yintercept = 0.5183008, color = "H1_control_SCC_RNA")) +
  geom_hline(aes(yintercept = 0.6339987, color = "H1_control_SCC_DNA")) +
  geom_hline(aes(yintercept = SCC_FL_control_eigen, color = "SCC_eigen_control")) +
  scale_color_manual(name = "", values = c(H1_control_SCC_RNA = "blue", H1_control_SCC_DNA = "purple", SCC_eigen_control = "red"))

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_barplot_across_conditions.png", width = 16, height = 4, res = 200, units = "in")
plot_grid(p1,p2,p3, nrow = 1)
dev.off()


### Plot by sample with error bars
corr_data_bio_rep = matrix(nrow=0,ncol=4)
for (i in 1:length(correlation_list_bio_rep)){
  corr_data_bio_rep = rbind(corr_data_bio_rep, correlation_list_bio_rep[[i]][1][[1]])
}
rownames(corr_data_bio_rep) = names(correlation_list_bio_rep)
write.table(corr_data_bio_rep, "/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_data_bio_replicate.txt", row.names = T, col.names = T, sep="\t", quote = F)

df_corr_bio_rep = as.data.frame(corr_data_bio_rep)
df_corr_bio_rep$sample = gsub("iMARGI_","",rownames(df_corr_bio_rep))
df_corr_bio_rep$sample_type = sapply(df_corr_bio_rep$sample, function(x){paste0("H1_",strsplit(x,"_")[[1]][2])})

df_corr_mean = aggregate(df_corr_bio_rep[,1:4],
                        by = list(sample_type = df_corr_bio_rep$sample_type),
                        FUN="mean")
df1=melt(df_corr_mean)

a = aggregate(df_corr_bio_rep[,1:4],
          by = list(sample_type = df_corr_bio_rep$sample_type),
          FUN="sd")
df_corr_error = cbind(a[,2:5] / sqrt(c(4,2,2,2)))
df_corr_error$sample_type = a[,1]
df2=melt(df_corr_error)

df = merge(df1,df2,by=c("sample_type","variable"))
colnames(df)[3:4] = c("value","error")

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/correlation_barplot_error_bars.png", width = 6, height = 5, res = 200, units = "in")
ggplot(df, aes(fill=sample_type, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.7)) +
  labs(x = "", y = "Correlation") +
  geom_errorbar(aes(ymin=value-error, ymax=value+error), width=.2, size = 0.5, position=position_dodge(.7)) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()


### Make scatter plots
p1 <- correlation_list[[1]][2][[1]]
p2 <- correlation_list[[2]][2][[1]]
p3 <- correlation_list[[3]][2][[1]]
p4 <- correlation_list[[4]][2][[1]]

p5 <- correlation_list[[1]][3][[1]]
p6 <- correlation_list[[2]][3][[1]]
p7 <- correlation_list[[3]][3][[1]]
p8 <- correlation_list[[4]][3][[1]]

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/all_scatter_plots_noLog.png", width = 15, height = 7, res = 200, units = "in")
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/all_scatter_plots_logRNA.png", width = 15, height = 7, res = 200, units = "in")
plot_grid(p5,p6,p7,p8, ncol = 2)
dev.off()

p1 <- correlation_list_compartment_control[[1]][3][[1]]
p2 <- correlation_list_compartment_control[[2]][3][[1]]
p3 <- correlation_list_compartment_control[[3]][3][[1]]
p4 <- correlation_list_compartment_control[[4]][3][[1]]

png("/dataOS/rcalandrelli/phase_separation/MARGI/compartments/all_scatter_plots_compartment_control_logRNA.png", width = 15, height = 7, res = 200, units = "in")
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()


######################## Calculate correlation between iMARGI and control HiC genome wise coverage
make_genome_wide_hic_coverage <- function(HiC_matrices_path,
                                          bin_size){
  
  options(scipen=999)
  ########## Load genome data
  chrSize = read.table(chrSize_file)
  chrSize$V1 = as.character(chrSize$V1)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ##########  Tile genome and make genomic windows for selected chromosomes
  
  ########## Load contact matrix
  genome_list = list()
  for (i in 1:(length(chrSize$V1)-1)){ # no chrY 
    chr_list = list()
    for (j in 1:(length(chrSize$V1)-1)){
      print(paste0(chrSize$V1[i],"_",chrSize$V1[j]))
      genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = bin_size, cut.last.tile.in.chrom = T)
      chr_row_genome_window = genome_window[seqnames(genome_window)==chrSize$V1[i]]
      chr_col_genome_window = genome_window[seqnames(genome_window)==chrSize$V1[j]]
      
      temp = data.frame(fread(paste0(HiC_matrices_path,"/",chrSize$V1[i],"_",chrSize$V1[j],"_",bin_size,".txt")))
      contact_matrix = matrix(0, nrow = length(chr_row_genome_window), ncol = length(chr_col_genome_window))
      if (i == j){
        contact_matrix[as.matrix(temp[,c("V1","V2")] / bin_size + 1)] = temp$V3
        contact_matrix[as.matrix(temp[,c("V2","V1")] / bin_size + 1)] = temp$V3
      } else if (i < j) {
        contact_matrix[as.matrix(temp[,c("V1","V2")] / bin_size + 1)] = temp$V3
      } else if (i > j) {
        contact_matrix[as.matrix(temp[,c("V2","V1")] / bin_size + 1)] = temp$V3
      }
      contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
      chr_list[[chrSize$V1[j]]] = contact_matrix
    }
    genome_list[[chrSize$V1[i]]] = do.call("cbind", chr_list)
  }
  
  genome_matrix = do.call("rbind", genome_list)
  
  genome_coverage = rowSums(genome_matrix)
  return(genome_coverage)
  
}

HiC_control_coverage_500000 = make_genome_wide_hic_coverage(HiC_matrices_path = "/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/H1_control_merged/500000/",
                                                            bin_size = 500000)


hic_margi_correlation <- function(HiC_coverage=HiC_control_coverage_500000,
                                  input_bedpe_margi,
                                  resolution,
                                  log_margi=F,
                                  fig_title){
  

  input_bedpe = data.frame(fread(input_bedpe_margi))
  input_data = input_bedpe[which(input_bedpe[,1] %in% hg38_chromosomes[1:23] & 
                                   input_bedpe[,4] %in% hg38_chromosomes[1:23]),1:6]
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
  
  cov_data_RNA = generate_karyoploteR_data(tags = list(Gr_RNA_end),
                                       genome_gr = Gr_hg38, 
                                       window_size = resolution,
                                       amplifier = c(1), 
                                       threshold = NA, 
                                       names = c("RNA_end"))
  rm(Gr_RNA_end)
  
  Gr_DNA_end = GRanges(
    seqnames = Rle(as.character(input_data[,4])),
    ranges = IRanges(input_data[,5], end = input_data[,6], names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  cov_data_DNA = generate_karyoploteR_data(tags = list(Gr_DNA_end),
                                       genome_gr = Gr_hg38, 
                                       window_size = resolution,
                                       amplifier = c(1), 
                                       threshold = NA, 
                                       names = c("DNA_end"))
  
  rm(Gr_DNA_end)
  rm(input_data)
  
  output = c()
  output = c(output, cor(HiC_coverage, cov_data_RNA$RNA_end, method = "pearson"))
  output = c(output, cor(HiC_coverage, cov_data_DNA$DNA_end, method = "pearson"))
  output = c(output, cor(HiC_coverage, cov_data_RNA$RNA_end, method = "spearman"))
  output = c(output, cor(HiC_coverage, cov_data_DNA$DNA_end, method = "spearman"))
  names(output) = c("PCC_RNA","PCC_DNA","SCC_RNA","SCC_DNA")
  
  df = data.frame(HiC_coverage, cov_data_RNA$RNA_end, cov_data_DNA$DNA_end)
  
  annotations_RNA <- data.frame(
    xpos = c(Inf),
    ypos =  c(-Inf),
    annotateText = c(paste0("PCC: ",round(output["PCC_RNA"],2),"\nSCC: ",round(output["SCC_RNA"],2))),
    hjustvar = c(1) ,
    vjustvar = c(-0.5))
  
  annotations_DNA <- data.frame(
    xpos = c(Inf),
    ypos =  c(-Inf),
    annotateText = c(paste0("PCC: ",round(output["PCC_DNA"],2),"\nSCC: ",round(output["SCC_DNA"],2))),
    hjustvar = c(1) ,
    vjustvar = c(-0.5))
  
  label_size = 12
  
  p1<-ggplot(df, aes(x=log(cov_data_RNA.RNA_end), y=HiC_coverage)) + 
    geom_point(color = "blue") +
    xlab("log(RNA end read coverage)") +
    ylab("HiC read coverage control") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size)) +
    ggtitle(fig_title) +
    geom_text(data=annotations_RNA,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  
  
  p2<-ggplot(df, aes(x=cov_data_DNA.DNA_end, y=HiC_coverage)) + 
    geom_point(color = "purple") +
    xlab("DNA end read coverage") +
    ylab("HiC read coverage control") +
    scale_x_continuous(labels = scientific, breaks = pretty_breaks(n=3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=label_size),
          axis.text.x = element_text(size=label_size),
          axis.text.y = element_text(size=label_size)) +
    ggtitle(fig_title) +
    geom_text(data=annotations_DNA,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  
  return(list(p1,p2))
  
}

temp = hic_margi_correlation(HiC_coverage=HiC_control_coverage_500000,
                            input_bedpe_margi="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe.gz",
                            resolution=500000,
                            fig_title="H1_control")
temp1 = hic_margi_correlation(HiC_coverage=HiC_control_coverage_500000,
                             input_bedpe_margi="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_FL/iMARGI_H1_FL.mapq30.1k.final.bedpe.gz",
                             resolution=500000,
                             fig_title="H1_FL")
temp2 = hic_margi_correlation(HiC_coverage=HiC_control_coverage_500000,
                              input_bedpe_margi="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_NH4OAc/iMARGI_H1_NH4OAc.mapq30.1k.final.bedpe.gz",
                              resolution=500000,
                              fig_title="H1_NH4OAc")
temp3 = hic_margi_correlation(HiC_coverage=HiC_control_coverage_500000,
                              input_bedpe_margi="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_RNase/iMARGI_H1_RNase.mapq30.1k.final.bedpe.gz",
                              resolution=500000,
                              fig_title="H1_RNase")

png("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/all_scatter_plots_logRNA.png", width = 15, height = 7, res = 200, units = "in")
plot_grid(plotlist=c(temp,temp1,temp2,temp3), nrow=2)
dev.off()



######################## Calculate euclidean distance between eigenvector and read coverage on the RNA and DNA ends
compartment_margi_euclidean_distance <- function(eigen_file,
                                          input_bedpe_file,
                                          resolution,
                                          log_margi=F,
                                          fig_title,
                                          force_compartment_control=F,
                                          force_imargi_control=F,
                                          MALAT1_cov=F){
  
  if (force_compartment_control == F){
    compartment = read.table(eigen_file, header = T)
    colnames(compartment) = c("chr","coord","eigen")
    compartment = compartment[which(compartment$chr != "chrY"),]
  } else {
    compartment = read.table("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_control_500000.txt", header = T)
    colnames(compartment) = c("chr","coord","eigen")
    compartment = compartment[which(compartment$chr != "chrY"),]
  }
  
  if (force_imargi_control == F){
    input_bedpe = data.frame(fread(input_bedpe_file))
  } else {
    input_bedpe = data.frame(fread("/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.final.bedpe.gz"))
  }
  
  input_data = input_bedpe[which(input_bedpe[,1] %in% hg38_chromosomes & 
                                   input_bedpe[,4] %in% hg38_chromosomes),1:6]
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
  
  if (MALAT1_cov == T){
    malat1 = annotation_hg38[which(annotation_hg38$gene_name=="MALAT1"),][1,]
    Gr_malat1 = GRanges(
      seqnames = Rle(as.character(malat1[1])),
      ranges = IRanges(as.numeric(malat1[2]), end = as.numeric(malat1[3]), names = 1),
      strand = Rle(strand('*')))
    
    overlap = countOverlaps(Gr_RNA_end, Gr_malat1, ignore.strand = T)
    mcols(Gr_RNA_end)["overlap_malat1"] = overlap
    Gr_RNA_end = Gr_RNA_end[mcols(Gr_RNA_end)[,"overlap_malat1"] >= 1]
    
    Gr_DNA_end = Gr_DNA_end[names(Gr_RNA_end)]
  }
  
  cov_data = generate_karyoploteR_data(tags = list(Gr_RNA_end,Gr_DNA_end),
                                       genome_gr = Gr_hg38, 
                                       window_size = resolution,
                                       amplifier = c(1,1), 
                                       threshold = NA, 
                                       names = c("RNA_end","DNA_end"))
  
  # rm(input_data)
  # rm(Gr_RNA_end)
  # rm(Gr_DNA_end)
  
  # if (log_margi == T){
  #   cov_data$RNA_end = log(cov_data$RNA_end)
  #   cov_data$RNA_end[which(is.infinite(cov_data$RNA_end))] = 0
  # 
  #   cov_data$DNA_end = log(cov_data$DNA_end)
  #   cov_data$DNA_end[which(is.infinite(cov_data$DNA_end))] = 0
  # }
  
  output = c()
  output = c(output, sqrt(sum((compartment$eigen - cov_data$RNA_end)^2)))
  output = c(output, sqrt(sum((compartment$eigen - cov_data$DNA_end)^2)))
  names(output) = c("euclidean_dist_RNA","euclidean_dist_DNA")
  
  return(output)
}


############# Merged replicates
euclidean_distance_list = list()

for (sample in all_imargi_samples){
  sample_name=gsub("iMARGI_","",sample)
  eigen_sample = strsplit(sample,"_")[[1]][3]
  euclidean_distance_list[[sample]] = compartment_margi_euclidean_distance(eigen_file=paste0("/dataOS/rcalandrelli/phase_separation/HiC/compartments/H1_", eigen_sample, "_500000.txt"),
                                                             input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/",sample,"/", sample, ".mapq30.1k.final.bedpe.gz"),
                                                             resolution=500000,
                                                             log_margi = F,
                                                             fig_title=sample_name,
                                                             force_compartment_control=F,
                                                             force_imargi_control=F)
}









