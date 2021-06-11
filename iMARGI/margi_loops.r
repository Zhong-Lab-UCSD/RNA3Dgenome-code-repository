############## MARGI analysis related to Hi-C loops
library(reshape2)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(pracma)
library(VennDiagram)

options(scipen=999)


############### Detect pattern 1 (high density regions) systematically using an "iMARGI insulation score" (NOT USED, see below)

calculate_insulation_score_and_HDR_margi <- function(contact_matrix_file,
                                             chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                             chromosome,
                                             sample,
                                             bin_size,
                                             template_dim,
                                             n_remove_diag=0, 
                                             non_zero_cutoff_score,
                                             bins_cutoff_HDR){
  
  print(paste0(sample,"_",chromosome))
  
  out_dir_score = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/insulation_score/", sample, "_", bin_size)
  if (!dir.exists(out_dir_score)){
    dir.create(out_dir_score)
  }
  
  out_dir_HDR = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/HDR/", sample, "_", bin_size)
  if (!dir.exists(out_dir_HDR)){
    dir.create(out_dir_HDR)
  }
  
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ##########  Tile genome and make genomic windows for selected chromosomes
  genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = bin_size, cut.last.tile.in.chrom = T)
  chr_genome_window = genome_window[seqnames(genome_window)==chromosome]

  ##########  Load contact matrix
  temp = data.frame(fread(contact_matrix_file))
  contact_matrix = matrix(0, nrow = length(chr_genome_window), ncol = length(chr_genome_window))
  contact_matrix[as.matrix(temp[,c("RNA_end_bins","DNA_end_bins")])] = temp$x

  chr = rep(chromosome, nrow(contact_matrix))
  bin = seq(1,nrow(contact_matrix))
  insulation_score_data = data.frame(chr, bin)
  insulation_score_data$insulation_score = 0
  
  ########## Calculate insulation score for each bin
  for (bin in 1:nrow(contact_matrix)){
    #print(paste0("Bin ", bin))
    shift = template_dim/2
    shift_bin = shift / bin_size
    upstream_bin = bin - shift_bin
    downstream_bin = bin + shift_bin - 1
    
    if (upstream_bin > 0 & downstream_bin <= chrSize[which(chrSize$V1 == chromosome),"V2"] %/% bin_size){
      
      L_matrix = contact_matrix[upstream_bin:bin, upstream_bin:bin]
      R_matrix = contact_matrix[bin:downstream_bin, bin:downstream_bin]
      X_up_matrix = contact_matrix[upstream_bin:bin, bin:downstream_bin]
      X_down_matrix = contact_matrix[bin:downstream_bin, upstream_bin:bin]
      
      non_zero_L_ratio = sum(L_matrix != 0) / shift_bin^2
      non_zero_R_ratio = sum(R_matrix != 0) / shift_bin^2
      
      ### Check that either L_matrix or R_matrix have a percentage of non-zero bins above non_zero_cutoff_score
      if (non_zero_L_ratio > non_zero_cutoff_score | non_zero_R_ratio > non_zero_cutoff_score){
        
        if (n_remove_diag == -1){
          mean_L = sum(L_matrix) / shift_bin^2
          mean_R = sum(R_matrix) / shift_bin^2
          mean_X_up = sum(X_up_matrix) / shift_bin^2
          mean_X_down = sum(X_down_matrix) / shift_bin^2
        } else {
          ### Removing elements on the diagonals
          num_L = sum(L_matrix)
          num_R = sum(R_matrix)
          denom = shift_bin^2
          for (i in 0:n_remove_diag){
            num_L = num_L - sum(Diag(L_matrix,i))
            num_R = num_R - sum(Diag(R_matrix,i))
            denom = denom - length(Diag(L_matrix,i))
          }
          mean_L = num_L / denom
          mean_R = num_R / denom
          mean_X_up = sum(X_up_matrix) / shift_bin^2
          mean_X_down = sum(X_down_matrix) / shift_bin^2
        }
        
        ### Insulation score calculation
        
        # if (non_zero_L_ratio > non_zero_cutoff & non_zero_R_ratio > non_zero_cutoff){
        #   insulation_score = (max(mean_L, mean_R)) / (max(mean_X_up, mean_X_down))
        # }
        # else if (non_zero_L_ratio > non_zero_cutoff & non_zero_R_ratio < non_zero_cutoff){
        #   insulation_score = mean_L / (max(mean_X_up, mean_X_down))
        # }
        # else if (non_zero_L_ratio < non_zero_cutoff & non_zero_R_ratio > non_zero_cutoff){
        #   insulation_score = mean_R / (max(mean_X_up, mean_X_down))
        # }
        
        insulation_score = (max(mean_L, mean_R)) / (max(mean_X_up, mean_X_down))
      
        
      } else {
        insulation_score = 0
      }
      
      
    } else {
      insulation_score = 0
    }
    
    insulation_score_data[bin,"insulation_score"] = insulation_score
    
  }
  
  write.table(insulation_score_data, paste0(out_dir_score, "/", chromosome, "_template_", template_dim, "_CS_", gsub("\\.","",as.character(non_zero_cutoff_score)), ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

  
  ################## Extract high density regions (HDR)
  
  peaks = findpeaks(insulation_score_data$insulation_score, nup = 1, minpeakdistance = 5, sortstr = F)
  insulation_score_data$peak = 0
  insulation_score_data[peaks[,2],"peak"] = 1
  
  ##### Plotting insulation score track with peaks
  
  df = insulation_score_data[325:375,]
  ggplot(df) +
    geom_point(aes(x=bin, y=insulation_score, color=peak), size = 3) +
    #scale_size(range=c(0,1.5)) +
    geom_line(aes(x=bin, y=insulation_score)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=24),
          axis.text.x = element_text(size = 24, color = "black"),
          axis.text.y = element_text(size = 24, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position="none")
  
  
  ##### Calculation of high density regions (HDR)
  
  peak_bins = unique(peaks[,2])
  peak_bins = peak_bins[order(peak_bins)]
  HDR_list = list()
  
  n = length(peak_bins)-1
  for (i in 1:n){
    start_bin = peak_bins[i]
    end_bin = peak_bins[i+1]
    temp_matrix = contact_matrix[start_bin:end_bin, start_bin:end_bin]
    k = sum(temp_matrix > 1) / dim(temp_matrix)[1]^2 # ratio of bins with > 1 interactions
    
    if (k > bins_cutoff_HDR){
      start = (start_bin * bin_size) - bin_size
      end = end_bin * bin_size
      out = c(chromosome, start, end, chromosome, start, end)
      HDR_list[[as.character(i)]] = out
    }
  }
  
  HDR = do.call("rbind", HDR_list)
  write.table(HDR, paste0(out_dir_HDR, "/", chromosome, "_template_", template_dim, "_CS_",gsub("\\.","",as.character(non_zero_cutoff_score)) , "_CHDR_" ,gsub("\\.","",as.character(bins_cutoff_HDR)) ,".bedpe"), row.names = F, col.names = F, quote = F, sep = "\t")
  
}



# my_file = "/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_control/10000_filter1k/chr1_chr1.txt"
# calculate_insulation_score_and_HDR_margi(contact_matrix_file=my_file,
#                                 chrSize_file="/dataOS/chenweizhong/Project____InterChrom/Data/4DN_files_reference/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#                                 chromosome="chr1",
#                                 sample="control",
#                                 bin_size=10000,
#                                 template_dim=500000,
#                                 n_remove_diag = 1,
#                                 non_zero_cutoff_score = 0.5,
#                                 bins_cutoff_HDR = 0.8)



for (my_sample in c("H1_control","H1_NH4OAc","H1_FL","H1_RNase","K562_control","HFF_control")){
  for (my_chr in hg38_chromosomes){
    my_file = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_",my_sample,"/10000_filter1k/", my_chr, "_", my_chr, ".txt")
    calculate_insulation_score_and_HDR_margi(contact_matrix_file=my_file,
                                             chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                             chromosome=my_chr,
                                             sample=my_sample,
                                             bin_size=10000,
                                             template_dim=300000,
                                             n_remove_diag = 1,
                                             non_zero_cutoff_score = 0.5,
                                             bins_cutoff_HDR = 0.7)
  }
}


###### Merge HDR together from all the chromosomes
HDR_list = list()
for (my_chr in hg38_chromosomes){
  temp = read.table(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/HDR/H1_control_10000/", my_chr, "_template_300000_CS_05_CHDR_07.bedpe"))
  HDR_list[[my_chr]] = temp
}

HDR_all_chromosomes = do.call("rbind", HDR_list)
write.table(HDR_all_chromosomes, "/dataOS/rcalandrelli/phase_separation/MARGI/HDR/H1_control_10000/all_chr_template_300000_CS_05_CHDR_07.bedpe", row.names = F, col.names = F, sep = "\t", quote = F)

Gr_HDR_all_chromosomes = GRanges(
  seqnames = Rle(HDR_all_chromosomes$V1),
  ranges = IRanges(HDR_all_chromosomes$V2+1, end = HDR_all_chromosomes$V3, names = c(1:nrow(HDR_all_chromosomes))),
  strand = Rle(strand('*')))

Gr_loop_control = GRanges(
  seqnames = Rle(paste0("chr",loops_control$V1)),
  ranges = IRanges(loops_control$V2+1, end = loops_control$V6, names = c(1:nrow(loops_control))),
  strand = Rle(strand('*')))

### Overlaps of HDR with loop regions
overlaps = countOverlaps(Gr_loop_control, Gr_HDR_all_chromosomes)
Gr_loop_control_HDR = Gr_loop_control[overlaps > 0]
Gr_loop_control_HDR = sort(Gr_loop_control_HDR)
names(Gr_loop_control_HDR) = seq(1:length(Gr_loop_control_HDR))

loop_control_HDR = data.frame(Gr_loop_control_HDR)
loop_control_HDR$seqnames = as.character(loop_control_HDR$seqnames)

add_ep_label_HDR_loops <- function(x, loop_ep_data)
{
  temp = loop_ep_data[which(loop_ep_data$V1 == gsub("chr","",as.character(x[1])) &
                              loop_ep_data$V2 == as.numeric(x[2]) - 1 &
                              loop_ep_data$V6 == as.numeric(x[3])),]
  if(nrow(temp)>0){
    return(temp$V25)
  } else {
    return("NA")
  } 
}


loop_control_HDR$EP = apply(loop_control_HDR, 1, add_ep_label_HDR_loops, loop_ep_data = loop_control_EP)
table(as.factor(loop_control_HDR$EP))

loop_control_HDR[which(loop_control_HDR$seqnames=="chr1" & loop_control_HDR$EP != "NA"),]

loop_control_HDR[which(loop_control_HDR$seqnames=="chr20" & 
                         loop_control_HDR$start >= 37150000 &
                         loop_control_HDR$end <= 37650000),]


HDR_all_chromosomes[which(HDR_all_chromosomes$V1=="chr20" & 
                         HDR_all_chromosomes$V2 >= 37150000 &
                         HDR_all_chromosomes$V3 <= 37650000),]





############### Detect pattern 1 (high density regions within loops) using a "in vs out" approach starting from loop coordinates (USE THIS!!!)

calculate_HDR_loop_margi <- function(loop_file,
                                     chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                     chromosomes,
                                     sample,
                                     bin_size,
                                     out_template_factor,
                                     shuffle_loops = F,
                                     out_dir_HDR=""
                                     )
  {
  
  if (out_dir_HDR == ""){
    out_dir_HDR = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_", sample)
  }
  
  chrSize = read.table(chrSize_file)
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ##########  Tile genome and make genomic windows for selected chromosomes
  genome_window <- tileGenome(seqinfo(Gr_genome), tilewidth = bin_size, cut.last.tile.in.chrom = T)
  
  loop_chr_list = list()
  for (my_chr in chromosomes){
    print(my_chr)
    chr_genome_window = genome_window[seqnames(genome_window)==my_chr]
    
    ##########  Load contact matrix
    temp = data.frame(fread(paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_",sample,"/10000_filter1k/", my_chr, "_", my_chr, ".txt")))
    contact_matrix = matrix(0, nrow = length(chr_genome_window), ncol = length(chr_genome_window))
    contact_matrix[as.matrix(temp[,c("RNA_end_bins","DNA_end_bins")])] = temp$x
    
    ########## Load loop data
    temp = read.table(loop_file)
    loops = temp[which(temp$V1 == gsub("chr","",my_chr)),]
    
    if (shuffle_loops == T){
      if (nrow(loops) > 0){
        for (l in 1:nrow(loops)){
          l_res = loops[l,3] - loops[l,2] # loop resolution
          l_size = loops[l,6] - loops[l,2] # loop_size 
          random_pos = runif(1, min=0, max=(chrSize[which(chrSize$V1==my_chr),2] - l_size)) # to avoid that the random loop goes beyond chromosome size
          loops[l,2] = random_pos
          loops[l,3] = random_pos + l_res
          loops[l,5] = random_pos + l_size - l_res
          loops[l,6] = random_pos + l_size
        }
      }
    }
    
    ########## Calculate density ratio for each bin
    if (nrow(loops) > 0){
      ratio_loops = list()
      for (i in 1:nrow(loops)){
        #print(paste0("Loop ", i))
        start_bin_loop = floor(loops[i,"V2"] / bin_size) + 1
        end_bin_loop = floor(loops[i,"V5"] / bin_size) + 1
        
        shift_bin = round((end_bin_loop - start_bin_loop + 1) * (out_template_factor - 1) / 2) 
        start_bin_out_template = start_bin_loop - shift_bin
        end_bin_out_template = end_bin_loop + shift_bin
        
        if (start_bin_out_template > 0 & end_bin_out_template <= chrSize[which(chrSize$V1 == my_chr),"V2"] %/% bin_size){
          
          ### Contacts within the loop region
          X_loop = contact_matrix[start_bin_loop:end_bin_loop, start_bin_loop:end_bin_loop]
          X_loop_mean = mean(X_loop)
          
          ### Contacts spanning outside of the loop region
          X_out_1 = contact_matrix[start_bin_out_template:(start_bin_loop - 1), start_bin_loop:end_bin_loop]
          X_out_2 = contact_matrix[start_bin_loop:end_bin_loop, (end_bin_loop + 1):end_bin_out_template]
          X_out_3 = contact_matrix[(end_bin_loop + 1):end_bin_out_template, start_bin_loop:end_bin_loop]
          X_out_4 = contact_matrix[start_bin_loop:end_bin_loop, start_bin_out_template:(start_bin_loop - 1)]
          
          X_out_mean = mean(c(X_out_1,X_out_2,X_out_3,X_out_4))
          
          in_out_ratio = ifelse(!is.na(X_loop_mean / X_out_mean), X_loop_mean / X_out_mean, 0)
          ratio_loops[[i]] = in_out_ratio
          
        } else {
          ratio_loops[[i]] = NaN
        } 
      }
      
      loops$density_ratio = do.call("rbind", ratio_loops)
      loop_chr_list[[my_chr]] = loops
      
    }
  }
  
  loops_out = do.call("rbind", loop_chr_list)
  if (shuffle_loops == F){
    write.table(loops_out, paste0(out_dir_HDR, "/loops_EP_ratio_template_factor_", out_template_factor,".bedpe"), row.names = F, col.names = F, sep= "\t", quote = F)
  } else {
    write.table(loops_out, paste0(out_dir_HDR, "/loops_EP_ratio_template_factor_", out_template_factor,"_shuffle.bedpe"), row.names = F, col.names = F, sep= "\t", quote = F)
  }
  return(loops_out)
}



loops_control_EP_ratio = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_control_EP.bedpe",
                                                  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                  chromosomes=hg38_chromosomes,
                                                  sample="H1_control",
                                                  bin_size=10000,
                                                  out_template_factor=2)

loops_NH4OAc_EP_ratio = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_NH4OAc_EP.bedpe",
                                                  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                  chromosomes=hg38_chromosomes,
                                                  sample="H1_NH4OAc",
                                                  bin_size=10000,
                                                  out_template_factor=2)

loops_FL_EP_ratio = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_FL_EP.bedpe",
                                                  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                  chromosomes=hg38_chromosomes,
                                                  sample="H1_FL",
                                                  bin_size=10000,
                                                  out_template_factor=2)

loops_RNase_EP_ratio = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_RNase_EP.bedpe",
                                                  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                  chromosomes=hg38_chromosomes,
                                                  sample="H1_RNase",
                                                  bin_size=10000,
                                                  out_template_factor=2)


loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio > 3 &
                               loops_control_EP_ratio$V1 == 10 & loops_control_EP_ratio$V2 > 69600000 & loops_control_EP_ratio$V6 < 70100000),]

loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio > 2.819 &
                               loops_control_EP_ratio$V1 == 1 & loops_control_EP_ratio$V2 > 3250000 & loops_control_EP_ratio$V6 < 3750000),]


loops_control_EP_ratio_shuffle = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_control_EP.bedpe",
                                                  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                  chromosomes=hg38_chromosomes[1:23],
                                                  sample="control",
                                                  bin_size=10000,
                                                  out_template_factor=2,
                                                  shuffle_loops = T)

loops_control_EP_ratio$density_ratio_shuffle = loops_control_EP_ratio_shuffle$density_ratio

wilcox.test(loops_control_EP_ratio$density_ratio, loops_control_EP_ratio$density_ratio_shuffle, alternative = "two.sided")

###### MicroC validation analysis
loops_MicroC_ratio = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_ES_MicroC/merged_loops.bedpe",
                                                chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                chromosomes=hg38_chromosomes,
                                                sample="H1_control",
                                                bin_size=10000,
                                                out_template_factor=2,
                                              out_dir_HDR="/dataOS/rcalandrelli/phase_separation/MARGI/loops/H1_MicroC")

loops_MicroC_ratio_shuffle = calculate_HDR_loop_margi(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_ES_MicroC/merged_loops.bedpe",
                                                          chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                          chromosomes=hg38_chromosomes,
                                                          sample="H1_control",
                                                          bin_size=10000,
                                                          out_template_factor=2,
                                                          shuffle_loops = T,
                                                          out_dir_HDR="/dataOS/rcalandrelli/phase_separation/MARGI/loops/H1_MicroC")

wilcox.test(loops_MicroC_ratio$density_ratio, loops_MicroC_ratio_shuffle$density_ratio, alternative = "two.sided")


### Distribution of density ratio
quantile(loops_control_EP_ratio$density_ratio, 0.85)
nrow(loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio > 3.7 &
                                    !is.infinite(loops_control_EP_ratio$density_ratio)),])

nrow(loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio > 3.7 &
                                    !is.infinite(loops_control_EP_ratio$density_ratio) &
                                    !is.na(loops_control_EP_ratio$V25)),])

nrow(loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio < 3.7 &
                                    !is.infinite(loops_control_EP_ratio$density_ratio) &
                                    !is.na(loops_control_EP_ratio$V25)),])

temp = as.table(matrix(c(71,239-71,491,1571-491), nrow=2, ncol=2, byrow = T))
chisq.test(temp)

nrow(loops_control_EP_ratio[which(loops_control_EP_ratio$density_ratio_shuffle > 3.7),])

p1<-ggplot(loops_control_EP_ratio, aes(x=density_ratio)) + 
  geom_histogram(color="black", fill="blue", bins=100) + 
  labs(x="MARGI density ratio", y="Number of loops") +
  scale_x_continuous(name="MARGI density ratio", breaks=c(0,5,10,15,20,25), labels=c("0","5","10","15","20","25")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none") +
  geom_vline(aes(xintercept=3.7),
             color="red", size=0.8)

##### Violin plots
temp_melt = reshape2::melt(loops_control_EP_ratio[,c("density_ratio", "density_ratio_shuffle")])
temp_melt$value = as.numeric(temp_melt$value)
temp_melt$variable = gsub("density_ratio_shuffle","Shuffle data",temp_melt$variable)
temp_melt$variable = gsub("density_ratio","Real data",temp_melt$variable)

p2<-ggplot(temp_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) +
  stat_compare_means(method = "wilcox", size = 6, label.y = 35) +
  labs(x="", y=expression("MARGI density ratio")) +
  ylim(c(0,36)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")

png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/MARGI_density_ratio_hist_box.png", width = 6, height = 8, res = 200, units = "in")
plot_grid(p1,p2,nrow=2, labels = c("E","F"), label_size = 18)
dev.off()





make_loop_EP_ratio_file_higlass <- function(loop_file_EP_ratio_path,
                                            threshold_ratio,
                                            out_template_factor){
  
  loop_file_EP_ratio = read.table(paste0(loop_file_EP_ratio_path, "/loops_EP_ratio_template_factor_",out_template_factor,".bedpe"))
  
  HDR = loop_file_EP_ratio[which(loop_file_EP_ratio$density_ratio > threshold_ratio),]
  HDR = HDR[,1:6]
  HDR[,1] = paste0("chr",HDR[,1])
  HDR[,4] = paste0("chr",HDR[,4])
  HDR[,3] = HDR[,6]
  HDR[,5] = HDR[,2]
  
  write.table(HDR, paste0(loop_file_EP_ratio_path, "/HDR_template_factor_",out_template_factor, "_thres_",threshold_ratio,".bedpe"), row.names = F, col.names = F, quote = F, sep="\t")
  
}                           


make_loop_EP_ratio_file_higlass(loop_file_EP_ratio_path="/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control",
                                            threshold_ratio=3,
                                            out_template_factor=2)


######### Merge the loop associated to pattern 1 in iMARGI with loop domains
make_loop_full <- function(loops_EP_ratio,
                           loop_domain,
                           loop_density_ratio_hic){
  temp_loop_margi = loops_EP_ratio[,c(1,2,3,4,5,6,25,26,27)]
  colnames(temp_loop_margi) = c("chr_1","start_1","end_1","chr_2","start_2","end_2","index","EP","density_ratio")
  
  temp_loop_hic = loop_domain[,c(1,2,3,4,5,6,25)]
  colnames(temp_loop_hic) = c("chr_1","start_1","end_1","chr_2","start_2","end_2","loop_domain")
  
  temp_loop_hic_2 = loop_density_ratio_hic[,c(1:6,26)]
  colnames(temp_loop_hic_2) = c("chr_1","start_1","end_1","chr_2","start_2","end_2","density_ratio_hic")
  
  loops_full = merge(temp_loop_margi, temp_loop_hic, by = c("chr_1","start_1","end_1","chr_2","start_2","end_2"))
  loops_full = merge(loops_full, temp_loop_hic_2, by = c("chr_1","start_1","end_1","chr_2","start_2","end_2"))
  return(loops_full)
}

loops_control_full = make_loop_full(loops_control_EP_ratio, loop_domain_control, loop_control_density_ratio_hic)
loops_control_full$density_ratio = as.numeric(loops_control_full$density_ratio)
loops_control_full$density_ratio_hic = as.numeric(loops_control_full$density_ratio_hic)

loops_NH4OAc_full = make_loop_full(loops_NH4OAc_EP_ratio, loop_domain_NH4OAc, loop_NH4OAc_density_ratio_hic)
loops_NH4OAc_full$density_ratio = as.numeric(loops_NH4OAc_full$density_ratio)
loops_NH4OAc_full$density_ratio_hic = as.numeric(loops_NH4OAc_full$density_ratio_hic)

loops_FL_full = make_loop_full(loops_FL_EP_ratio, loop_domain_FL, loop_FL_density_ratio_hic)
loops_FL_full$density_ratio = as.numeric(loops_FL_full$density_ratio)
loops_FL_full$density_ratio_hic = as.numeric(loops_FL_full$density_ratio_hic)

loops_RNase_full = make_loop_full(loops_RNase_EP_ratio, loop_domain_RNase, loop_RNase_density_ratio_hic)
loops_RNase_full$density_ratio = as.numeric(loops_RNase_full$density_ratio)
loops_RNase_full$density_ratio_hic = as.numeric(loops_RNase_full$density_ratio_hic)

loops_control_full = loops_control_full[which(!is.infinite(loops_control_full$density_ratio)),]
wilcox.test(loops_control_full$density_ratio_hic, loops_control_full$density_ratio, paired = T)

##### Scatter plot

thres_density_ratio = 3.7
thres_hic_density_ratio = 3.7
labelsize = 12
p1<-ggplot(loops_control_full, aes(x=density_ratio_hic, y=density_ratio)) +
  geom_point(aes(color=loop_domain)) +
  geom_smooth(method='lm', se = F, color="black") +
  geom_text(x=25, y=5, label="r = 0.21", size = 5) +
  scale_color_manual(labels = c("Loop domain", "Loop"), values = c("red", "blue")) +
  labs(x="HiC density ratio", y="MARGI density ratio") +
  xlim(0,30) + ylim(0,30) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=labelsize),
        axis.text.x = element_text(size = labelsize, color = "black"),
        axis.text.y = element_text(size = labelsize, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank(),
        legend.position = "top")


x = as.numeric(loops_control_full$density_ratio)
y = as.numeric(loops_control_full$density_ratio_hic)
cor(x,y, method = "pearson")

table(as.factor(loops_control_full[which(loops_control_full$density_ratio>3),"loop_domain"]))
hist(loops_control_full[which(loops_control_full$loop_domain=="loop_domain"),"density_ratio"])

##### Violin plots to compare HiC and MARGI density ratios
temp_melt = reshape2::melt(loops_control_full[,c("density_ratio", "density_ratio_hic")])
temp_melt$value = as.numeric(temp_melt$value)
temp_melt$variable = gsub("density_ratio_hic","HiC",temp_melt$variable)
temp_melt$variable = gsub("density_ratio","MARGI",temp_melt$variable)

p2<-ggplot(temp_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) +
  stat_compare_means(method = "wilcox", size = 6, label.y = 22, label = "p.signif") +
  ylim(c(0,23)) +
  labs(x="", y=expression("Density ratio")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=labelsize),
        axis.text.x = element_text(size = labelsize, color = "black"),
        axis.text.y = element_text(size = labelsize, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")

# Find p3 in script loop_analysis.r
png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/density_ratio_violin_HiC_MARGI.png", width = 10, height = 4, res = 200, units = "in")
plot_grid(p1,p2,p3, nrow=1, labels = c("A","B","C"), label_size = 18)
dev.off()



##### Venn plots
temp1 = rownames(loops_control_full[which(loops_control_full$density_ratio>thres_density_ratio),]) # Loops with high MARGI density ratio
temp2 = rownames(loops_control_full[which(loops_control_full$density_ratio_hic>thres_hic_density_ratio),]) # Loops with high HiC density ratio
temp3 = rownames(loops_control_full[which(loops_control_full$loop_domain=="loop_domain"),]) # Loop domains

### All together
venn.diagram(
  x = list(temp1, temp2, temp3),
  category.names = c("MARGI_HDR","HiC_HDR","LD"),
  filename = "/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/venn_diagram_3.png",
  output=T,
  
  # Main title
  main = "",
  main.fontfamily = "sans",
  main.fontface = "bold",
  main.cex = .8,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("red", "blue", "darkgreen"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.06, 0.06, 0.06),
  cat.fontfamily = "sans"
  #rotation = 1
)


make_venn <- function(x1, x2, categories, colors, cat_dist){
  venn.diagram(
    x = list(x1, x2),
    category.names = rev(categories),
    filename = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/", categories[1], "_", categories[2], "_venn.png"),
    output=T,
    
    # Main title
    main = "",
    main.fontfamily = "sans",
    main.fontface = "bold",
    main.cex = .8,
    
    # Output features
    imagetype="png" ,
    height = 500 , 
    width = 600 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = colors,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(130, -130),
    cat.dist = cat_dist,
    cat.fontfamily = "sans"
    #rotation = 1
  )
}

make_venn(temp1, temp2, c("MARGI_HDR","HiC_HDR"), c("red","blue"), c(0.22, 0.1))
make_venn(temp1, temp3, c("MARGI_HDR","LD"), c("red","darkgreen"), c(0.23, 0.08))
make_venn(temp2, temp3, c("HiC_HDR","LD"), c("blue","darkgreen"), c(0.16, 0.1))


figures_list = list()
figures_list[[1]] = readPNG("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/MARGI_HDR_HiC_HDR_venn.png")
figures_list[[2]] = readPNG("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/MARGI_HDR_LD_venn.png")
figures_list[[3]] = readPNG("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/HiC_HDR_LD_venn.png")
figures_list[[4]] = readPNG("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/venn_diagram_3.png")


pdf("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/venn_diagrams_separated.pdf", width = 10, height = 4)
plot_grid(rasterGrob(figures_list[[4]]),
          rasterGrob(figures_list[[1]]),
          rasterGrob(figures_list[[2]]), ncol=3, labels = c("D","E","F"), label_size = 18)
dev.off()


###### the high-HDR loops were associated with the high-MDR loops (odds ratio = xx, p-value = xx, Chi-square test) 
nrow(loops_control_full[which(loops_control_full$density_ratio>3.7),])
nrow(loops_control_full[which(loops_control_full$density_ratio_hic>3.7),])
nrow(loops_control_full[which(loops_control_full$density_ratio_hic<3.7 & loops_control_full$density_ratio<3.7),])
nrow(loops_control_full[which(loops_control_full$density_ratio_hic>3.7 & loops_control_full$density_ratio>3.7),])

temp = as.table(matrix(c(100,342,139,1229), nrow=2, ncol=2, byrow = T))
chisq.test(temp)


###### the loop domains were associated with the high-MDR loops (odds ratio = xx, p-value = xx, Chi-square test) 
nrow(loops_control_full[which(loops_control_full$density_ratio>3.7),])
nrow(loops_control_full[which(loops_control_full$loop_domain=="loop_domain"),])
nrow(loops_control_full[which(loops_control_full$density_ratio<3.7 & loops_control_full$loop_domain=="NA"),])
nrow(loops_control_full[which(loops_control_full$density_ratio>3.7 & loops_control_full$loop_domain=="loop_domain"),])

temp = as.table(matrix(c(106,514,133,1057), nrow=2, ncol=2, byrow = T))
chisq.test(temp)




################# Loop CTCF analysis

loops_ctcf_analysis <- function(loops_file,
                                ctcf_file_path){
  
  ### Loop data
  # loops_file = read.table(loops_file_path)[,1:6]
  # colnames(loops_file)[1] = "V1"
  loops_file[,1] = paste0("chr", loops_file[,1])
  loops_file[,4] = paste0("chr", loops_file[,4])
  
  Gr_loop_upstream <- GRanges(
    seqnames = Rle(loops_file$chr_1),
    ranges = IRanges(as.numeric(loops_file$start_1), end = as.numeric(loops_file$end_1), names = c(1:nrow(loops_file))),
    strand = Rle(strand('*')))
  
  Gr_loop_downstream <- GRanges(
    seqnames = Rle(loops_file$chr_2),
    ranges = IRanges(as.numeric(loops_file$start_2), end = as.numeric(loops_file$end_2), names = c(1:nrow(loops_file))),
    strand = Rle(strand('*')))
  
  ### CTCF data
  ctcf_file = data.frame(fread(ctcf_file_path))
  Gr_ctcf <- GRanges(
    seqnames = Rle(as.character(ctcf_file[,1])),
    ranges = IRanges(as.numeric(ctcf_file[,2]), end = as.numeric(ctcf_file[,3]), names = c(1:nrow(ctcf_file))),
    strand = Rle(strand('*')))
  
  ### Calculating overlap
  overlap_upstream <- countOverlaps(Gr_loop_upstream, Gr_ctcf)
  overlap_downstream <- countOverlaps(Gr_loop_downstream, Gr_ctcf)
  
  loops_file$ctcf_upstream <- overlap_upstream
  loops_file$ctcf_downstream <- overlap_downstream
  
  # total_loops = nrow(loops_file)
  # loops_both_end_ctcf <- nrow(loops_file[which(loops_file$ctcf_upstream > 0 & loops_file$ctcf_downstream > 0),])
  # loops_one_end_ctcf <- nrow(loops_file[which((loops_file$ctcf_upstream > 0 & loops_file$ctcf_downstream == 0) |
  #                                               (loops_file$ctcf_downstream > 0 & loops_file$ctcf_upstream == 0)),])
  # loops_no_ctcf <- nrow(loops_file[which(loops_file$ctcf_upstream == 0 & loops_file$ctcf_downstream == 0),])
  # 
  # out = c(total_loops, loops_both_end_ctcf, loops_one_end_ctcf, loops_no_ctcf)
  
  return(loops_file)
  
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

ctcf_file_path = "/dataOS/rcalandrelli/phase_separation/MARGI/loops/H1____GSM822297_hg38_wgEncodeOpenChromChipH1hescCtcfPk.bed"
ctcf_file_path = "/dataOS/wenxingzhao/database/ENCODE/H1_chipseq/MACS2_narrow_peaks/CTCF_MACS2.peak"

loops_control_full = loops_ctcf_analysis(loops_control_full, ctcf_file_path)
loops_control_full$ctcf_both = loops_control_full$ctcf_upstream + loops_control_full$ctcf_downstream

loops_NH4OAc_full = loops_ctcf_analysis(loops_NH4OAc_full, ctcf_file_path)
loops_NH4OAc_full$ctcf_both = loops_NH4OAc_full$ctcf_upstream + loops_NH4OAc_full$ctcf_downstream

loops_FL_full = loops_ctcf_analysis(loops_FL_full, ctcf_file_path)
loops_FL_full$ctcf_both = loops_FL_full$ctcf_upstream + loops_FL_full$ctcf_downstream

loops_RNase_full = loops_ctcf_analysis(loops_RNase_full, ctcf_file_path)
loops_RNase_full$ctcf_both = loops_RNase_full$ctcf_upstream + loops_RNase_full$ctcf_downstream


### HiC HDR
df_hic_hdr = loops_control_full[,c("density_ratio_hic","ctcf_upstream","ctcf_downstream","ctcf_both")]
df_hic_hdr$density_ratio_hic = ifelse(df_hic_hdr$density_ratio_hic > 3.7, "HiC_HDR","Not HiC_HDR")

df_hic_hdr_melt = reshape2::melt(df_hic_hdr)
df_hic_hdr_melt$title = "HiC density ratio"
colnames(df_hic_hdr_melt)[1] = "type"

p1<-ggplot(df_hic_hdr_melt, aes(x=variable, y=value, fill = type)) + 
  geom_split_violin() +
  stat_compare_means(method = "wilcox", size = 3, label.y = 15, label = "p.format") +
  labs(x="", y=expression("CTCF enrichment")) +
  scale_x_discrete(labels=c("ctcf_upstream" = "Up anchor", "ctcf_downstream" = "Down anchor", "ctcf_both" = "Both anchors")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank(),
        legend.position = "top")
  
  
### MARGI HDR
df_margi_hdr = loops_control_full[,c("density_ratio","ctcf_upstream","ctcf_downstream","ctcf_both")]
df_margi_hdr$density_ratio = ifelse(df_margi_hdr$density_ratio > 3.7, "MARGI_HDR","Not MARGI_HDR")

df_margi_hdr_melt = reshape2::melt(df_margi_hdr)
df_margi_hdr_melt$title = "MARGI density ratio"
colnames(df_margi_hdr_melt)[1] = "type"

p2<-ggplot(df_margi_hdr_melt, aes(x=variable, y=value, fill = type)) + 
  geom_split_violin() +
  stat_compare_means(method = "wilcox", size = 3, label.y = 15, label = "p.format") +
  labs(x="", y=expression("CTCF enrichment")) +
  scale_x_discrete(labels=c("ctcf_upstream" = "Up anchor", "ctcf_downstream" = "Down anchor", "ctcf_both" = "Both anchors")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank(),
        legend.position = "top")

### Loop domains
df_LD = loops_control_full[,c("loop_domain","ctcf_upstream","ctcf_downstream","ctcf_both")]
df_LD$loop_domain = ifelse(df_LD$loop_domain =="loop_domain", "Loop domain","Not loop domain")

df_LD_melt = reshape2::melt(df_LD)
df_LD_melt$title = "Loop domains"
colnames(df_LD_melt)[1] = "type"

p3<-ggplot(df_LD_melt, aes(x=variable, y=value, fill = type)) + 
  #facet_grid(cols=vars(title)) +
  geom_split_violin() +
  stat_compare_means(method = "wilcox", size = 3, label.y = 15, label = "p.format") +
  labs(x="", y=expression("CTCF enrichment")) +
  scale_x_discrete(labels=c("ctcf_upstream" = "Up anchor", "ctcf_downstream" = "Down anchor", "ctcf_both" = "Both anchors")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank(),
        legend.position = "top")


png("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/CTCF_violin_plots.png", width = 10, height = 4, res = 200, units = "in")
plot_grid(p1,p2,p3, ncol = 3, labels = c("G","",""), label_size = 18)
dev.off()



# loop_ctcf_summary_control = loops_ctcf_analysis(loops_file_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_control_merged/merged_loops.bedpe",
#                                                 ctcf_file_path)
# loop_ctcf_summary_NH4OAc = loops_ctcf_analysis(loops_file_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_NH4OAc_merged/merged_loops.bedpe",
#                                                ctcf_file_path)
# loop_ctcf_summary_FL = loops_ctcf_analysis(loops_file_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_FL_merged/merged_loops.bedpe",
#                                            ctcf_file_path)
# loop_ctcf_summary_RNase = loops_ctcf_analysis(loops_file_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_RNaseTreat_merged/merged_loops.bedpe",
#                                               ctcf_file_path)
# 
# df_loop_ctcf_summary = as.data.frame(rbind(loop_ctcf_summary_control,
#                                            loop_ctcf_summary_NH4OAc,
#                                            loop_ctcf_summary_FL,
#                                            loop_ctcf_summary_RNase))
# rownames(df_loop_ctcf_summary) = c("Control","NH4OAc","FL","RNase")
# colnames(df_loop_ctcf_summary) = c("Total_loops", "CTCF_both_sides", "CTCF_one_side", "CTCF_none")
# df_loop_ctcf_summary$sample = c("Control","NH4OAc","FL","RNase")
# df_loop_ctcf_summary[,2:4] = round(df_loop_ctcf_summary[,2:4] / df_loop_ctcf_summary[,1], 3)
# df_loop_ctcf_summary = df_loop_ctcf_summary[,-1]
# 
# df = melt(df_loop_ctcf_summary)
# 
# ggplot(df, aes(x=sample, y=value, fill=variable)) +
#   geom_bar(stat="identity", position=position_dodge())
# 
# 
# my_data = data.frame(
#   y=c(rnorm(1000), rnorm(1000, 0.5), rnorm(1000, 1), rnorm(1000, 1.5)),
#   x=c(rep('a', 2000), rep('b', 2000)),
#   m=c(rep('i', 1000), rep('j', 2000), rep('i', 1000))
# )






################# Loops associated with enriched iMARGI regions (computation in margi_loops_run_fun.r)

loops_label = read.table("/dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/loops_RNA_enrich_label.bedpe", header = T)
table(as.factor(loops_label$label))
loops_label$seqnames = as.character(loops_label$seqnames)

### Adding enhancer_promoter labeling
add_ep_label <- function(x,
                         loop_ep_data)
  {
  temp = loop_ep_data[which(loop_ep_data[,1] == gsub("chr","",as.character(x[1])) &
                              loop_ep_data[,2] == as.numeric(x[2]) &
                              loop_ep_data[,6] == as.numeric(x[3])),]
  if(nrow(temp)>0){
    return(temp$V25)
  } else {
    return("NA")
  } 
}

loops_label$EP = apply(loops_label, 1, add_ep_label, loop_ep_data = loop_control_EP)

enriched_loops = loops_label[which(loops_label$label=="enrichment"),]
table(as.factor(enriched_loops$EP))


# temp_list = list()
# for (i in 1:nrow(loops_label)){
#   temp1 = loop_ep_data[which(loop_ep_data$V1 == gsub("chr","",as.character(loops_label[i,1])) &
#                                loop_ep_data$V2 == as.numeric(loops_label[i,2]) &
#                                loop_ep_data$V6 == as.numeric(loops_label[i,3])),]
#   
#   if (nrow(temp1)>0){
#     temp_list[[i]] = as.character(temp1$V25)
#   } else {
#     temp_list[[i]] = "NA"
#   }
# }
# loops_label$EP = do.call("rbind", temp_list)



enriched_loops[which(enriched_loops$EP=="both"),]

loops_label[which(loops_label$seqnames=="chr10" & 
                    loops_label$start==73530000 &
                    loops_label$end == 73650000),]

loops_label[which(loops_label$seqnames=="chr1" & 
                    loops_label$start==204075000 &
                    loops_label$end == 204130000),]

loops_label[which(loops_label$seqnames=="chr3" & 
                    loops_label$start==38510000 &
                    loops_label$end == 38660000),]

enriched_loops[which(enriched_loops$seqnames=="chr10" & 
                       enriched_loops$start==73530000 &
                       enriched_loops$end == 73650000),]



##########################################
###### pattern 2 analysis

### Loading stripe data
H1_RNAstripes = read.table("/dataOS/rcalandrelli/phase_separation/MARGI/loops/RNAStripes_H1.bedpe.gz")

H1_RNAstripes[which(H1_RNAstripes$V1 == "chr20" &
                      H1_RNAstripes$V2 > 37000000 &
                      H1_RNAstripes$V3 < 38000000),]

Gr_H1_RNAstripes_RNA <- GRanges(
  seqnames = Rle(as.character(H1_RNAstripes$V1)),
  ranges = IRanges(as.numeric(H1_RNAstripes$V2) + 1, end = as.numeric(H1_RNAstripes$V3), names = c(1:nrow(H1_RNAstripes))),
  strand = Rle(strand('*')))

Gr_H1_RNAstripes_DNA <- GRanges(
  seqnames = Rle(as.character(H1_RNAstripes$V4)),
  ranges = IRanges(as.numeric(H1_RNAstripes$V5) + 1, end = as.numeric(H1_RNAstripes$V6), names = c(1:nrow(H1_RNAstripes))),
  strand = Rle(strand('*')))

### Make loop control dataset from union
loops_control_from_upset = df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 1),c(1,2,3,5,6,7)]
colnames(loops_control_from_upset) = c("V1","V2","V3","V4","V5","V6")
loops_control_from_upset$index = seq(1:nrow(loops_control_from_upset))

Gr_loops_control_from_upset_1 = GRanges(
  seqnames = Rle(loops_control_from_upset$V1),
  ranges = IRanges(loops_control_from_upset$V2, end = loops_control_from_upset$V3, names = c(1:nrow(loops_control_from_upset))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control_from_upset)))
Gr_loops_control_from_upset_2 = GRanges(
  seqnames = Rle(loops_control_from_upset$V1),
  ranges = IRanges(loops_control_from_upset$V5, end = loops_control_from_upset$V6, names = c(1:nrow(loops_control_from_upset))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control_from_upset)))

loops_control_from_upset_EP = find_promoter_enhancer_loop(Gr_loop_1=Gr_loops_control_from_upset_1, 
                                                          Gr_loop_2=Gr_loops_control_from_upset_2,
                                              loop_file=loops_control_from_upset,
                                              sample_name="H1_control_upset",
                                              Gr_enhancer=Gr_H1_enhancers,
                                              Gr_promoter=gencode.24_promoters)

### Analysis for upper triangular matrix
Gr_loop_control_RNA = GRanges(
  seqnames = Rle(loops_control_from_upset$V1),
  ranges = IRanges(loops_control_from_upset$V2, end = loops_control_from_upset$V3, names = c(1:nrow(loops_control_from_upset))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control_from_upset)))
Gr_loop_control_DNA = GRanges(
  seqnames = Rle(loops_control_from_upset$V1),
  ranges = IRanges(loops_control_from_upset$V5, end = loops_control_from_upset$V6, names = c(1:nrow(loops_control_from_upset))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control_from_upset)))

overlaps_RNA = data.frame(findOverlaps(Gr_H1_RNAstripes_RNA, Gr_loop_control_RNA))

out_list = list()
for (i in 1:nrow(overlaps_RNA)){
  Gr_loop_control_RNA_stripe_DNA = Gr_loop_control_DNA[overlaps_RNA[i,2]]
  overlaps_DNA = countOverlaps(Gr_H1_RNAstripes_DNA[overlaps_RNA[i,1]], Gr_loop_control_RNA_stripe_DNA)
  if (overlaps_DNA == 1){
    distance = end(ranges(Gr_loop_control_RNA_stripe_DNA)) - end(ranges(Gr_H1_RNAstripes_DNA[overlaps_RNA[i,1]]))
    out_list[[i]] = data.frame(stripe=overlaps_RNA[i,1], loop=overlaps_RNA[i,2], distance=distance)
  }
}

stripe_loop_df_upper = do.call("rbind", out_list)

### Analysis for lower triangular matrix
Gr_loop_control_DNA = GRanges(
  seqnames = Rle(loops_control_from_upset$V1),
  ranges = IRanges(loops_control_from_upset$V2, end = loops_control_from_upset$V3, names = c(1:nrow(loops_control_from_upset))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control_from_upset)))
Gr_loop_control_RNA = GRanges(
  seqnames = Rle(loops_control_from_upset$V1),
  ranges = IRanges(loops_control_from_upset$V5, end = loops_control_from_upset$V6, names = c(1:nrow(loops_control_from_upset))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control_from_upset)))

overlaps_RNA = data.frame(findOverlaps(Gr_H1_RNAstripes_RNA, Gr_loop_control_RNA))

out_list = list()
for (i in 1:nrow(overlaps_RNA)){
  Gr_loop_control_RNA_stripe_DNA = Gr_loop_control_DNA[overlaps_RNA[i,2]]
  overlaps_DNA = countOverlaps(Gr_H1_RNAstripes_DNA[overlaps_RNA[i,1]], Gr_loop_control_RNA_stripe_DNA)
  if (overlaps_DNA == 1){
    distance = start(ranges(Gr_loop_control_RNA_stripe_DNA)) - start(ranges(Gr_H1_RNAstripes_DNA[overlaps_RNA[i,1]]))
    out_list[[i]] = data.frame(stripe=overlaps_RNA[i,1], loop=overlaps_RNA[i,2], distance=distance)
  }
}

stripe_loop_df_lower = do.call("rbind", out_list)

stripe_loop_df_lower[which(stripe_loop_df_lower$distance == min(stripe_loop_df_lower$distance)),]

length(intersect(stripe_loop_df_upper$loop, stripe_loop_df_lower$loop))


###### Merge the data

stripe_loop_df = rbind(stripe_loop_df_upper, stripe_loop_df_lower)[,1:2]
stripe_loop_df = stripe_loop_df[!duplicated(stripe_loop_df),]

length(unique(stripe_loop_df$loop))
length(unique(stripe_loop_df$stripe))

temp = as.table(matrix(c(1067,2473-1067), nrow=1, ncol=2, byrow = T))
chisq.test(temp)

loops_control_stripe = loops_control_from_upset_EP[unique(stripe_loop_df$loop),]
#temp = merge(loops_control_stripe, loops_control_EP_ratio, by = c("V1","V2","V3","V4","V5","V6"), all.x = T)
table(as.factor(loops_control_stripe$V26))


loops_control_no_stripe = loops_control_from_upset_EP[setdiff(c(1:2473),unique(stripe_loop_df$loop)),]
#temp = merge(loops_control_no_stripe, loops_control_EP_ratio, by = c("V1","V2","V3","V4","V5","V6"), all.x = T)
table(as.factor(loops_control_no_stripe$V26))


temp = as.table(matrix(c(600,1067-600,693,1406-693), nrow=2, ncol=2, byrow = T))
chisq.test(temp)




  




