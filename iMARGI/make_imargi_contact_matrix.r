######### Script to make iMARGI contact matrix from BEDPE file
library(GenomicRanges)
library(GenomicAlignments)
library(data.table)
library(GenomicInteractions)
library(InteractionSet)
library(GenomicFeatures)
library(AnnotationDbi)

options(scipen=999)


################################### Function called internally by the function "contact_matrix_from_bedpe" (NOT TO BE DIRECTLY RUN)
make_imargi_contact_matrix <- function(chr_row,
                                       chr_col,
                                       input_bedpe,
                                       chrSize,
                                       output_path,
                                       resolution,
                                       dist_filter,
                                       dense_matrix,
                                       only_RNA_starting_from_gene,
                                       only_RNA_within_gene
                                       ){
  
  print(paste0("Making ", chr_row, "_", chr_col, " at ", resolution, " ..."))
  
  ##### Making Granges genome
  Gr_genome <- GRanges(
    seqnames = Rle(as.character(chrSize$V1)),
    ranges = IRanges(rep(1,nrow(chrSize)), end = as.numeric(chrSize$V2), names = c(1:nrow(chrSize))),
    strand = Rle(strand('*')))
  seqlengths(Gr_genome) <- as.numeric(chrSize$V2)
  
  ##### Reading input bedpe file and select only valid chromosomes
  input_data = input_bedpe[which(input_bedpe[,1] == chr_row & input_bedpe[,4] == chr_col),1:6]
  input_data[,2] = input_data[,2] + 1 # 1-based 
  input_data[,5] = input_data[,5] + 1 # 1-based 
  input_data[,1] = as.character(input_data[,1])
  input_data[,4] = as.character(input_data[,4])
  colnames(input_data) = paste0("V",seq(1,6))
  # Select only intra-chromosomal read pairs over "dist_filter" distance if parameter is declared
  # if (!is.na(dist_filter) & chr_row == chr_col){
  #   input_data = input_data[which((input_data[,1] == input_data[,4] &
  #                        abs(input_data[,2]-input_data[,5]) >= dist_filter) |
  #                          input_data[,1] != input_data[,4]),]
  # }

  ##### Make Granges for RNA and DNA ends
  Gr_RNA_end = GRanges(
    seqnames = Rle(input_data$V1),
    ranges = IRanges(input_data$V2, end = input_data$V3, names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  Gr_DNA_end = GRanges(
    seqnames = Rle(input_data$V4),
    ranges = IRanges(input_data$V5, end = input_data$V6, names = c(1:nrow(input_data))),
    strand = Rle(strand('*')))
  
  if (only_RNA_starting_from_gene == T){
    
    # Make Granges with start coordinate of each RNA end read
    Gr_RNA_end_start = GRanges(
      seqnames = Rle(input_data$V1),
      ranges = IRanges(input_data$V2, end = input_data$V2+1, names = c(1:nrow(input_data))),
      strand = Rle(strand('*')))
    
    # Overlap with genes
    overlaps = countOverlaps(Gr_RNA_end_start, gencode.24_genes, ignore.strand = T)
    
    Gr_RNA_end = Gr_RNA_end[overlaps > 0]
    Gr_DNA_end = Gr_DNA_end[overlaps > 0]
    
    names(Gr_RNA_end) = seq(1:length(Gr_RNA_end))
    names(Gr_DNA_end) = seq(1:length(Gr_DNA_end))
    
  }
  
  if (only_RNA_within_gene == T){
    
    # Overlap with genes
    overlaps = countOverlaps(Gr_RNA_end, gencode.24_genes, ignore.strand = T, type = "within")
    
    Gr_RNA_end = Gr_RNA_end[overlaps > 0]
    Gr_DNA_end = Gr_DNA_end[overlaps > 0]
    
    names(Gr_RNA_end) = seq(1:length(Gr_RNA_end))
    names(Gr_DNA_end) = seq(1:length(Gr_DNA_end))
    
  }
  
  ##### Make genomic windows at "resolution" and calculating overlap
  genome_window_row <- tileGenome(seqinfo(Gr_genome)[chr_row], tilewidth = resolution, cut.last.tile.in.chrom = T)
  genome_window_col <- tileGenome(seqinfo(Gr_genome)[chr_col], tilewidth = resolution, cut.last.tile.in.chrom = T)
  
  RNA_end_bins = data.frame(findOverlaps(Gr_RNA_end, genome_window_row))
  RNA_end_bins = RNA_end_bins[!duplicated(RNA_end_bins$queryHits),"subjectHits"]
  
  DNA_end_bins = data.frame(findOverlaps(Gr_DNA_end, genome_window_col))
  DNA_end_bins = DNA_end_bins[!duplicated(DNA_end_bins$queryHits),"subjectHits"]
  
  ##### Make sparse contact matrix
  temp_contact_matrix = data.frame(cbind(RNA_end_bins,DNA_end_bins,1))
  colnames(temp_contact_matrix)[3] = "counts"
  
  contact_matrix_sparse = aggregate(temp_contact_matrix$counts, 
                                    by=list(RNA_end_bins=temp_contact_matrix$RNA_end_bins, DNA_end_bins=temp_contact_matrix$DNA_end_bins), 
                                    FUN=sum)
  write.table(contact_matrix_sparse, paste0(output_path,"/",chr_row, "_", chr_col, ".txt"), row.names = F, col.names = T, sep = "\t", quote = F) 
  
  ##### Make dense matrix
  if (dense_matrix == T){
    contact_matrix_full = matrix(0, nrow = length(genome_window_row), ncol = length(genome_window_col))
    contact_matrix_full[as.matrix(contact_matrix_sparse[,c("RNA_end_bins","DNA_end_bins")])] = contact_matrix_sparse$x
    rownames(contact_matrix_full) = paste0(chr_row,":",start(ranges(genome_window_row)),"-",end(ranges(genome_window_row)))
    colnames(contact_matrix_full) = paste0(chr_col,":",start(ranges(genome_window_col)),"-",end(ranges(genome_window_col)))
    write.table(contact_matrix_full, paste0(output_path,"/",chr_row, "_", chr_col, "_dense.txt"), row.names = T, col.names = T, sep = "\t", quote = F) 
  }

}


################################### Function to be run which calls the function"make_imargi_contact_matrix" (TO BE RUN)
contact_matrix_from_bedpe <- function(input_bedpe_file, # input bedpe file as output from imargi pipeline (can be directly inputted in gz format)
                                      chrSize_file, # chrSize file with two tab separated column: chr, size
                                      output_path, # output path to save contact matrices (without trailing slash in the end)
                                      all_chr_pairs=F, # set to "T" to compute all the contact matrices
                                      chrs_row=c(), # to choose only specific chromosomes on the rows (note only explicit pairs with chrs_col following the order in this vector will be considered)
                                      chrs_col=c(), # to choose only specific chromosomes on the columns (note only explicit pairs with chrs_row following the order in this vector will be considered)
                                      resolution_vector=c(), # vector of resolutions or bin size
                                      dist_filter="", # in order to add the label of the distance filtering to the resolution folder
                                      dense_matrix=F, # set to "T" to save also the contact matrices in dense format
                                      only_RNA_starting_from_gene=F,
                                      only_RNA_within_gene=F
                                      ) # if true, uses only RNA-DNA pairs where the RNA end read starts within a gene
  {
  print(paste0("Loading ", input_bedpe_file, " ..."))
  input_bedpe = data.frame(fread(input_bedpe_file))
  print("Done!")
  chrSize = read.table(chrSize_file)
  
  for (resolution in resolution_vector){
    resolution_dir = paste0(output_path, "/", resolution, "_", dist_filter)
    if (!dir.exists(resolution_dir)){
      dir.create(resolution_dir)
    }
    if (all_chr_pairs == T){
      temp1 = t(combn(chrSize$V1, 2))
      temp2 = temp1[,c(2,1)]
      temp3 = cbind(as.character(chrSize$V1), as.character(chrSize$V1))
      
      for (i in 1:nrow(temp1)){
        make_imargi_contact_matrix(chr_row = as.character(temp1[i,][1]), 
                                   chr_col = as.character(temp1[i,][2]),
                                   input_bedpe,
                                   chrSize,
                                   resolution_dir,
                                   resolution,
                                   dist_filter,
                                   dense_matrix,
                                   only_RNA_starting_from_gene,
                                   only_RNA_within_gene)
      }
      
      for (i in 1:nrow(temp2)){
        make_imargi_contact_matrix(chr_row = as.character(temp2[i,][1]), 
                                   chr_col = as.character(temp2[i,][2]),
                                   input_bedpe,
                                   chrSize,
                                   resolution_dir,
                                   resolution,
                                   dist_filter,
                                   dense_matrix,
                                   only_RNA_starting_from_gene,
                                   only_RNA_within_gene)
      }
      
      for (i in 1:nrow(temp3)){
        make_imargi_contact_matrix(chr_row = as.character(temp3[i,][1]), 
                                   chr_col = as.character(temp3[i,][2]),
                                   input_bedpe,
                                   chrSize,
                                   resolution_dir,
                                   resolution,
                                   dist_filter,
                                   dense_matrix,
                                   only_RNA_starting_from_gene,
                                   only_RNA_within_gene)
      }
      
    } else {
      for (i in 1:length(chrs_row)){
        make_imargi_contact_matrix(chr_row = chrs_row[i], 
                                   chr_col = chrs_col[i],
                                   input_bedpe,
                                   chrSize,
                                   resolution_dir,
                                   resolution,
                                   dist_filter,
                                   dense_matrix,
                                   only_RNA_starting_from_gene,
                                   only_RNA_within_gene)
      }
    }
  }
}


#################################################
### Run the function for all the samples

#all_imargi_samples = list.dirs(path = "/dataOS/rcalandrelli/phase_separation/MARGI/data", full.names = F, recursive = F)
all_imargi_samples = list.dirs(path = "/dataOS/frankyan/phase_separation/std_results/merged_iMARGI", full.names = F, recursive = F)
all_imargi_bio_replicates = list.dirs(path = "/dataOS/frankyan/phase_separation/std_results/iMARGI", full.names = F, recursive = F)

all_resolutions = c(1000000, 500000, 100000, 50000, 10000)
hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY')) # UCSC

TxDb.Hsapiens.gencode.24 = makeTxDbFromGFF('/home/frankyan/research/refGenome/standard4DN/GRCh38/gencode.v24.primary_assembly.annotation.gtf')
gencode.24_genes = genes(TxDb.Hsapiens.gencode.24)
gencode.24_genes = gencode.24_genes[seqnames(gencode.24_genes) %in% hg38_chromosomes]
#gencode.24_promoters = promoters(gencode.24_genes, upstream = 2000, downstream = 2000, use.names = T)

for (i in all_imargi_samples[1:4]){
  sample_contact_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/", i)
  if (!dir.exists(sample_contact_dir)){
    dir.create(sample_contact_dir)
  }
  # contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.1k.bedpe.gz"),
  #                           chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  #                           output_path=sample_contact_dir,
  #                           all_chr_pairs=F,
  #                           chrs_row=hg38_chromosomes,
  #                           chrs_col=hg38_chromosomes,
  #                           resolution_vector=all_resolutions,
  #                           dist_filter="filter1k",
  #                           dense_matrix=F,
  #                           only_RNA_starting_from_gene=F)
  # 
  # contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.200k.bedpe.gz"),
  #                           chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  #                           output_path=sample_contact_dir,
  #                           all_chr_pairs=F,
  #                           chrs_row=hg38_chromosomes,
  #                           chrs_col=hg38_chromosomes,
  #                           resolution_vector=all_resolutions,
  #                           dist_filter="filter200k",
  #                           dense_matrix=F)
  # 
  # contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.1k_200k.bedpe.gz"),
  #                           chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  #                           output_path=sample_contact_dir,
  #                           all_chr_pairs=F,
  #                           chrs_row=hg38_chromosomes,
  #                           chrs_col=hg38_chromosomes,
  #                           resolution_vector=all_resolutions,
  #                           dist_filter="filter1k_200k",
  #                           dense_matrix=F)
  
  # contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.1k.bedpe.gz"),
  #                           chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
  #                           output_path=sample_contact_dir,
  #                           all_chr_pairs=F,
  #                           chrs_row=hg38_chromosomes,
  #                           chrs_col=hg38_chromosomes,
  #                           resolution_vector=all_resolutions,
  #                           dist_filter="filter1k_RNA_from_gene",
  #                           dense_matrix=F,
  #                           only_RNA_starting_from_gene=T)
  
  contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.1k.bedpe.gz"),
                            chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                            output_path=sample_contact_dir,
                            all_chr_pairs=F,
                            chrs_row=hg38_chromosomes,
                            chrs_col=hg38_chromosomes,
                            resolution_vector=all_resolutions,
                            dist_filter="filter1k_RNA_within_gene",
                            dense_matrix=F,
                            only_RNA_within_gene=T)
}

for (i in c("iMARGI_HFF_control")){
  sample_contact_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/", i)
  if (!dir.exists(sample_contact_dir)){
    dir.create(sample_contact_dir)
  }
  contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.200k.bedpe.gz"),
                            chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                            output_path=sample_contact_dir,
                            all_chr_pairs=F,
                            chrs_row=hg38_chromosomes,
                            chrs_col=hg38_chromosomes,
                            resolution_vector=all_resolutions,
                            dist_filter="filter200k",
                            dense_matrix=F)
}

# for (i in all_imargi_bio_replicates[2:11]){
#   sample_contact_dir = paste0("/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/", i)
#   if (!dir.exists(sample_contact_dir)){
#     dir.create(sample_contact_dir)
#   }
#   contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.1k.bedpe.gz"),
#                             chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#                             output_path=sample_contact_dir,
#                             all_chr_pairs=F,
#                             chrs_row=hg38_chromosomes,
#                             chrs_col=hg38_chromosomes,
#                             resolution_vector=all_resolutions,
#                             dist_filter="filter1k",
#                             dense_matrix=F)
# 
#   contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.200k.bedpe.gz"),
#                             chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#                             output_path=sample_contact_dir,
#                             all_chr_pairs=F,
#                             chrs_row=hg38_chromosomes,
#                             chrs_col=hg38_chromosomes,
#                             resolution_vector=all_resolutions,
#                             dist_filter="filter200k",
#                             dense_matrix=F)
#   
#   contact_matrix_from_bedpe(input_bedpe_file=paste0("/dataOS/rcalandrelli/phase_separation/MARGI/data/", i, "/", i, ".mapq30.1k_200k.bedpe.gz"),
#                             chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
#                             output_path=sample_contact_dir,
#                             all_chr_pairs=F,
#                             chrs_row=hg38_chromosomes,
#                             chrs_col=hg38_chromosomes,
#                             resolution_vector=all_resolutions,
#                             dist_filter="filter1k_200k",
#                             dense_matrix=F)
# }


###### chrM contact matrix
contact_matrix_from_bedpe(input_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_RNase/iMARGI_H1_RNase.mapq30.1k.bedpe.gz",
                          chrSize_file="/dataOS/rcalandrelli/mitoRNA/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                          output_path="/dataOS/rcalandrelli/phase_separation/MARGI/contact_data/iMARGI_H1_RNase/",
                          all_chr_pairs=F,
                          chrs_row=c("chrM","chr17"),
                          chrs_col=c("chr17","chrM"),
                          resolution_vector=c(1000,500,100,50,10),
                          dist_filter="filter1k",
                          dense_matrix=F)





#################################################
### Testing examples
# temp = makeGenomicInteractionsFromFile("/dataOS/frankyan/phase_separation/std_results/iMARGI/iMARGI_H1_control_4/strict_MAPQ30/iMARGI_H1_control_4.mapq30.bedpe.gz",
#                                        type="bedpe")
# 
# temp_file <- fread("/dataOS/frankyan/phase_separation/std_results/iMARGI/iMARGI_H1_control_4/strict_MAPQ30/iMARGI_H1_control_4.mapq30.bedpe.gz")



# input_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control_4/filter1k_iMARGI_H1_control_4.bedpe.gz"
# #input_bedpe_file="/dataOS/rcalandrelli/phase_separation/MARGI/data/iMARGI_H1_control_4/temp.bedpe"
# chrSize_file = "/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____InterChrom/Data/4DN_files_reference/4DNFI823LSII.hg38.mainonly.chrom.sizes"
# output_path = "/dataOS/rcalandrelli/phase_separation/MARGI"
# hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY')) # UCSC
# 
# input_bedpe = read.table("/dataOS/wenxingzhao/project/13_rna-dna/SPIN/1____method1_two_end_bed_intersect/anno_iMARGI_H1_RNaseTreat.bedpe.gz")
# chrSize = read.table(chrSize_file)
# # chr_row = "chr11"
# # chr_col = "chr11"
# # resolution = 1000000
# 
# make_imargi_contact_matrix(chr_row = "chr6", 
#                            chr_col = "chr6",
#                            input_bedpe,
#                            chrSize,
#                            output_path,
#                            resolution=2000,
#                            dist_filter=NA,
#                            dense_matrix=T)
# 
# contact_matrix_from_bedpe(input_bedpe_file=input_bedpe_file,
#                           chrSize_file=chrSize_file,
#                           output_path="/dataOS/rcalandrelli/phase_separation/MARGI/H1_4_10kb/",
#                           all_chr_pairs=F,
#                           chrs_row=hg38_chromosomes,
#                           chrs_col=hg38_chromosomes,
#                           resolution_vector=c(10000),
#                           dist_filter=NA,
#                           dense_matrix=T)



library(Rpairix)
filename = "/dataOS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.pairs.gz"
px_build_index(filename, preset = "pairs", force=TRUE)

px_seq1list("/dataOS/rcalandrelli/phase_separation/PLACseq/4DNFILZBUG96.pairs.gz") # list of first chromosomes

px_build_index("inst/test_4dn.pairs.gz", force=T)

px_chr1_col(filename)
px_chr2_col(filename)
px_startpos1_col(filename)
px_startpos2_col(filename)
px_endpos1_col(filename)
px_endpos2_col(filename)

filename="/home/rcalandrelli/R/x86_64-pc-linux-gnu-library/3.6/Rpairix/test_4dn.pairs.gz"

gr <- GRanges(
     seqnames = Rle(c("chr10", "chr20", "chr21", "chr22"), c(1, 2, 1, 2)),
     ranges = IRanges((0:5*1000000)+1, end = (0:5*1000000)+13000000))
grl <- split(gr, rep(1:2,3))
px_query(filename,query=grl)


