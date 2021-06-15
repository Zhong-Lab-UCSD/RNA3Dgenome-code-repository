# 1____ratio_matrix_figure.R
require(GenomicRanges)
require(tidyverse)
require(tibble)
require(data.table)
require(magrittr)
require(InteractionSet)
library(parallel)
require(karyoploteR)
require(regioneR)
require(GenomicFeatures)
`%ni%` = Negate(`%in%`)

# ---------- load data ----------

# granges only at gene level, genecode v24
gtfgr24 <- rtracklayer::import("/home/frankyan/research/refGenome/standard4DN/GRCh38/gencode.v24.primary_assembly.annotation.gtf")
gtfgr24 <- gtfgr24[gtfgr24$type=="gene"]


# prepare SPIN
H1_bed <- fread("/dataOS/wenxingzhao/database/4DN/SPIN/H1_new.SPIN.JAWG.25kb.9_state.bed")

STATES <- c("Speckle",  "Interior_Act1","Interior_Act2","Interior_Act3","Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina")

color_ <- c("#8b254a","#c14e4c","#ec7b57","#f2b579","#dbd291","#a8d29f", "#5fbba2","#7d9a98","#54508b")

H1_spin <- H1_bed %>% dplyr::filter(!grepl("NAN",V4)) %>%
  makeGRangesFromDataFrame(.,seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)


# exclude less than 10M away RD interactions 
ctrl_gi <- readRDS("/dataOS/wenxingzhao/project/13_rna-dna/RData/iMARGI_H1_control.mapq30.1k_all_gi_nochrM.rds") # Ginteraction obj for all valid RD pairs

pairdist_ <- pairdist(ctrl_gi)

ctrl_gi_10m <- ctrl_gi[c(which(pairdist_ >= 10000000), which(is.na(pairdist_)))] # keep interactions with pairdist>10mb or interchromosomal 


# ---------- generate gene-SPIN count table ----------


# gene-SPIN count table 

findoverlap_with_perc_overlap <- function(querygr, subgr, QueryMinOverlapPerc, StrandIgnore){

  # an improved version of findOverlaps, enable min percentage
  hits <- findOverlaps(querygr, subgr, ignore.strand=StrandIgnore )

  # intersect region
  overlaps <- IRanges::pintersect(querygr[queryHits(hits)], subgr[subjectHits(hits)])

  # minPercentage
  percentOverlap <- width(overlaps) / width(querygr[queryHits(hits)])

  hits <- hits[percentOverlap > QueryMinOverlapPerc]

  return(as.data.frame(hits))

}

# RNA end mapping
RNA_end_anno <- findoverlap_with_perc_overlap(querygr = anchors(ctrl_gi_10m)$first,
                                              subgr = gtfgr24,
                                              QueryMinOverlapPerc = 0.51,
                                              StrandIgnore = F) # iMARGI RNA end is strand specific 

# if one read overlapped with two gene, chose the one with smaller index
colnames(RNA_end_anno) <- c("gi_idx", "gtfidx")
RNA_end_anno <- RNA_end_anno %>% dplyr::group_by(gi_idx) %>% dplyr::summarise(gtfidx_min = min(gtfidx))


# DNA end mapping
DNA_end_anno <- findoverlap_with_perc_overlap(querygr = anchors(ctrl_gi_10m)$second,
                                              subgr = H1_spin,
                                              QueryMinOverlapPerc = 0.51,
                                              StrandIgnore = T) # iMARGI DNA end is not strand specific 

colnames(DNA_end_anno) <- c("gi_idx", "spinidx")

# merge two
merged <- data.frame(gi_idx = seq(1,length(ctrl_gi_10m))) %>% dplyr::left_join(RNA_end_anno) %>% dplyr::left_join(DNA_end_anno)

#
gene_SPIN <- data.frame(gene = gtfgr24_df$gene_name[merged$gtfidx_min],
                        spin = H1_spin$V4[merged$spinidx]) %>%
  dplyr::filter(!is.na(gene)) %>% dplyr::filter(!is.na(spin)) %>%
  dplyr::group_by(gene,spin) %>% dplyr::summarise(ct = n())   %>%
  tidyr::spread(.,key=spin,value=ct,fill=0)

count <- gene_SPIN %>% tibble::column_to_rownames("gene")

##### check outliers ######
# pca <- prcomp(count)
# plot(pca$x[,1],pca$x[,2],xlab = "PC1",ylab = "PC2")

count <- count[rownames(count) %ni% names(which(pca$x[,1]>200000)),]

thre=1
filter_count <- count[apply(count,1,function(x){min(x)>thre}),]

null_ct <- apply(filter_count,2,sum)

pval_ <- apply(filter_count,1,function(x){
  t <- cbind(b=null_ct-unlist(x), a=unlist(x)) %>% chisq.test()
  return(t$p.value)
})

padj <- p.adjust(pval_ , method = "BH")
sig_genes <- filter_count %>% rownames_to_column("gene") %>% dplyr::mutate(pval=padj) %>% dplyr::filter(pval<0.01) %>% column_to_rownames("gene")

ratiomtx <- sig_genes[,c(1:length(STATES))] %>% sweep(.,2,apply(.,2,sum),"/") %>% sweep(.,1,apply(.,1,sum),"/") %>% .[,STATES]


# plot heatmap 

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(as.matrix(xs), probs = seq(0, 1, length.out=n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(ratiomtx, n = 100)


require(ComplexHeatmap)
require(viridis)
colorr <- circlize::colorRamp2(mat_breaks,viridis(length(mat_breaks)))

png(filename="/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Paper_V6_code_0408_2021_H1_9/FigureEigen____C____RatioMtxHeatmap_no_gene.png",
    res = 600, units = "in",
    width = 4.5, height = 4.5)

Heatmap(ratiomtx,
        col = colorr,
        name = "Ratio",
        show_row_names = F,
        cluster_columns = F,
        row_dend_width = unit(20, "pt"),
        column_names_gp = gpar(fontsize = 10),
        column_names_side = "t",
        #right_annotation = anno_apex_lamina_genes,
        # column_names_rot = 45,
        #top_annotation = column_an,
        clustering_distance_rows = "pearson")

dev.off()



# ---------- lincRNA, exon/intron ratio project on the UMAP ----------

#UMAP embedding
embd <- xzwensc::Embedding(t(ratiomtx),method = "umap",n_neighours = 20,min_dist = 0.5,rdm = 1117)

# intron. exon 
    
    full_gtf24 <- makeTxDbFromGFF("/dataOS/wenxingzhao/database/Human/genome/gencode.v24.primary_assembly.annotation.gtf")
    
    exons <- exonsBy(full_gtf24,by="gene")
    
    exons <- GenomicRanges::reduce(exons) # reduce overlap within each transcript
    
    #calculate exon intron ratio 
    eiratio <- mclapply(exons, function(x) {
      
      # if the gene only has one exons
      if(length(x)==1){
        ie_ratio = Inf 
        return(ie_ratio)
      } else {
        exon_length <- sum(width(x))
        intro_length <- max(IRanges::end(x)) - min(IRanges::start(x)) - exon_length
        ie_ratio <- exon_length/intro_length
        
        return(ie_ratio)
      }
      
    }, mc.cores = 10)
    
    
    gene_eiratio <- data.frame(gene_id = names(eiratio), eiratio = sapply(eiratio,function(x){x[[1]]})) %>% 
      left_join(gtfdf) %>%
      dplyr::select(gene_id, gene_name, eiratio, gene_type) %>% 
      dplyr::filter(!is.infinite(eiratio))
    

    # two feature (e/i ratio and lincRNA) df 
    feature_df <- gene_eiratio %>% mutate(lincRNA = ifelse(gene_type=="lincRNA",1,0)) %>% dplyr::select(-gene_id, -gene_type) %>% 
      dplyr::group_by(gene_name) %>% dplyr::summarise(eiratio = mean(eiratio), lincRNA = mean(lincRNA)) %>% 
      column_to_rownames("gene_name") %>% mutate(eiratio = log2(eiratio+1))
    

# ---------- other gene features ----------
    
# gene repeat table on gene level

hg38_human_rpt <- fread("/dataOS/wenxingzhao/database/Repeat_masker_human/hg38.fa.out.bed", header = F) %>% # repeatmasker 
                  makeGRangesFromDataFrame(., seqnames.field = "V1", start.field = "V2", end.field = "V3", strand.field = "V6", keep.extra.columns = T)

gene_rpt_ovlp <- findOverlaps(gtfgr24, hg38_human_rpt, ignore.strand = F, minoverlap = 10)

REPEAT_STATES <- c('LINE/L1','SINE/Alu','LINE/L2','SINE/MIR','LTR/ERVL-MaLR','DNA/hAT-Charlie','DNA/TcMar-Tigger','LTR/ERV1','LTR/ERVL','Simple_repeat')
top10_repreat <- data.frame(gene_name = gtfgr24$gene_name[queryHits(gene_rpt_ovlp)], repeats = hg38_human_rpt$V4[subjectHits(gene_rpt_ovlp)]) %>%
                   group_by(gene_name, repeats) %>% summarise( sum_ = n() ) %>%
                   spread(., key=repeats, value = sum_, fill=0) %>% .[, c("gene_name",REPEAT_STATES)]


# gene length and gene type 
gen_len_type <- gtfdf %>% dplyr::select(gene_name,width,gene_type)


# 3'UTR and 5'UTR length 
gr_3UTR <- fread("/dataOS/wenxingzhao/database/Human/genome/genecode_v24_primary_split_bed/3UTRs.bed") %>% mutate(width = V3-V2) %>% 
  dplyr::group_by(V4) %>% dplyr::summarise(utr3_mean_width = mean(width)) %>% as.data.frame()

gr_5UTR <- fread("/dataOS/wenxingzhao/database/Human/genome/genecode_v24_primary_split_bed/5UTRs.bed") %>% mutate(width = V3-V2) %>% 
  dplyr::group_by(V4) %>% dplyr::summarise(utr5_mean_width = mean(width))%>% as.data.frame()

ensg_symbol <- 
  gtfdf %>% dplyr::select(gene_id,gene_name) %>% dplyr::mutate(gene_id = gsub("\\.\\d+$","",gene_id,perl = TRUE)) 

UTR3_symbol <- gr_3UTR %>% left_join(.,ensg_symbol,by = c("V4"="gene_id")) %>% dplyr::select(gene_name,utr3_mean_width) %>% 
  dplyr::group_by(gene_name) %>% dplyr::summarise(utr3_mean_width = mean(utr3_mean_width)) 

UTR5_symbol <- gr_5UTR %>% left_join(.,ensg_symbol,by = c("V4"="gene_id")) %>% dplyr::select(gene_name,utr5_mean_width) %>% 
  dplyr::group_by(gene_name) %>% dplyr::summarise(utr5_mean_width = mean(utr5_mean_width))

    

X_all_genes <- gen_len_type %>% 
  
  dplyr::left_join(UTR3_symbol,by=c("gene_name"="gene_name")) %>% 
  dplyr::distinct(gene_name,.keep_all=T) %>% 
  
  dplyr::left_join(UTR5_symbol,by=c("gene_name"="gene_name")) %>% 
  dplyr::distinct(gene_name,.keep_all=T) %>% 
  
  dplyr::left_join(gene_ieratio_twocol,by=c("gene_name"="gene_name")) %>% 
  dplyr::distinct(gene_name,.keep_all=T) %>% 
  
  dplyr::left_join(top10_repreat,by=c("gene_name"="gene_name")) %>% 
  dplyr::distinct(gene_name,.keep_all=T) %>% 
  
  dplyr::mutate(utr3_mean_width=replace_na(utr3_mean_width,0)) %>% 
  
  dplyr::mutate(utr5_mean_width=replace_na(utr5_mean_width,0)) 


# add DNA end SPIN info 
X <- DNA_SPIN_df %>% left_join(X_all_genes,by=c("gene"="gene_name")) %>% 
     dplyr::filter(ieratio>0) %>% 
     tibble::column_to_rownames("gene")


X$lincRNA <- factor(ifelse(X$gene_type=="lincRNA",1,0))


# for every gene feature we use RPKM stats
X_norm <- X %>% dplyr::mutate(
  `LINE/L1` =                1000*`LINE/L1`/width,
  `LINE/L2` =                1000*`LINE/L2`/width,
  `SINE/Alu` =               1000*`SINE/Alu`/width,
  `SINE/MIR` =               1000*`SINE/MIR`/width,
  `LTR/ERVL-MaLR` =          1000*`LTR/ERVL-MaLR`/width,
  `DNA/hAT-Charlie` =        1000*`DNA/hAT-Charlie`/width,
  `DNA/TcMar-Tigger` =       1000*`DNA/TcMar-Tigger`/width,
  `LTR/ERV1` =               1000*`LTR/ERV1`/width,
  `LTR/ERVL` =               1000*`LTR/ERVL`/width,
  `Simple_repeat` =          1000*`Simple_repeat`/width,
  width = width/1000, 
  utr3_mean_width = utr3_mean_width/1000,
  utr5_mean_width = utr5_mean_width/1000) %>% dplyr::mutate_at("dna_spin",factor)


# all numeric value log transform
X_numeric <- select_if(X_norm,is.numeric) %>% sapply(.,function(x){log2(x+1)}) %>% scale()

X_cat <- select_if(X_norm,is.factor)
X_scale <- cbind(X_numeric,X_cat) %>% 
  dplyr::rename(Gene_Length=width) %>% 
  dplyr::rename(`Intron/Exon_Ratio`=ieratio)

rownames(X_scale) <- rownames(X)


STATES <- c("Speckle",  "Interior_Act1","Interior_Act2","Interior_Act3","Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina")
INOUT_STATES <- c(rep("Interior",6),rep("Exterior",3))
RAW_BINARY_SPIN <- X_scale %>% mutate(bin_dna_spin = factor(plyr::mapvalues(dna_spin, STATES, INOUT_STATES),
                                                            levels=c("Interior","Exterior"))) %>% dplyr::select(-dna_spin)


RAW_BINARY_SPIN <- RAW_BINARY_SPIN %>% drop_na

BINARY_SPIN <- RAW_BINARY_SPIN %>% group_by(bin_dna_spin) %>% do(sample_n(.,table(RAW_BINARY_SPIN$bin_dna_spin)["Exterior"])) # balanced training set 


# ----------- GLM train and test ------------------
# we will use logistic regression to evaluate the performance and plot the PRC curve.
set.seed(1117)
require(ROCR)
acc_one <- function(train_X, test_X, idx){
  
  train_log <- glm(bin_dna_spin ~ `SINE/Alu` + `Intron/Exon_Ratio` + utr3_mean_width + lincRNA,
                   data=train_X, maxit=500, family = "binomial" )
  
  
  predict_ <- predict(train_log, test_X, type = 'response')
  
  pred <- prediction(predictions = predict_,
                     plyr::mapvalues(test_X$bin_dna_spin,c("Interior","Exterior"),c(0,1)))
  perf <- performance(pred,"tpr","fpr")
  print(performance(pred,"auc")@y.values)
  DF <- data.frame(tpr = perf@y.values %>% unlist, #add the initial point
                   fpr = perf@x.values %>% unlist,
                   batch = idx)
  return(DF)
}







