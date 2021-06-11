library(ggtern)
library(plotly)
library(readr)
library(dplyr)
library(tidyr)


###### Loop analysis RNA-genome interaction project

chrSize = read.table("/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes")
hg38_chromosomes = as.character(chrSize$V1)
hg38_lengths = as.numeric(chrSize$V2)
names(hg38_lengths) = hg38_chromosomes

loops_control = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_control_merged/merged_loops.bedpe")
loops_NH4OAc = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_NH4OAc_merged/merged_loops.bedpe")
loops_FL = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_FL_merged/merged_loops.bedpe")
loops_RNase = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_RNase_merged/merged_loops.bedpe")

############### Statistical difference loop numbers

require(dplyr)
tab_loop_control = loops_control %>% count(V1)
tab_loop_NH4OAc = loops_NH4OAc %>% count(V1)
tab_loop_FL = loops_FL %>% count(V1)
tab_loop_RNase = loops_RNase %>% count(V1)

t.test(tab_loop_control[which(tab_loop_control$V1 %in% tab_loop_NH4OAc$V1),"n"], tab_loop_NH4OAc$n, paired = T)
t.test(tab_loop_control[which(tab_loop_control$V1 %in% tab_loop_FL$V1),"n"], tab_loop_FL[which(tab_loop_FL$V1 != "Y"),"n"], paired = T)
t.test(tab_loop_control[which(tab_loop_control$V1 %in% tab_loop_RNase$V1),"n"], tab_loop_RNase[which(tab_loop_RNase$V1 != "Y"),"n"], paired = T)

############### Summary of distance between loop anchor points
distance_anchor_control = loops_control$V5 - loops_control$V2
distance_anchor_FL = loops_FL$V5 - loops_FL$V2
distance_anchor_NH4OAc = loops_NH4OAc$V5 - loops_NH4OAc$V2
distance_anchor_RNase = loops_RNase$V5 - loops_RNase$V2


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

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/anchor_points_distance_boxplots_log10.png", width = 6, height = 6, units = "in", res = 200)
ggplot(df, aes(x=sample, y=log10(distances))) + 
  geom_boxplot() +
  labs(x="", y="Log10 (Anchor points distance)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()

wilcox.test(distance_anchor_control,distance_anchor_FL)
wilcox.test(distance_anchor_control,distance_anchor_NH4OAc)
wilcox.test(distance_anchor_control,distance_anchor_RNase)

mean(distance_anchor_control)
mean(distance_anchor_FL)
mean(distance_anchor_NH4OAc)
mean(distance_anchor_RNase)

hist(distance_anchor_control)
hist(distance_anchor_FL)
hist(distance_anchor_NH4OAc)
hist(distance_anchor_RNase)


############### Aggregate peak analysis (APA)

## Summary statistics
make_summary_apa <- function(directory="/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/",
                 samples=c("control","NH4OAc","FL","RNase"),
                 resolution){
  sample_summary = list()
  for (i in samples){
    sample_summary[[i]] = read.table(paste0(directory,"H1_",i,"_merged/",resolution,"/gw/measures.txt"))
  }
  df_sample_summary = do.call("cbind", sample_summary)
  df_sample_summary = df_sample_summary[,c(1,2,4,6,8)]
  colnames(df_sample_summary) = c("Index","Control","NH4OAc","FL","RNase")
  df_sample_summary[,2:5] = round(df_sample_summary[,2:5],2)
  write.table(df_sample_summary,paste0(directory,"APA_statistics_summary_",resolution,".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  return(df_sample_summary)
}

summary_apa_5000 = melt(make_summary_apa(resolution="5000"))
summary_apa_10000 = melt(make_summary_apa(resolution="10000"))

df = summary_apa_5000[which(summary_apa_5000$Index == "P2LL"),]
df$Index = 5000
df = rbind(df, summary_apa_10000[which(summary_apa_10000$Index == "P2LL"),])
df[which(df$Index=="P2LL"),"Index"] = 10000
df$Index = factor(df$Index, c(5000,10000))
colnames(df)[1] = "Resolution"
df$variable = factor(df$variable, levels = c("Control", "FL", "NH4OAc", "RNase"))

p1<-ggplot(df, aes(fill=Resolution, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "P2LL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))

df = summary_apa_5000[which(summary_apa_5000$Index == "ZscoreLL"),]
df$Index = 5000
df = rbind(df, summary_apa_10000[which(summary_apa_10000$Index == "ZscoreLL"),])
df[which(df$Index=="ZscoreLL"),"Index"] = 10000
df$Index = factor(df$Index, c(5000,10000))
colnames(df)[1] = "Resolution"
df$variable = factor(df$variable, levels = c("Control", "FL", "NH4OAc", "RNase"))

p2<-ggplot(df, aes(fill=Resolution, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "ZscoreLL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/barplot_P2LL_ZscoreLL.png", width = 16, height = 8, units = "in", res = 200)
plot_grid(p1,p2)
dev.off()


## Plotting figures derived from Juicer apa
library(png)
library(grid)
library(gridExtra)

directory = "/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/"
samples = paste0("H1_",c("control","FL", "NH4OAc", "RNase"),"_merged")
resolutions = c("5000","10000")
figures_list = list()
k = 1
for (i in samples){
  for (j in resolutions){
    figures_list[[k]] = readPNG(paste0(directory,i,"/",j,"/gw/APA.png"))
    k = k + 1
  }
}

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/APA_global_plots.png", width = 5, height = 7, units = "in", res = 500)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]),
             rasterGrob(figures_list[[5]]),
             rasterGrob(figures_list[[6]]),
             rasterGrob(figures_list[[7]]),
             rasterGrob(figures_list[[8]]), ncol=2)
             # rasterGrob(figures_list[[9]]),
             # rasterGrob(figures_list[[10]]),
             # rasterGrob(figures_list[[11]]),
             # rasterGrob(figures_list[[12]]), ncol=3)
dev.off()

### Plot P2LL scores across all the chromosomes, for each sample, at the 3 resolutions

make_P2LL_ZscoreLL_df_plot <- function(directory="/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/",
                         samples = c("control","FL","NH4OAc","RNase"),
                         resolution){
  p2ll_vector = c()
  Zscore_vector = c()
  for (i in samples){
    for (j in gsub("chr","",hg38_chromosomes)[1:23]){
      infile = paste0(directory,"H1_",i,"_merged/",resolution,"/",j,"v",j,"/measures.txt")
      if (file.exists(infile)){
        p2ll_temp = read.table(infile)[4,2]
        Zscore_temp = read.table(infile)[6,2]
      } else {
        p2ll_temp = NA
        Zscore_temp = NA
      }
      p2ll_vector = c(p2ll_vector, p2ll_temp)
      Zscore_vector = c(Zscore_vector, Zscore_temp)
    }
  }
  sample = c(rep("Control",23),
             rep("NH4OAc",23),
             rep("FL",23),
             rep("RNase",23))
  #chromosomes = rep(hg38_chromosomes[1:23],4)
  chromosomes = rep(c(1:23),4)
  
  # P2LL data and plot
  df.p2ll = data.frame(sample,chromosomes,p2ll_vector)
  df.p2ll$p2ll_vector[which(is.na(df.p2ll$p2ll_vector))] = 0
  df.p2ll$sample = factor(df.p2ll$sample, c("Control","FL","NH4OAc","RNase"))
  
  p<-ggplot(df.p2ll, aes(x=chromosomes, y=p2ll_vector, color=sample)) + 
    geom_point() + geom_line() +
    labs(x = "", y = "P2LL") + 
    scale_x_discrete(breaks=c(1:23), limits = c(1:23), labels=hg38_chromosomes[1:23]) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20),
          axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.title = element_blank()) +
    ggtitle(paste0("Resolution: ", resolution))
  
  # ZscoreLL data and plot
  df.Zscore = data.frame(sample,chromosomes,Zscore_vector)
  df.Zscore$Zscore_vector[which(is.na(df.Zscore$Zscore_vector))] = 0
  df.Zscore$sample = factor(df.Zscore$sample, c("Control","FL","NH4OAc","RNase"))
  
  p1<-ggplot(df.Zscore, aes(x=chromosomes, y=Zscore_vector, color=sample)) + 
    geom_point() + geom_line() +
    labs(x = "", y = "ZscoreLL") + 
    scale_x_discrete(breaks=c(1:23), limits = c(1:23), labels=hg38_chromosomes[1:23]) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20),
          axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.title = element_blank()) +
    ggtitle(paste0("Resolution: ", resolution))
  
  return(list(df.p2ll,p,df.Zscore,p1))
  
}

P2LL_Zscore_list_5000 = make_P2LL_ZscoreLL_df_plot(resolution = "5000")
P2LL_Zscore_list_10000 = make_P2LL_ZscoreLL_df_plot(resolution = "10000")

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2LL_ZscoreLL_tracks.png", width = 8, height = 14, units = "in", res = 200)
plot_grid(P2LL_Zscore_list_5000[[2]],P2LL_Zscore_list_5000[[4]],
          P2LL_Zscore_list_10000[[2]],P2LL_Zscore_list_10000[[4]],
          nrow=4)
dev.off()

# Wilcoxon test 
P2LL_5000_df = P2LL_Zscore_list_5000[[1]]
wilcox.test(P2LL_5000_df[which(P2LL_5000_df$sample=="Control"),"p2ll_vector"][c(1:21,23)],P2LL_5000_df[which(P2LL_5000_df$sample=="NH4OAc"),"p2ll_vector"][c(1:21,23)])
wilcox.test(P2LL_5000_df[which(P2LL_5000_df$sample=="Control"),"p2ll_vector"][c(1:21,23)],P2LL_5000_df[which(P2LL_5000_df$sample=="FL"),"p2ll_vector"][c(1:21,23)])
wilcox.test(P2LL_5000_df[which(P2LL_5000_df$sample=="Control"),"p2ll_vector"][c(1:21,23)],P2LL_5000_df[which(P2LL_5000_df$sample=="RNase"),"p2ll_vector"][c(1:21,23)])

Zscore_5000_df = P2LL_Zscore_list_5000[[3]]
wilcox.test(Zscore_5000_df[which(Zscore_5000_df$sample=="Control"),"Zscore_vector"][c(1:21,23)],Zscore_5000_df[which(Zscore_5000_df$sample=="NH4OAc"),"Zscore_vector"][c(1:21,23)])
wilcox.test(Zscore_5000_df[which(Zscore_5000_df$sample=="Control"),"Zscore_vector"][c(1:21,23)],Zscore_5000_df[which(Zscore_5000_df$sample=="FL"),"Zscore_vector"][c(1:21,23)])
wilcox.test(Zscore_5000_df[which(Zscore_5000_df$sample=="Control"),"Zscore_vector"][c(1:21,23)],Zscore_5000_df[which(Zscore_5000_df$sample=="RNase"),"Zscore_vector"][c(1:21,23)])


P2LL_10000_df = P2LL_Zscore_list_10000[[1]]
wilcox.test(P2LL_10000_df[which(P2LL_10000_df$sample=="Control"),"p2ll_vector"][c(1:20,23)],P2LL_10000_df[which(P2LL_10000_df$sample=="NH4OAc"),"p2ll_vector"][c(1:20,23)])
wilcox.test(P2LL_10000_df[which(P2LL_10000_df$sample=="Control"),"p2ll_vector"][c(1:20,23)],P2LL_10000_df[which(P2LL_10000_df$sample=="FL"),"p2ll_vector"][c(1:20,23)])
wilcox.test(P2LL_10000_df[which(P2LL_10000_df$sample=="Control"),"p2ll_vector"][c(1:20,23)],P2LL_10000_df[which(P2LL_10000_df$sample=="RNase"),"p2ll_vector"][c(1:20,23)])

Zscore_10000_df = P2LL_Zscore_list_10000[[3]]
wilcox.test(Zscore_10000_df[which(Zscore_10000_df$sample=="Control"),"Zscore_vector"][c(1:20,23)],Zscore_10000_df[which(Zscore_10000_df$sample=="NH4OAc"),"Zscore_vector"][c(1:20,23)])
wilcox.test(Zscore_10000_df[which(Zscore_10000_df$sample=="Control"),"Zscore_vector"][c(1:20,23)],Zscore_10000_df[which(Zscore_10000_df$sample=="FL"),"Zscore_vector"][c(1:20,23)])
wilcox.test(Zscore_10000_df[which(Zscore_10000_df$sample=="Control"),"Zscore_vector"][c(1:20,23)],Zscore_10000_df[which(Zscore_10000_df$sample=="RNase"),"Zscore_vector"][c(1:20,23)])


##### P2M values for each loop: the ratio of the central pixel to the mean of the remaining pixels
p2m_control_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_control_merged/5000/gw/enhancement.txt")
p2m_NH4OAc_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_NH4OAc_merged/5000/gw/enhancement.txt")
p2m_FL_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_FL_merged/5000/gw/enhancement.txt")
p2m_RNase_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_RNase_merged/5000/gw/enhancement.txt")

wilcox.test(as.numeric(p2m_control_5000),as.numeric(p2m_NH4OAc_5000))
wilcox.test(as.numeric(p2m_control_5000),as.numeric(p2m_FL_5000))
wilcox.test(as.numeric(p2m_control_5000),as.numeric(p2m_RNase_5000))

p2m_control_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_control_merged/10000/gw/enhancement.txt")
p2m_NH4OAc_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_NH4OAc_merged/10000/gw/enhancement.txt")
p2m_FL_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_FL_merged/10000/gw/enhancement.txt")
p2m_RNase_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_RNase_merged/10000/gw/enhancement.txt")

wilcox.test(as.numeric(p2m_control_10000),as.numeric(p2m_NH4OAc_10000))
wilcox.test(as.numeric(p2m_control_10000),as.numeric(p2m_FL_10000))
wilcox.test(as.numeric(p2m_control_10000),as.numeric(p2m_RNase_10000))

# Distribution plots
df_p2m_5000 = data.frame(rbind(
  cbind("Control",t(p2m_control_5000)),
  cbind("NH4OAc",t(p2m_NH4OAc_5000)),
  cbind("FL",t(p2m_FL_5000)),
  cbind("RNase",t(p2m_RNase_5000))
  ))
df_p2m_5000$X1 = factor(df_p2m_5000$X1, levels = c("Control","FL","NH4OAc","RNase"))
df_p2m_5000$X2 = as.numeric(as.character(df_p2m_5000$X2))
colnames(df_p2m_5000) = c("sample","P2M")

df_p2m_10000 = data.frame(rbind(
  cbind("Control",t(p2m_control_10000)),
  cbind("NH4OAc",t(p2m_NH4OAc_10000)),
  cbind("FL",t(p2m_FL_10000)),
  cbind("RNase",t(p2m_RNase_10000))
))
df_p2m_10000$X1 = factor(df_p2m_10000$X1, levels = c("Control","FL","NH4OAc","RNase"))
df_p2m_10000$X2 = as.numeric(as.character(df_p2m_10000$X2))
colnames(df_p2m_10000) = c("sample","P2M")

df_p2m_10000$sample = factor(df_p2m_10000$sample, levels = c("Control","NH4OAc", "FL","RNase"))

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2M_boxplots.png", width = 12, height = 12, units = "in", res = 200)
p1<-ggplot(df_p2m_5000, aes(x=sample, y=P2M)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 5000")
p2<-ggplot(df_p2m_10000, aes(x=sample, y=P2M)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 10000")
p3<-ggplot(df_p2m_5000, aes(x=sample, y=P2M)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  ylim(0,20) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 5000")
p4<-ggplot(df_p2m_10000, aes(x=sample, y=P2M)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  ylim(0,20) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 10000")
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()


################################ Union of loops analysis
loops_control = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_control_merged/merged_loops.bedpe")
loops_NH4OAc = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_NH4OAc_merged/merged_loops.bedpe")
loops_FL = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_FL_merged/merged_loops.bedpe")
loops_RNase = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_RNase_merged/merged_loops.bedpe")

table_loops = as.data.frame(rbind(table(as.factor(loops_control$V3-loops_control$V2)),
                                  table(as.factor(loops_NH4OAc$V3-loops_NH4OAc$V2)),
                                  table(as.factor(loops_FL$V3-loops_FL$V2)),
                                  table(as.factor(loops_RNase$V3-loops_RNase$V2))))
table_loops$sample = c("Control","NH4OAc","FL","RNase")
df = melt(table_loops)
#df$sample = factor(df$sample, levels = c("Control","NH4OAc","FL","RNase"))
colnames(df)[2] = "Resolution"

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/barplot_loops.png", width = 8, height = 8, units = "in", res = 200)
ggplot(df, aes(fill=Resolution, y=value, x=sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "Number of loops") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()

### Union of loops
temp = rbind(loops_control,loops_NH4OAc,loops_FL,loops_RNase)
loops_union = temp[!duplicated(temp[,c(1,2,3,4,5,6)]),]
loops_union$V1 = as.character(loops_union$V1)
loops_union$V4 = as.character(loops_union$V4)

# Header to this file same as "merged_loops.bedpe" is manually added after saving
write.table(loops_union,"/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe", row.names = F, col.names = F, sep = "\t", quote = F)



##### Specific loop in Figure LOOPD B
loops_FL[which(loops_FL$V1 == "12" &
                 loops_FL$V2 == 53030000 &
                 loops_FL$V6 == 53480000),]

loops_RNase[which(loops_RNase$V1 == "12" &
                 loops_RNase$V2 == 53020000 &
                 loops_RNase$V6 == 53480000),]

temp = loops_union[which(loops_union$V1 == "12" &
                    loops_union$V2 %in% c(53020000,53030000) &
                    loops_union$V6 == 53480000),]

write.table(temp[,1:6],"/dataOS/rcalandrelli/phase_separation/HiC/loops/loops_LOOPD_B.bedpe", row.names = F, col.names = F, sep = "\t", quote = F)




out_list = list()
for (i in 1:nrow(loops_union)){
  out_list[[i]] = c(nrow(loops_control[do.call(paste, loops_control[1:6]) %in% paste(loops_union[i,1:6], collapse = " "), ]),
                  nrow(loops_NH4OAc[do.call(paste, loops_NH4OAc[1:6]) %in% paste(loops_union[i,1:6], collapse = " "), ]),
                  nrow(loops_FL[do.call(paste, loops_FL[1:6]) %in% paste(loops_union[i,1:6], collapse = " "), ]),
                  nrow(loops_RNase[do.call(paste, loops_RNase[1:6]) %in% paste(loops_union[i,1:6], collapse = " "), ]))
}

### Upset plot of union of loops
df_loops_union_upset = do.call("rbind", out_list)
colnames(df_loops_union_upset) = c("Control","NH4OAc","FL","RNase")
df_loops_union_upset = data.frame(df_loops_union_upset)

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/union_loops_upset_plot.png", width = 9, height = 5, units = "in", res = 200)
upset(df_loops_union_upset, 
      sets = c("RNase","NH4OAc","FL","Control"),
      sets.x.label = "Number of loops",
      order.by = "freq",
      keep.order = T,
      text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
      point.size = 3, line.size = 1)
dev.off()

### APA summaries
summary_apa_5000_union = melt(make_summary_apa(directory="/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/",
                samples=c("control","NH4OAc","FL","RNase"),
                resolution="5000"))
summary_apa_10000_union = melt(make_summary_apa(directory="/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/",
                 samples=c("control","NH4OAc","FL","RNase"),
                 resolution="10000"))

df = summary_apa_5000_union[which(summary_apa_5000_union$Index == "P2LL"),]
df$Index = 5000
df = rbind(df, summary_apa_10000_union[which(summary_apa_10000_union$Index == "P2LL"),])
df[which(df$Index=="P2LL"),"Index"] = 10000
df$Index = factor(df$Index, c(5000,10000))
colnames(df)[1] = "Resolution"
df$variable = factor(df$variable, levels = c("Control", "FL", "NH4OAc", "RNase"))

p1<-ggplot(df, aes(fill=Resolution, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "P2LL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))

df = summary_apa_5000_union[which(summary_apa_5000_union$Index == "ZscoreLL"),]
df$Index = 5000
df = rbind(df, summary_apa_10000_union[which(summary_apa_10000_union$Index == "ZscoreLL"),])
df[which(df$Index=="ZscoreLL"),"Index"] = 10000
df$Index = factor(df$Index, c(5000,10000))
colnames(df)[1] = "Resolution"
df$variable = factor(df$variable, levels = c("Control", "FL", "NH4OAc", "RNase"))

p2<-ggplot(df, aes(fill=Resolution, y=value, x=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "", y = "ZscoreLL") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/barplot_P2LL_ZscoreLL_union.png", width = 16, height = 8, units = "in", res = 200)
plot_grid(p1,p2)
dev.off()

## Plotting figures derived from Juicer apa
library(png)
library(grid)
library(gridExtra)

directory = "/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/"
samples = paste0("H1_",c("control","FL", "NH4OAc", "RNase"),"_merged")
resolutions = c("5000","10000")
figures_list = list()
k = 1
for (i in samples){
  for (j in resolutions){
    figures_list[[k]] = readPNG(paste0(directory,i,"/",j,"/gw/APA.png"))
    k = k + 1
  }
}

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/APA_union_global_plots.png", width = 5, height = 7, units = "in", res = 500)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]),
             rasterGrob(figures_list[[4]]),
             rasterGrob(figures_list[[5]]),
             rasterGrob(figures_list[[6]]),
             rasterGrob(figures_list[[7]]),
             rasterGrob(figures_list[[8]]), ncol=2)
dev.off()

### P2LL and ZscoreLL analysis

P2LL_Zscore_union_list_5000 = make_P2LL_ZscoreLL_df_plot(directory="/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/",
                                                   samples = c("control","NH4OAc","FL","RNase"),
                                                   resolution="5000")
P2LL_Zscore_union_list_10000 = make_P2LL_ZscoreLL_df_plot(directory="/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/",
                                                    samples = c("control","NH4OAc","FL","RNase"),
                                                    resolution="10000")

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2LL_union_ZscoreLL_tracks.png", width = 8, height = 14, units = "in", res = 200)
plot_grid(P2LL_Zscore_union_list_5000[[2]],P2LL_Zscore_union_list_5000[[4]],
          P2LL_Zscore_union_list_10000[[2]],P2LL_Zscore_union_list_10000[[4]],
          nrow=4)
dev.off()

# Wilcoxon test 
P2LL_union_5000_df = P2LL_Zscore_union_list_5000[[1]]
wilcox.test(P2LL_union_5000_df[which(P2LL_union_5000_df$sample=="Control"),"p2ll_vector"][c(1:21,23)],P2LL_union_5000_df[which(P2LL_union_5000_df$sample=="NH4OAc"),"p2ll_vector"][c(1:21,23)])
wilcox.test(P2LL_union_5000_df[which(P2LL_union_5000_df$sample=="Control"),"p2ll_vector"][c(1:21,23)],P2LL_union_5000_df[which(P2LL_union_5000_df$sample=="FL"),"p2ll_vector"][c(1:21,23)])
wilcox.test(P2LL_union_5000_df[which(P2LL_union_5000_df$sample=="Control"),"p2ll_vector"][c(1:21,23)],P2LL_union_5000_df[which(P2LL_union_5000_df$sample=="RNase"),"p2ll_vector"][c(1:21,23)])

Zscore_union_5000_df = P2LL_Zscore_union_list_5000[[3]]
wilcox.test(Zscore_union_5000_df[which(Zscore_union_5000_df$sample=="Control"),"Zscore_vector"][c(1:21,23)],Zscore_union_5000_df[which(Zscore_union_5000_df$sample=="NH4OAc"),"Zscore_vector"][c(1:21,23)])
wilcox.test(Zscore_union_5000_df[which(Zscore_union_5000_df$sample=="Control"),"Zscore_vector"][c(1:21,23)],Zscore_union_5000_df[which(Zscore_union_5000_df$sample=="FL"),"Zscore_vector"][c(1:21,23)])
wilcox.test(Zscore_union_5000_df[which(Zscore_union_5000_df$sample=="Control"),"Zscore_vector"][c(1:21,23)],Zscore_union_5000_df[which(Zscore_union_5000_df$sample=="RNase"),"Zscore_vector"][c(1:21,23)])

P2LL_union_10000_df = P2LL_Zscore_union_list_10000[[1]]
wilcox.test(P2LL_union_10000_df[which(P2LL_union_10000_df$sample=="Control"),"p2ll_vector"][c(1:20,23)],P2LL_union_10000_df[which(P2LL_union_10000_df$sample=="NH4OAc"),"p2ll_vector"][c(1:20,23)])
wilcox.test(P2LL_union_10000_df[which(P2LL_union_10000_df$sample=="Control"),"p2ll_vector"][c(1:20,23)],P2LL_union_10000_df[which(P2LL_union_10000_df$sample=="FL"),"p2ll_vector"][c(1:20,23)])
wilcox.test(P2LL_union_10000_df[which(P2LL_union_10000_df$sample=="Control"),"p2ll_vector"][c(1:20,23)],P2LL_union_10000_df[which(P2LL_union_10000_df$sample=="RNase"),"p2ll_vector"][c(1:20,23)])

Zscore_union_10000_df = P2LL_Zscore_union_list_10000[[3]]
wilcox.test(Zscore_union_10000_df[which(Zscore_union_10000_df$sample=="Control"),"Zscore_vector"][c(1:20,23)],Zscore_union_10000_df[which(Zscore_union_10000_df$sample=="NH4OAc"),"Zscore_vector"][c(1:20,23)])
wilcox.test(Zscore_union_10000_df[which(Zscore_union_10000_df$sample=="Control"),"Zscore_vector"][c(1:20,23)],Zscore_union_10000_df[which(Zscore_union_10000_df$sample=="FL"),"Zscore_vector"][c(1:20,23)])
wilcox.test(Zscore_union_10000_df[which(Zscore_union_10000_df$sample=="Control"),"Zscore_vector"][c(1:20,23)],Zscore_union_10000_df[which(Zscore_union_10000_df$sample=="RNase"),"Zscore_vector"][c(1:20,23)])


##### P2M values for each loop: the ratio of the central pixel to the mean of the remaining pixels

### 5 kb
p2m_union_control_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_control_merged/5000/gw/enhancement.txt")
p2m_union_NH4OAc_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_NH4OAc_merged/5000/gw/enhancement.txt")
p2m_union_FL_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_FL_merged/5000/gw/enhancement.txt")
p2m_union_RNase_5000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_RNase_merged/5000/gw/enhancement.txt")

# Remove NA entries
na.entries = unique(c(which(is.na(p2m_union_control_5000[1,])),
                      which(is.na(p2m_union_NH4OAc_5000[1,])),
                      which(is.na(p2m_union_FL_5000[1,])),
                      which(is.na(p2m_union_RNase_5000[1,]))))

p2m_union_control_5000 = p2m_union_control_5000[1,-na.entries]
p2m_union_NH4OAc_5000 = p2m_union_NH4OAc_5000[1,-na.entries]
p2m_union_FL_5000 = p2m_union_FL_5000[1,-na.entries]
p2m_union_RNase_5000 = p2m_union_RNase_5000[1,-na.entries]

wilcox.test(as.numeric(p2m_union_control_5000),as.numeric(p2m_union_NH4OAc_5000))
wilcox.test(as.numeric(p2m_union_control_5000),as.numeric(p2m_union_FL_5000))
wilcox.test(as.numeric(p2m_union_control_5000),as.numeric(p2m_union_RNase_5000))

### 10 kb
p2m_union_control_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_control_merged/10000/gw/enhancement.txt")
p2m_union_NH4OAc_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_NH4OAc_merged/10000/gw/enhancement.txt")
p2m_union_FL_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_FL_merged/10000/gw/enhancement.txt")
p2m_union_RNase_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_RNase_merged/10000/gw/enhancement.txt")

p2m_union_control_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_control_merged/10000/gw/enhancement.txt")
p2m_union_NH4OAc_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_NH4OAc_merged/10000/gw/enhancement.txt")
p2m_union_FL_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_FL_merged/10000/gw/enhancement.txt")
p2m_union_RNase_10000 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_RNase_merged/10000/gw/enhancement.txt")

# Remove NA entries
na.entries = unique(c(which(is.na(p2m_union_control_10000[1,])),
                      which(is.na(p2m_union_NH4OAc_10000[1,])),
                      which(is.na(p2m_union_FL_10000[1,])),
                      which(is.na(p2m_union_RNase_10000[1,]))))

p2m_union_control_10000 = p2m_union_control_10000[1,-na.entries]
p2m_union_NH4OAc_10000 = p2m_union_NH4OAc_10000[1,-na.entries]
p2m_union_FL_10000 = p2m_union_FL_10000[1,-na.entries]
p2m_union_RNase_10000 = p2m_union_RNase_10000[1,-na.entries]

wilcox.test(as.numeric(p2m_union_control_10000),as.numeric(p2m_union_NH4OAc_10000))
wilcox.test(as.numeric(p2m_union_control_10000),as.numeric(p2m_union_FL_10000))
wilcox.test(as.numeric(p2m_union_control_10000),as.numeric(p2m_union_RNase_10000))

# Distribution plots
df_p2m_union_5000 = data.frame(Control = t(p2m_union_control_5000),
                               FL = t(p2m_union_FL_5000),
                               NH4OAc = t(p2m_union_NH4OAc_5000),
                               RNase = t(p2m_union_RNase_5000))
colnames(df_p2m_union_5000) = c("Control","FL","NH4OAc","RNase")

df_p2m_union_10000 = data.frame(Control = t(p2m_union_control_10000),
                                FL = t(p2m_union_FL_10000),
                               NH4OAc = t(p2m_union_NH4OAc_10000),
                               RNase = t(p2m_union_RNase_10000))
colnames(df_p2m_union_10000) = c("Control","FL","NH4OAc","RNase")
df_10000 = melt(df_p2m_union_10000)
df_10000$variable = factor(df_10000$variable, levels = c("Control", "NH4OAc", "FL", "RNase"))

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2M_union_boxplots.png", width = 12, height = 12, units = "in", res = 200)
p1<-ggplot(melt(df_p2m_union_5000), aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 5000")
p2<-ggplot(melt(df_p2m_union_10000), aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 10000")
p3<-ggplot(melt(df_p2m_union_5000), aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(x="", y="P2M") +
  ylim(0,20) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 5000")

my_comparisons <- list(c("Control","NH4OAc"), c("Control","FL"), c("Control","RNase"))
p4<-ggplot(df_10000[which(df_10000$value <= 20),], aes(x=variable, y=value)) + 
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.x = c(2,3,4), label.y = c(20,22,24), size = 6) +
  labs(x="", y="P2M") +
  ylim(0,25) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.ticks.length = unit(0.2, "cm")) +
  ggtitle("Resolution: 10000")
plot_grid(p1,p2,p3,p4,ncol=2)
dev.off()


#### Scatter plots of the P2M values for each samples vs control

make_scatter_P2M <- function(df_p2m_union,
                             sample1="Control",
                             sample2,
                             resolution,
                             log_log=F){
  
  df_p2m_union =  melt(df_p2m_union[,c("Control","NH4OAc","FL","RNase")])
  x = df_p2m_union[which(df_p2m_union$variable==sample1),"value"]
  y = df_p2m_union[which(df_p2m_union$variable==sample2),"value"]
  
  if (log_log==T){
    x = log(x+1)
    y = log(y+1)
  }
  
  #print(cor(x,y))
  df = data.frame(x,y)
  linearMod <- lm(y ~ x, data = df)
  regression_label = paste0("y=",round(linearMod$coefficients[2],2),"x+",round(linearMod$coefficients[1],3))
  
  annotations <- data.frame(
    xpos = c(Inf,Inf),
    ypos =  c(-Inf,-Inf),
    annotateText = c(paste0("r=",round(cor(x,y),2)), regression_label),
    hjustvar = c(1.5, 1) ,
    vjustvar = c(-1, -3))
  
  p<-ggplot(df, aes(y=y, x=x)) +
    geom_point() +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=5) +
    labs(x = ifelse(log_log==F,paste0("P2M_",sample1),paste0("log(P2M+1)_",sample1)),
         y = ifelse(log_log==F,paste0("P2M_",sample2),paste0("log(P2M+1)_",sample2))) + 
    geom_smooth(method=lm) +
    # xlim(0,20) +
    # ylim(0,20) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.title = element_blank()) +
    ggtitle(paste0(sample2,"_",resolution))
  return(p)
}

p1 <- make_scatter_P2M(df_p2m_union_5000,
                       sample1="Control",
                       sample2="NH4OAc",
                       resolution="5000")
p2 <- make_scatter_P2M(df_p2m_union_5000,
                       sample1="Control",
                       sample2="FL",
                       resolution="5000")
p3 <- make_scatter_P2M(df_p2m_union_5000,
                       sample1="Control",
                       sample2="RNase",
                       resolution="5000")
p4 <- make_scatter_P2M(df_p2m_union_10000,
                       sample1="Control",
                       sample2="NH4OAc",
                       resolution="10000")
p5 <- make_scatter_P2M(df_p2m_union_10000,
                       sample1="Control",
                       sample2="FL",
                       resolution="10000")
p6 <- make_scatter_P2M(df_p2m_union_10000,
                       sample1="Control",
                       sample2="RNase",
                       resolution="10000")
png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2M_union_scatter_plots.png", width = 12, height = 8, units = "in", res = 200)
plot_grid(p2,p1,p3,p5,p4,p6, ncol = 3)
dev.off()


p1 <- make_scatter_P2M(df_p2m_union_5000,
                       sample1="Control",
                       sample2="NH4OAc",
                       resolution="5000",
                       log_log=T)
p2 <- make_scatter_P2M(df_p2m_union_5000,
                       sample1="Control",
                       sample2="FL",
                       resolution="5000",
                       log_log=T)
p3 <- make_scatter_P2M(df_p2m_union_5000,
                       sample1="Control",
                       sample2="RNase",
                       resolution="5000",
                       log_log=T)
p4 <- make_scatter_P2M(df_p2m_union_10000,
                       sample1="Control",
                       sample2="NH4OAc",
                       resolution="10000",
                       log_log=T)
p5 <- make_scatter_P2M(df_p2m_union_10000,
                       sample1="Control",
                       sample2="FL",
                       resolution="10000",
                       log_log=T)
p6 <- make_scatter_P2M(df_p2m_union_10000,
                       sample1="Control",
                       sample2="RNase",
                       resolution="10000",
                       log_log=T)
png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2M_union_scatter_plots_log_log.png", width = 12, height = 8, units = "in", res = 200)
plot_grid(p2,p1,p3,p5,p4,p6, ncol = 3)
dev.off()


### Scatter plot P2M with outliers that could be identified as differential loops
make_scatter_P2M_outliers <- function(df_p2m_union,
                                      sample1="Control",
                                      sample2,
                                      resolution,
                                      sd_factor=2,
                                      multiple_factor = 3){
  
  df_p2m_union_melt = melt(df_p2m_union[,c("Control","NH4OAc","FL","RNase")])
  
  x = df_p2m_union_melt[which(df_p2m_union_melt$variable==sample1),"value"]
  y = df_p2m_union_melt[which(df_p2m_union_melt$variable==sample2),"value"]
  
  #print(cor(x,y))
  df = data.frame(x,y)
  
  # df$diff = df$y - df$x
  # my_sd = sd(abs(df$diff))
  # df$position = "within"
  # df[which(abs(df$diff) > sd_factor*my_sd &
  #                         df$diff > 0),"position"] = "above"
  # df[which(abs(df$diff) > sd_factor*my_sd &
  #                         df$diff < 0),"position"] = "below"
  
  df$position = "within"
  df[which(df$y > multiple_factor*df$x),"position"] = "above"
  df[which(df$x > multiple_factor*df$y),"position"] = "below"
  
  above_num = as.character(table(as.factor(df$position))["above"])
  below_num = as.character(table(as.factor(df$position))["below"])
  above_perc = paste0(round(table(as.factor(df$position))["above"] / nrow(df) * 100,2),"%")
  below_perc = paste0(round(table(as.factor(df$position))["below"] / nrow(df) * 100,2),"%")
  
  annotations <- data.frame(
    xpos = c(Inf,-Inf),
    ypos =  c(-Inf,Inf),
    annotateText = c(paste0(below_num," (",below_perc,")"),
                     paste0(above_num," (",above_perc,")")),
    hjustvar = c(1.2,-0.2),
    vjustvar = c(-2,2))
  
  p <- ggplot(df) +
    geom_point(aes(y=y, x=x, color=position), size = 2, alpha = 0.5) +
    geom_line(aes(x=x, y=x)) +
    geom_line(aes(x=x, y=x*multiple_factor), linetype = "dashed") +
    geom_line(aes(x=x, y=x/multiple_factor), linetype = "dashed") +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=5) +
    labs(x = paste0("P2M_",sample1), y = paste0("P2M_",sample2)) +
    # xlim(0,20) +
    ylim(0,max(df$y)) +
    scale_color_manual(values=c('blue','red','gray')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=16),
          axis.text.x = element_text(size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none") +
    ggtitle(paste0(sample2,"_",resolution))
  return(p)
  
}

p1 <- make_scatter_P2M_outliers(df_p2m_union_5000,
                       sample1="Control",
                       sample2="NH4OAc",
                       resolution="5000")
p2 <- make_scatter_P2M_outliers(df_p2m_union_5000,
                       sample1="Control",
                       sample2="FL",
                       resolution="5000")
p3 <- make_scatter_P2M_outliers(df_p2m_union_5000,
                       sample1="Control",
                       sample2="RNase",
                       resolution="5000")
p4 <- make_scatter_P2M_outliers(df_p2m_union_10000,
                       sample1="Control",
                       sample2="NH4OAc",
                       resolution="10000")
p5 <- make_scatter_P2M_outliers(df_p2m_union_10000,
                       sample1="Control",
                       sample2="FL",
                       resolution="10000")
p6 <- make_scatter_P2M_outliers(df_p2m_union_10000,
                       sample1="Control",
                       sample2="RNase",
                       resolution="10000")
png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/P2M_union_scatter_plots_outliers.png", width = 12, height = 8, units = "in", res = 200)
plot_grid(p2,p1,p3,p5,p4,p6, ncol = 3)
dev.off()

### Add columns to label "differential peaks" based on P2M values
multiple_factor = 3

df_p2m_union_5000$NH4OAc_diff = ifelse(df_p2m_union_5000$NH4OAc > multiple_factor*df_p2m_union_5000$Control, "Up", ifelse(df_p2m_union_5000$NH4OAc < 1/multiple_factor*df_p2m_union_5000$Control,"Down","Not"))
df_p2m_union_5000$FL_diff = ifelse(df_p2m_union_5000$FL > multiple_factor*df_p2m_union_5000$Control, "Up", ifelse(df_p2m_union_5000$FL < 1/multiple_factor*df_p2m_union_5000$Control,"Down","Not"))
df_p2m_union_5000$RNase_diff = ifelse(df_p2m_union_5000$RNase > multiple_factor*df_p2m_union_5000$Control, "Up", ifelse(df_p2m_union_5000$RNase < 1/multiple_factor*df_p2m_union_5000$Control,"Down","Not"))

df_p2m_union_10000$NH4OAc_diff = ifelse(df_p2m_union_10000$NH4OAc > multiple_factor*df_p2m_union_10000$Control, "Up", ifelse(df_p2m_union_10000$NH4OAc < 1/multiple_factor*df_p2m_union_10000$Control,"Down","Not"))
df_p2m_union_10000$FL_diff = ifelse(df_p2m_union_10000$FL > multiple_factor*df_p2m_union_10000$Control, "Up", ifelse(df_p2m_union_10000$FL < 1/multiple_factor*df_p2m_union_10000$Control,"Down","Not"))
df_p2m_union_10000$RNase_diff = ifelse(df_p2m_union_10000$RNase > multiple_factor*df_p2m_union_10000$Control, "Up", ifelse(df_p2m_union_10000$RNase < 1/multiple_factor*df_p2m_union_10000$Control,"Down","Not"))


library(VennDiagram)
make_venn_diagram <- function(df_p2m_union,
                              resolution,
                              diff_type){
  
  NH4OAc <- which(df_p2m_union$NH4OAc_diff == diff_type)
  FL <- which(df_p2m_union$FL_diff == diff_type)
  RNase <- which(df_p2m_union$RNase_diff == diff_type)
  
  library(RColorBrewer)
  myCol <- brewer.pal(3, "Pastel2")
  
  # Chart
  venn.diagram(
    x = list(NH4OAc, FL, RNase),
    category.names = c(paste0("NH4OAc (",length(NH4OAc),")"),
                       paste0("FL (",length(FL),")"),
                       paste0("RNase (",length(RNase),")")),
    filename = paste0('/dataOS/rcalandrelli/phase_separation/loops/result/venn_diagrams/venn_diagram_',diff_type,"_",resolution,'.png'),
    output=T,
    
    # Main title
    main = diff_type,
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

make_venn_diagram(df_p2m_union_5000,"5000","Up")
make_venn_diagram(df_p2m_union_5000,"5000","Down")
make_venn_diagram(df_p2m_union_5000,"5000","Not")

make_venn_diagram(df_p2m_union_10000,"10000","Up")
make_venn_diagram(df_p2m_union_10000,"10000","Down")
make_venn_diagram(df_p2m_union_10000,"10000","Not")


diff_types = c("Up","Down","Not")
resolutions = c("5000","10000")
figures_list = list()
k = 1
for (i in resolutions){
  for (j in diff_types){
    figures_list[[k]] = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/loops/result/venn_diagrams/venn_diagram_",j,"_",i,".png"))
    k = k + 1
  }
}

png("/dataOS/rcalandrelli/phase_separation/loops/result/venn_diagrams_5000.png", width = 8, height = 2.2, units = "in", res = 200)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]),
             rasterGrob(figures_list[[3]]), ncol=3, left = "Res: 5000")
dev.off()

png("/dataOS/rcalandrelli/phase_separation/loops/result/venn_diagrams_10000.png", width = 8, height = 2.2, units = "in", res = 200)
grid.arrange(rasterGrob(figures_list[[4]]),
             rasterGrob(figures_list[[5]]),
             rasterGrob(figures_list[[6]]), ncol=3, left = "Res: 10000")
dev.off()

### 3D scatter plots relative P2M versus control
df_p2m_union_5000_rel = df_p2m_union_5000
df_p2m_union_5000_rel$NH4OAc = df_p2m_union_5000_rel$NH4OAc / df_p2m_union_5000_rel$Control
df_p2m_union_5000_rel$FL = df_p2m_union_5000_rel$FL / df_p2m_union_5000_rel$Control
df_p2m_union_5000_rel$RNase = df_p2m_union_5000_rel$RNase / df_p2m_union_5000_rel$Control

df_p2m_union_10000_rel = df_p2m_union_10000
df_p2m_union_10000_rel$NH4OAc = df_p2m_union_10000_rel$NH4OAc / df_p2m_union_10000_rel$Control
df_p2m_union_10000_rel$FL = df_p2m_union_10000_rel$FL / df_p2m_union_10000_rel$Control
df_p2m_union_10000_rel$RNase = df_p2m_union_10000_rel$RNase / df_p2m_union_10000_rel$Control

library(scatterplot3d)
png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/scatter3D_P2M_5000.png", width = 5, height = 5, units = "in", res = 200)
scatterplot3d(df_p2m_union_5000_rel[,2:4], main = "Resolution: 5000", cex.axis = 1, cex.lab = 1,
              xlab = "P2M_NH4OAc / P2M_control", ylab = "P2M_FL / P2M_control", zlab = "P2M_RNase / P2M_control")
dev.off()

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/scatter3D_P2M_10000.png", width = 5, height = 5, units = "in", res = 200)
scatterplot3d(df_p2m_union_10000_rel[,2:4], main = "Resolution: 10000", cex.axis = 1, cex.lab = 1,
              xlab = "P2M_NH4OAc / P2M_control", ylab = "P2M_FL / P2M_control", zlab = "P2M_RNase / P2M_control")
dev.off()

figures_list = list()
k = 1
for (i in resolutions){
  figures_list[[k]] = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/scatter3D_P2M_",i,".png"))
  k = k + 1
}

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/scatter3D_plots_grid.png", width = 10, height = 5, units = "in", res = 200)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]), ncol=2)
dev.off()


###### Same information but on ternary plots
temp = df_p2m_union_5000_rel[,2:4]
temp = temp / rowSums(temp)

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/ternary_P2M_5000.png", width = 5, height = 5, units = "in", res = 200)
ggtern(data=df_p2m_union_5000_rel, aes(x=FL,y=NH4OAc, z=RNase)) +
  geom_point(size=0.5) +
  theme_rgbw() +
  labs(title="Resolution: 5000") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/ternary_P2M_10000.png", width = 5, height = 5, units = "in", res = 200)
ggtern(data=df_p2m_union_10000_rel, aes(x=FL,y=NH4OAc, z=RNase)) +
  geom_point(size=0.5) +
  theme_rgbw() +
  labs(title="Resolution: 10000") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

figures_list = list()
k = 1
for (i in c(5000,10000)){
  figures_list[[k]] = readPNG(paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/ternary_P2M_",i,".png"))
  k = k + 1
}

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/result/ternary_plots_grid.png", width = 10, height = 5, units = "in", res = 200)
grid.arrange(rasterGrob(figures_list[[1]]),
             rasterGrob(figures_list[[2]]), ncol=2)
dev.off()


######################################### Enhancer-promoter loops

##### Enhancers
H1_enhancers = read.table("/dataOS/rcalandrelli/phase_separation/H1_enhancers_hg38.bed")
H1_enhancers$V2 = H1_enhancers$V2
H1_enhancers = H1_enhancers[which(H1_enhancers$V1 %in% hg38_chromosomes),] # 58596 enhancers 
Gr_H1_enhancers = GRanges(
  seqnames = Rle(H1_enhancers$V1),
  ranges = IRanges(H1_enhancers$V2 + 1, end = H1_enhancers$V3, names = c(1:nrow(H1_enhancers))),
  strand = Rle(strand('*')))

H1_enhancers[which(H1_enhancers$V1 == "chr20" & 
                     H1_enhancers$V2 > 37260000 & 
                     H1_enhancers$V3 < 37290000),]

##### Promoters
annotations_edb_hg38 <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                              keys = keys(EnsDb.Hsapiens.v86),
                                              columns = c("SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND", "SYMBOL", "ENTREZID","GENEBIOTYPE"),
                                              keytype = "GENEID")

annotations_edb_hg38$SEQNAME = paste0("chr", annotations_edb_hg38$SEQNAME)
annotations_edb_hg38 = annotations_edb_hg38[which(annotations_edb_hg38$SEQNAME %in% hg38_chromosomes),]

TxDb.Hsapiens.gencode.24 = makeTxDbFromGFF('/home/frankyan/research/refGenome/standard4DN/GRCh38/gencode.v24.primary_assembly.annotation.gtf')
gencode.24_genes = genes(TxDb.Hsapiens.gencode.24)
gencode.24_genes = gencode.24_genes[seqnames(gencode.24_genes) %in% hg38_chromosomes]
gencode.24_promoters = promoters(gencode.24_genes, upstream = 2000, downstream = 2000, use.names = T)


##### Loop data
Gr_loop_control_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_control$V1)),
  ranges = IRanges(loops_control$V2+1, end = loops_control$V3, names = c(1:nrow(loops_control))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control)))
Gr_loop_control_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_control$V1)),
  ranges = IRanges(loops_control$V5+1, end = loops_control$V6, names = c(1:nrow(loops_control))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_control)))

Gr_loop_NH4OAc_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_NH4OAc$V1)),
  ranges = IRanges(loops_NH4OAc$V2+1, end = loops_NH4OAc$V3, names = c(1:nrow(loops_NH4OAc))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_NH4OAc)))
Gr_loop_NH4OAc_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_NH4OAc$V1)),
  ranges = IRanges(loops_NH4OAc$V5+1, end = loops_NH4OAc$V6, names = c(1:nrow(loops_NH4OAc))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_NH4OAc)))

Gr_loop_FL_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_FL$V1)),
  ranges = IRanges(loops_FL$V2+1, end = loops_FL$V3, names = c(1:nrow(loops_FL))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_FL)))
Gr_loop_FL_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_FL$V1)),
  ranges = IRanges(loops_FL$V5+1, end = loops_FL$V6, names = c(1:nrow(loops_FL))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_FL)))

Gr_loop_RNase_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_RNase$V1)),
  ranges = IRanges(loops_RNase$V2+1, end = loops_RNase$V3, names = c(1:nrow(loops_RNase))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_RNase)))
Gr_loop_RNase_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_RNase$V1)),
  ranges = IRanges(loops_RNase$V5+1, end = loops_RNase$V6, names = c(1:nrow(loops_RNase))),
  strand = Rle(strand('*')),
  index=c(1:nrow(loops_RNase)))


##### Detect promoter-enhancer (or vice versa) loops

find_promoter_enhancer_loop <- function(Gr_loop_1, # loop anchor 1
                                        Gr_loop_2, # loop anchor 2
                                        loop_file,
                                        sample_name,
                                        Gr_enhancer=Gr_H1_enhancers,
                                        Gr_promoter=gencode.24_promoters){
  
  loop_file$index = 1:nrow(loop_file)
  
  ### Promoter-enhancer
  overlaps = countOverlaps(Gr_loop_1, Gr_promoter)
  Gr_loop_1_promoter = Gr_loop_1[overlaps>0]
  
  overlaps = countOverlaps(Gr_loop_2, Gr_enhancer)
  Gr_loop_2_enhancer = Gr_loop_2[overlaps>0]
  
  loop_promoter_enhancer = loop_file[intersect(names(Gr_loop_1_promoter), names(Gr_loop_2_enhancer)),]
  loop_promoter_enhancer$V25 = "promoter_enhancer"
  
  ### Enhancer-promoter
  overlaps = countOverlaps(Gr_loop_1, Gr_enhancer)
  Gr_loop_1_enhancer = Gr_loop_1[overlaps>0]
  
  overlaps = countOverlaps(Gr_loop_2, Gr_promoter)
  Gr_loop_2_promoter = Gr_loop_2[overlaps>0]
  
  loop_enhancer_promoter = loop_file[intersect(names(Gr_loop_1_enhancer), names(Gr_loop_2_promoter)),]
  loop_enhancer_promoter$V25 = "enhancer_promoter"
  
  ### Concatenate data
  temp1 = intersect(names(Gr_loop_1_promoter), names(Gr_loop_2_enhancer))
  temp2 = intersect(names(Gr_loop_1_enhancer), names(Gr_loop_2_promoter))
  temp = intersect(temp1,temp2)
  
  if (length(temp) > 0){
    loop_out = rbind(loop_promoter_enhancer[setdiff(rownames(loop_promoter_enhancer),temp),],
                     loop_enhancer_promoter)
    loop_out[temp,"V26"] = "both"
  } else {
    loop_out = rbind(loop_promoter_enhancer, loop_enhancer_promoter)
  }
  
  loops_control_full = loops_control
  loop_file$V26 = "NA"
  loop_file[rownames(loop_out), "V26"] = loop_out$V25

  
  write.table(loop_file, paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/",sample_name,"_EP.bedpe"), row.names = F, col.names = F, sep = "\t", quote = F)
  return(loop_file)
}

loop_control_EP = find_promoter_enhancer_loop(Gr_loop_control_1, 
                                              Gr_loop_control_2,
                                              loops_control,
                                              "H1_control",
                                              Gr_enhancer=Gr_H1_enhancers,
                                              Gr_promoter=gencode.24_promoters)

loop_NH4OAc_EP = find_promoter_enhancer_loop(Gr_loop_NH4OAc_1, 
                                              Gr_loop_NH4OAc_2,
                                              loops_NH4OAc,
                                              "H1_NH4OAc",
                                              Gr_enhancer=Gr_H1_enhancers,
                                              Gr_promoter=gencode.24_promoters)

loop_FL_EP = find_promoter_enhancer_loop(Gr_loop_FL_1, 
                                              Gr_loop_FL_2,
                                              loops_FL,
                                              "H1_FL",
                                              Gr_enhancer=Gr_H1_enhancers,
                                              Gr_promoter=gencode.24_promoters)

loop_RNase_EP = find_promoter_enhancer_loop(Gr_loop_RNase_1, 
                                              Gr_loop_RNase_2,
                                              loops_RNase,
                                              "H1_RNase",
                                              Gr_enhancer=Gr_H1_enhancers,
                                              Gr_promoter=gencode.24_promoters)


loop_control_EP[which(loop_control_EP$V1 == 20 & loop_control_EP$V2 > 37000000 & loop_control_EP$V3 < 38000000),]

# loops_control_full = loops_control
# loops_control_full$V25 = "NA"
# loops_control_full[rownames(loop_control_EP), "V25"] = loop_control_EP$V25
# write.table(loops_control_full, "/dataOS/rcalandrelli/phase_separation/HiC/loops/loops_with_EP_info/loops_control_EP_label.bedpe", row.names = F, col.names = F, quote = F, sep ="\t")

loops_control_higlass = loops_control[,1:6]
loops_control_higlass$V1 = paste0("chr",loops_control_higlass$V1)
loops_control_higlass$V4 = paste0("chr",loops_control_higlass$V4)
loops_control_higlass$V3 = loops_control_higlass$V6
loops_control_higlass$V5 = loops_control_higlass$V2
write.table(loops_control_higlass, "/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_control_merged/control_merged_loops_higlass.bedpe", row.names = F, col.names = F, sep = "\t", quote = F)



############### Detect pattern 1 (high density regions within loops) using a "in vs out" approach starting from loop coordinates

calculate_HDR_loop_hic <- function(loop_file,
                                   chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                   chromosomes,
                                   sample,
                                   bin_size,
                                   out_template_factor,
                                   shuffle_loops = F
)
{
  
  out_dir_HDR = paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_", sample, "_merged")
  
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
    if (file.info(paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/", sample, "_merged/", bin_size, "/", my_chr, "_", my_chr, "_", bin_size, ".txt"))$size != 0){
      temp = data.frame(fread(paste0("/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/", sample, "_merged/", bin_size, "/", my_chr, "_", my_chr, "_", bin_size, ".txt")))
      contact_matrix = matrix(0, nrow = length(chr_genome_window), ncol = length(chr_genome_window))
      contact_matrix[as.matrix(temp[,c("V1","V2")] / bin_size + 1)] = temp$V3
      contact_matrix[as.matrix(temp[,c("V2","V1")] / bin_size + 1)] = temp$V3
      contact_matrix[which(is.na(contact_matrix), arr.ind = T)] = 0
      
      ########## Load loop data
      loops = read.table(loop_file)
      loops = loops[which(loops$V1 == gsub("chr","",my_chr)),]
      
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
          } 
          else {
            ratio_loops[[i]] = NA
          }
        }
        
        loops$density_ratio_hic = do.call("rbind", ratio_loops)
        loop_chr_list[[my_chr]] = loops
      }
    }
  }
  
  loops_out = do.call("rbind", loop_chr_list)
  write.table(loops_out, paste0(out_dir_HDR, "/loops_EP_ratio_template_factor_", out_template_factor,".bedpe"), row.names = F, col.names = F, sep= "\t", quote = F)
  return(loops_out)
}


loop_control_density_ratio_hic = calculate_HDR_loop_hic(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_control_EP.bedpe",
                                                  chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                  chromosomes=hg38_chromosomes[1:23],
                                                  sample="H1_control",
                                                  bin_size=10000,
                                                  out_template_factor=2)

loop_NH4OAc_density_ratio_hic = calculate_HDR_loop_hic(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_NH4OAc_EP.bedpe",
                                                        chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                        chromosomes=hg38_chromosomes[1:23],
                                                        sample="H1_NH4OAc",
                                                        bin_size=10000,
                                                        out_template_factor=2)

loop_FL_density_ratio_hic = calculate_HDR_loop_hic(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_FL_EP.bedpe",
                                                        chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                        chromosomes=hg38_chromosomes[1:23],
                                                        sample="H1_FL",
                                                        bin_size=10000,
                                                        out_template_factor=2)

loop_RNase_density_ratio_hic = calculate_HDR_loop_hic(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_RNase_EP.bedpe",
                                                        chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                        chromosomes=hg38_chromosomes[1:23],
                                                        sample="H1_RNase",
                                                        bin_size=10000,
                                                        out_template_factor=2)


loop_control_density_ratio_hic_shuffle = calculate_HDR_loop_hic(loop_file="/dataOS/rcalandrelli/phase_separation/HiC/loops/enhancer_promoter_loops/H1_control_EP.bedpe",
                                                        chrSize_file="/dataOS/rcalandrelli/phase_separation/4DNFI823LSII.hg38.mainonly.chrom.sizes",
                                                        chromosomes=hg38_chromosomes[1:23],
                                                        sample="control",
                                                        bin_size=10000,
                                                        out_template_factor=2,
                                                        shuffle_loops = T)



loop_control_density_ratio_hic$density_ratio_hic_shuffle = loop_control_density_ratio_hic_shuffle$density_ratio_hic


loop_control_density_ratio_hic[which(loop_control_density_ratio_hic$V1 == 1 & loop_control_density_ratio_hic$V2 > 3250000 & loops_control_EP_ratio$V6 < 3750000),]
head(loop_control_density_ratio_hic[which(loop_control_density_ratio_hic$density_ratio_hic > 3.683),],20)


### Distribution of density ratio
quantile(loop_control_density_ratio_hic$density_ratio_hic, 0.75)
nrow(loop_control_density_ratio_hic[which(loop_control_density_ratio_hic$density_ratio_hic > 3.7),])
nrow(loop_control_density_ratio_hic_shuffle[which(loop_control_density_ratio_hic_shuffle$density_ratio_hic > 3.7),])

wilcox.test(loop_control_density_ratio_hic$density_ratio_hic, loop_control_density_ratio_hic_shuffle$density_ratio_hic, alternative = "greater")

p1<-ggplot(loop_control_density_ratio_hic, aes(x=density_ratio_hic)) + 
  geom_histogram(color="black", fill="blue", bins = 100) + 
  labs(x="HiC density ratio", y=expression("Number of loops")) +
  #scale_x_continuous(name="Density ratio", breaks=c(0,5,10,15,20,25), labels=c("0","5","10","15","20","25")) +
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
temp_melt = reshape2::melt(loop_control_density_ratio_hic[,c("density_ratio_hic", "density_ratio_hic_shuffle")])
temp_melt$value = as.numeric(temp_melt$value)
temp_melt$variable = gsub("density_ratio_hic_shuffle","Shuffle data",temp_melt$variable)
temp_melt$variable = gsub("density_ratio_hic","Real data",temp_melt$variable)

summary(loop_control_density_ratio_hic$density_ratio_hic)

p2<-ggplot(temp_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) +
  stat_compare_means(method = "wilcox", size = 6, label.y = 25) +
  ylim(c(0,26)) +
  labs(x="", y=expression("HiC density ratio")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=16),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_control_merged/HiC_density_ratio_hist_box.png", width = 6, height = 8, res = 200, units = "in")
plot_grid(p1,p2, nrow=2, labels = c("B","C"), label_size = 18)
dev.off()



############### Detect loop domains using Rao et al. method

find_loop_domain <- function(Gr_loop_1, # loop anchor 1
                             Gr_loop_2, # loop anchor 2
                             loop_file,
                             sample_name,
                             Gr_tad){
  tad_file = data.frame(Gr_tad)
  
  ### TAD start coordinate (width + 10000 because the tad spans until the end of the last bin)
  tad_1 = tad_file
  tad_1$start = ifelse(tad_file$width + 10000 > 250000, tad_file$start + 1 - 25000, tad_file$start + 1 - 0.1*(tad_file$width-1))
  tad_1$end = ifelse(tad_file$width + 10000 > 250000, tad_file$start + 1 + 25000, tad_file$start + 1 + 0.1*(tad_file$width-1))
  Gr_tad_1 = GRanges(
    seqnames = Rle(as.character(tad_1$seqnames)),
    ranges = IRanges(as.numeric(tad_1$start), end = as.numeric(tad_1$end), names = c(1:nrow(tad_1))),
    strand = Rle(strand('*')))
  
  ### TAD end coordinate
  tad_2 = tad_file
  tad_2$start = ifelse(tad_file$width + 10000 > 250000, tad_file$end - 25000, tad_file$end - 0.1*(tad_file$width-1))
  tad_2$end = ifelse(tad_file$width + 10000 > 250000, tad_file$end + 25000, tad_file$end + 0.1*(tad_file$width-1))
  Gr_tad_2 = GRanges(
    seqnames = Rle(as.character(tad_2$seqnames)),
    ranges = IRanges(as.numeric(tad_2$start), end = as.numeric(tad_2$end), names = c(1:nrow(tad_2))),
    strand = Rle(strand('*')))
  
  ### Loop domain
  overlaps = countOverlaps(Gr_loop_1, Gr_tad_1)
  Gr_loop_1_tad = Gr_loop_1[overlaps>0]
  
  overlaps = countOverlaps(Gr_loop_2, Gr_tad_2)
  Gr_loop_2_tad = Gr_loop_2[overlaps>0]
  
  loop_file$V25 = "NA"
  loop_file[intersect(names(Gr_loop_1_tad), names(Gr_loop_2_tad)),25] = "loop_domain"
  write.table(loop_file, paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/loop_domain/",sample_name,"_loop_domain.bedpe"), row.names = F, col.names = F, sep = "\t", quote = F)
  
  ### Contact domains associated with loop domains
  overlaps = data.frame(findOverlaps(Gr_tad_1, Gr_loop_1))
  tads_with_loop = overlaps[which(overlaps$subjectHits %in% temp),"queryHits"]
  tad_file_loop = tad_file[tads_with_loop,]
  write.table(tad_file_loop, paste0("/dataOS/rcalandrelli/phase_separation/HiC/loops/loop_domain/",sample_name,"_domains_with_corner_loop.bed"), row.names = F, col.names = T, sep = "\t", quote = F)
  
  return(loop_file)
}

tad_file[which(tad_file$seqnames == "chr1" &
                 tad_file$start > 44600000 &
                 tad_file$end < 45600000), ]


loop_domain_control = find_loop_domain(Gr_loop_control_1, 
                                       Gr_loop_control_2,
                                       loops_control,
                                       "H1_control",
                                       Gr_tads_control_10000)

loop_domain_NH4OAc = find_loop_domain(Gr_loop_NH4OAc_1, 
                                       Gr_loop_NH4OAc_2,
                                       loops_NH4OAc,
                                       "H1_NH4OAc",
                                       Gr_tads_NH4OAc_10000)

loop_domain_FL = find_loop_domain(Gr_loop_FL_1, 
                                       Gr_loop_FL_2,
                                       loops_FL,
                                       "H1_FL",
                                       Gr_tads_FL_10000)

loop_domain_RNase = find_loop_domain(Gr_loop_RNase_1, 
                                       Gr_loop_RNase_2,
                                       loops_RNase,
                                       "H1_RNase",
                                       Gr_tads_RNase_10000)

nrow(loop_domain_control[which(loop_domain_control$V25 == "loop_domain"),])

loop_domain_control[which(loop_domain_control$V25 == "loop_domain" & loop_domain_control$V1 == 1),]


df = data.frame(group = c("Loops", "Loop domains"),
                value = c(1814-620, 620))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/loop_domain/pie_chart_loop_domains.png", width = 5, height = 5, units = "in", res = 200)
ggplot(df, aes(x="", y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  blank_theme +
  theme(axis.text.x=element_blank(),
        legend.title = element_blank(),
        text = element_text(size=16),
        legend.position = "top") +
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
      label = value), size=6)
dev.off()


######## TADs with corner loop insulation score

tad_with_corner_loop = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/loop_domain/H1_control_domains_with_corner_loop.bed", header = T)

boundaries_with_scores_control$corner_loop = "No"
for (i in 1:nrow(boundaries_with_scores_control)){
  i_chr = boundaries_with_scores_control[i,"chr"]
  i_coord = boundaries_with_scores_control[i,"coord"]
  if (nrow(tad_with_corner_loop[which((tad_with_corner_loop$seqnames == i_chr & tad_with_corner_loop$start == i_coord) |
                                 (tad_with_corner_loop$seqnames == i_chr & tad_with_corner_loop$end == i_coord)),]) > 0){
    boundaries_with_scores_control[i,"corner_loop"] = "Yes"
  } 
}

# Insulation score domains with corner loop
x = boundaries_with_scores_control[which(boundaries_with_scores_control$corner_loop=="Yes"),"insulation_score"]
# Insulation score domains without corner loop
y = boundaries_with_scores_control[which(boundaries_with_scores_control$corner_loop=="No"),"insulation_score"]
y = y[!is.infinite(y)]

wilcox.test(x, y)


df1 = data.frame(x,"Corner loop")
colnames(df1) = c("insulation_score","type")
df2 = data.frame(y, "No corner loop")
colnames(df2) = c("insulation_score","type")

df = rbind(df1,df2)

# Saved into /dataOS/rcalandrelli/phase_separation/MARGI/loops/iMARGI_H1_control/density_ratio_violin_HiC_MARGI.png in script margi_loops.r
p3 <- ggplot(df, aes(x=type, y=insulation_score, fill=type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) +
  stat_compare_means(method = "wilcox", size = 6, label.y = 175, label = "p.signif") +
  ylim(c(0,180)) +
  labs(x="", y="TAD boundary insulation score") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=labelsize),
        axis.text.x = element_text(size = labelsize, color = "black"),
        axis.text.y = element_text(size = labelsize, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")




#########################
### Analysis of differential loops in FL and RNase

nrow(loops_control[which(loops_control$V3- loops_control$V2 == 10000),])
nrow(loops_NH4OAc[which(loops_NH4OAc$V3- loops_NH4OAc$V2 == 10000),])
nrow(loops_FL[which(loops_FL$V3- loops_FL$V2 == 10000),])
nrow(loops_RNase[which(loops_RNase$V3- loops_RNase$V2 == 10000),])

temp = rbind(loops_control,loops_NH4OAc,loops_FL,loops_RNase)
loops_union = temp[!duplicated(temp[,c(1,2,3,4,5,6)]),]
loops_union$V1 = as.character(loops_union$V1)
loops_union$V4 = as.character(loops_union$V4)

loops_union_5000 = loops_union[which(loops_union$V3 - loops_union$V2 == 5000),]
loops_union_10000 = loops_union[which(loops_union$V3 - loops_union$V2 == 10000),]

Gr_loop_union_10000_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_10000$V1)),
  ranges = IRanges(loops_union_10000$V2+1-10000, end = loops_union_10000$V3+10000, names = c(1:nrow(loops_union_10000))),
  strand = Rle(strand('*')))

Gr_loop_union_10000_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_10000$V1)),
  ranges = IRanges(loops_union_10000$V5+1-10000, end = loops_union_10000$V6+10000, names = c(1:nrow(loops_union_10000))),
  strand = Rle(strand('*')))

Gr_loop_union_5000_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_5000$V1)),
  ranges = IRanges(loops_union_5000$V2+1-5000, end = loops_union_5000$V3+5000, names = c(1:nrow(loops_union_5000))),
  strand = Rle(strand('*')))

Gr_loop_union_5000_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_5000$V1)),
  ranges = IRanges(loops_union_5000$V5+1-5000, end = loops_union_5000$V6+5000, names = c(1:nrow(loops_union_5000))),
  strand = Rle(strand('*')))


##### 10000
make_upset_list_10000 <- function(i){
  # We need to make sure that both anchors come from the same loop
  temp_control = ifelse(length(intersect(which(countOverlaps(Gr_loop_control_1, Gr_loop_union_10000_1[i])>0), which(countOverlaps(Gr_loop_control_2, Gr_loop_union_10000_2[i])>0))) >= 1,
                        1,
                        0)
  
  temp_FL = ifelse(length(intersect(which(countOverlaps(Gr_loop_FL_1, Gr_loop_union_10000_1[i])>0), which(countOverlaps(Gr_loop_FL_2, Gr_loop_union_10000_2[i])>0))) >= 1,
                   1,
                   0)
  
  temp_NH4OAc = ifelse(length(intersect(which(countOverlaps(Gr_loop_NH4OAc_1, Gr_loop_union_10000_1[i])>0), which(countOverlaps(Gr_loop_NH4OAc_2, Gr_loop_union_10000_2[i])>0))) >= 1,
                       1,
                       0)
  
  temp_RNase = ifelse(length(intersect(which(countOverlaps(Gr_loop_RNase_1, Gr_loop_union_10000_1[i])>0), which(countOverlaps(Gr_loop_RNase_2, Gr_loop_union_10000_2[i])>0))) >= 1,
                      1,
                      0)
  
  out_vect = as.numeric(c(temp_control,
                          temp_FL,
                          temp_NH4OAc,
                          temp_RNase))
  return(out_vect)
}

df_upset_H1_loops_10000_list = mclapply(c(1:length(Gr_loop_union_10000_1)), make_upset_list_10000, mc.cores = 32)
df_upset_H1_loops_10000 = data.frame(do.call("rbind", df_upset_H1_loops_10000_list))
colnames(df_upset_H1_loops_10000) = c("Control","FL","NH4OAc","RNase")
df_upset_H1_loops_10000_full = cbind(data.frame(Gr_loop_union_10000_1)[,c(1,2,3,4)], data.frame(Gr_loop_union_10000_2)[,c(1,2,3,4)], df_upset_H1_loops_10000)
colnames(df_upset_H1_loops_10000_full)[1:8] = c("chr1","start1","end1","width1","chr2","start2","end2","width2")

library(UpSetR)
png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/loop_H1_10000_extended_upset.png", width = 5, height = 5.2, units = "in", res = 400)
upset(df_upset_H1_loops_10000_full, 
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


##### 5000

make_upset_list_5000 <- function(i){
  # temp_control = ifelse(sum(countOverlaps(Gr_loop_control_1, Gr_loop_union_5000_1[i])) > 0 & sum(countOverlaps(Gr_loop_control_2, Gr_loop_union_5000_2[i])) > 0,
  #                       1,
  #                       0)
  # 
  # temp_FL = ifelse(sum(countOverlaps(Gr_loop_FL_1, Gr_loop_union_5000_1[i])) > 0 & sum(countOverlaps(Gr_loop_FL_2, Gr_loop_union_5000_2[i])) > 0,
  #                  1,
  #                  0)
  # 
  # temp_NH4OAc = ifelse(sum(countOverlaps(Gr_loop_NH4OAc_1, Gr_loop_union_5000_1[i])) > 0 & sum(countOverlaps(Gr_loop_NH4OAc_2, Gr_loop_union_5000_2[i])) > 0,
  #                      1,
  #                      0)
  # 
  # temp_RNase = ifelse(sum(countOverlaps(Gr_loop_RNase_1, Gr_loop_union_5000_1[i])) > 0 & sum(countOverlaps(Gr_loop_RNase_2, Gr_loop_union_5000_2[i])) > 0,
  #                     1,
  #                     0)
  
  # We need to make sure that both anchors come from the same loop
  temp_control = ifelse(length(intersect(which(countOverlaps(Gr_loop_control_1, Gr_loop_union_5000_1[i])>0), which(countOverlaps(Gr_loop_control_2, Gr_loop_union_5000_2[i])>0))) >= 1,
                        1,
                        0)
  
  temp_FL = ifelse(length(intersect(which(countOverlaps(Gr_loop_FL_1, Gr_loop_union_5000_1[i])>0), which(countOverlaps(Gr_loop_FL_2, Gr_loop_union_5000_2[i])>0))) >= 1,
                   1,
                   0)
  
  temp_NH4OAc = ifelse(length(intersect(which(countOverlaps(Gr_loop_NH4OAc_1, Gr_loop_union_5000_1[i])>0), which(countOverlaps(Gr_loop_NH4OAc_2, Gr_loop_union_5000_2[i])>0))) >= 1,
                       1,
                       0)
  
  temp_RNase = ifelse(length(intersect(which(countOverlaps(Gr_loop_RNase_1, Gr_loop_union_5000_1[i])>0), which(countOverlaps(Gr_loop_RNase_2, Gr_loop_union_5000_2[i])>0))) >= 1,
                      1,
                      0)
  
  out_vect = as.numeric(c(temp_control,
                          temp_FL,
                          temp_NH4OAc,
                          temp_RNase))
  return(out_vect)
}

df_upset_H1_loops_5000_list = mclapply(c(1:length(Gr_loop_union_5000_1)), make_upset_list_5000, mc.cores = 32)
df_upset_H1_loops_5000 = data.frame(do.call("rbind", df_upset_H1_loops_5000_list))
colnames(df_upset_H1_loops_5000) = c("Control","FL","NH4OAc","RNase")
df_upset_H1_loops_5000_full = cbind(data.frame(Gr_loop_union_5000_1)[,c(1,2,3,4)], data.frame(Gr_loop_union_5000_2)[,c(1,2,3,4)], df_upset_H1_loops_5000)
colnames(df_upset_H1_loops_5000_full)[1:8] = c("chr1","start1","end1","width1","chr2","start2","end2","width2")


png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/loop_H1_5000_extended_upset.png", width = 5, height = 5.2, units = "in", res = 400)
upset(df_upset_H1_loops_5000, 
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


################# Compare differential loops in FL and RNase with other cell lines
loops_HFF = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_HFF/merged_loops.bedpe")
Gr_loop_HFF_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_HFF$V1)),
  ranges = IRanges(loops_HFF$V2+1, end = loops_HFF$V3, names = c(1:nrow(loops_HFF))),
  strand = Rle(strand('*')))
Gr_loop_HFF_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_HFF$V1)),
  ranges = IRanges(loops_HFF$V5+1, end = loops_HFF$V6, names = c(1:nrow(loops_HFF))),
  strand = Rle(strand('*')))

loops_K562 = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_K562/merged_loops.bedpe")
Gr_loop_K562_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_K562$V1)),
  ranges = IRanges(loops_K562$V2+1, end = loops_K562$V3, names = c(1:nrow(loops_K562))),
  strand = Rle(strand('*')))
Gr_loop_K562_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_K562$V1)),
  ranges = IRanges(loops_K562$V5+1, end = loops_K562$V6, names = c(1:nrow(loops_K562))),
  strand = Rle(strand('*')))

loops_H1_ES_MicroC = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_ES_MicroC/merged_loops.bedpe")
Gr_loop_H1_ES_MicroC_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_H1_ES_MicroC$V1)),
  ranges = IRanges(loops_H1_ES_MicroC$V2+1, end = loops_H1_ES_MicroC$V3, names = c(1:nrow(loops_H1_ES_MicroC))),
  strand = Rle(strand('*')))
Gr_loop_H1_ES_MicroC_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_H1_ES_MicroC$V1)),
  ranges = IRanges(loops_H1_ES_MicroC$V5+1, end = loops_H1_ES_MicroC$V6, names = c(1:nrow(loops_H1_ES_MicroC))),
  strand = Rle(strand('*')))

loops_H1_endoderm_MicroC = read.table("/dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_endoderm_MicroC/merged_loops.bedpe")
Gr_loop_H1_endoderm_MicroC_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_H1_endoderm_MicroC$V1)),
  ranges = IRanges(loops_H1_endoderm_MicroC$V2+1, end = loops_H1_endoderm_MicroC$V3, names = c(1:nrow(loops_H1_endoderm_MicroC))),
  strand = Rle(strand('*')))
Gr_loop_H1_endoderm_MicroC_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_H1_endoderm_MicroC$V1)),
  ranges = IRanges(loops_H1_endoderm_MicroC$V5+1, end = loops_H1_endoderm_MicroC$V6, names = c(1:nrow(loops_H1_endoderm_MicroC))),
  strand = Rle(strand('*')))

print_loop_numbers(loops_HFF)
print_loop_numbers(loops_K562)
print_loop_numbers(loops_H1_ES_MicroC)
print_loop_numbers(loops_H1_endoderm_MicroC)


temp = rbind(loops_control,loops_NH4OAc,loops_FL,loops_RNase,loops_HFF,loops_K562,loops_H1_ES_MicroC,loops_H1_endoderm_MicroC)
loops_union_ALL = temp[!duplicated(temp[,c(1,2,3,4,5,6)]),]
loops_union_ALL$V1 = as.character(loops_union_ALL$V1)
loops_union_ALL$V4 = as.character(loops_union_ALL$V4)

loops_union_ALL_5000 = loops_union_ALL[which(loops_union_ALL$V3 - loops_union_ALL$V2 == 5000),]
loops_union_ALL_10000 = loops_union_ALL[which(loops_union_ALL$V3 - loops_union_ALL$V2 == 10000),]

Gr_loop_union_10000_ALL_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_ALL_10000$V1)),
  ranges = IRanges(loops_union_ALL_10000$V2+1-10000, end = loops_union_ALL_10000$V3+10000, names = c(1:nrow(loops_union_ALL_10000))),
  strand = Rle(strand('*')))

Gr_loop_union_10000_ALL_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_ALL_10000$V1)),
  ranges = IRanges(loops_union_ALL_10000$V5+1-10000, end = loops_union_ALL_10000$V6+10000, names = c(1:nrow(loops_union_ALL_10000))),
  strand = Rle(strand('*')))

Gr_loop_union_5000_ALL_1 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_ALL_5000$V1)),
  ranges = IRanges(loops_union_ALL_5000$V2+1-5000, end = loops_union_ALL_5000$V3+5000, names = c(1:nrow(loops_union_ALL_5000))),
  strand = Rle(strand('*')))

Gr_loop_union_5000_ALL_2 = GRanges(
  seqnames = Rle(paste0("chr",loops_union_ALL_5000$V1)),
  ranges = IRanges(loops_union_ALL_5000$V5+1-5000, end = loops_union_ALL_5000$V6+5000, names = c(1:nrow(loops_union_ALL_5000))),
  strand = Rle(strand('*')))


make_upset_list_ALL <- function(i){
  
  temp_control = ifelse(length(intersect(which(countOverlaps(Gr_loop_control_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_control_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                        1,
                        0)
  
  temp_FL = ifelse(length(intersect(which(countOverlaps(Gr_loop_FL_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_FL_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                   1,
                   0)
  
  temp_NH4OAc = ifelse(length(intersect(which(countOverlaps(Gr_loop_NH4OAc_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_NH4OAc_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                       1,
                       0)
  
  temp_RNase = ifelse(length(intersect(which(countOverlaps(Gr_loop_RNase_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_RNase_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                      1,
                      0)
  
  temp_HFF = ifelse(length(intersect(which(countOverlaps(Gr_loop_HFF_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_HFF_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                    1,
                    0)
  
  temp_K562 = ifelse(length(intersect(which(countOverlaps(Gr_loop_K562_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_K562_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                     1,
                     0)
  
  temp_H1_ES_MicroC = ifelse(length(intersect(which(countOverlaps(Gr_loop_H1_ES_MicroC_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_H1_ES_MicroC_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                             1,
                             0)
  
  temp_H1_endoderm_MicroC = ifelse(length(intersect(which(countOverlaps(Gr_loop_H1_endoderm_MicroC_1, Gr_loop_union_10000_ALL_1[i])>0), which(countOverlaps(Gr_loop_H1_endoderm_MicroC_2, Gr_loop_union_10000_ALL_2[i])>0))) >= 1,
                                   1,
                                   0)
  
  out_vect = as.numeric(c(temp_control,
                          temp_FL,
                          temp_NH4OAc,
                          temp_RNase,
                          temp_HFF,
                          temp_K562,
                          temp_H1_ES_MicroC,
                          temp_H1_endoderm_MicroC))
  return(out_vect)
}
  

temp_list = mclapply(c(1:length(Gr_loop_union_10000_ALL_1)), make_upset_list_ALL, mc.cores = 40)
df_upset_loops_10000_ALL = data.frame(do.call("rbind", temp_list))
colnames(df_upset_loops_10000_ALL) = c("Control","FL","NH4OAc","RNase","HFF","K562","H1_ES_MicroC","H1_endoderm_MicroC")
df_upset_loops_10000_ALL_full = cbind(data.frame(Gr_loop_union_10000_ALL_1)[,c(1,2,3,4)], data.frame(Gr_loop_union_10000_ALL_2)[,c(1,2,3,4)], df_upset_loops_10000_ALL)
colnames(df_upset_loops_10000_ALL_full)[1:8] = c("chr1","start1","end1","width1","chr2","start2","end2","width2")
write.table(df_upset_loops_10000_ALL_full, "/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_loops_10000_ALL_full.txt", row.names = F, col.names = T, sep = "\t", quote = F)


png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_loops_10000_ALL_full.png", width = 9, height = 7, units = "in", res = 400)
upset(df_upset_loops_10000_ALL_full, 
      sets = c("H1_endoderm_MicroC","H1_ES_MicroC","K562","HFF","RNase","FL","NH4OAc","Control"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.5, 0),
      #text.scale = c(1.8, 1.8, 0, 0, 1.6, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()

df_upset_FL_vs_control = df_upset_loops_10000_ALL_full[which(df_upset_loops_10000_ALL_full$Control == 0 &
                                                               df_upset_loops_10000_ALL_full$FL == 1),]

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_FL_vs_control.png", width = 9, height = 7, units = "in", res = 400)
upset(df_upset_FL_vs_control, 
      sets = c("H1_endoderm_MicroC","H1_ES_MicroC","K562","HFF","RNase","FL","NH4OAc"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.5, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()

df_upset_RNase_vs_control = df_upset_loops_10000_ALL_full[which(df_upset_loops_10000_ALL_full$Control == 0 &
                                                               df_upset_loops_10000_ALL_full$RNase == 1),]

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_RNase_vs_control.png", width = 9, height = 7, units = "in", res = 400)
upset(df_upset_RNase_vs_control, 
      sets = c("H1_endoderm_MicroC","H1_ES_MicroC","K562","HFF","RNase","FL","NH4OAc"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.5, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()

df_upset_FL_and_RNase_vs_control = df_upset_loops_10000_ALL_full[which(df_upset_loops_10000_ALL_full$Control == 0 &
                                                                  (df_upset_loops_10000_ALL_full$FL == 1 &
                                                                  df_upset_loops_10000_ALL_full$RNase == 1)),]

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_FL_and_RNase_vs_control.png", width = 9, height = 7, units = "in", res = 400)
upset(df_upset_FL_and_RNase_vs_control, 
      sets = c("H1_endoderm_MicroC","H1_ES_MicroC","K562","HFF","RNase","FL","NH4OAc"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.5, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()


df_upset_FL_or_RNase_vs_control = df_upset_loops_10000_ALL_full[which(df_upset_loops_10000_ALL_full$Control == 0 &
                                                                         (df_upset_loops_10000_ALL_full$FL == 1 |
                                                                            df_upset_loops_10000_ALL_full$RNase == 1)),]

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_FL_or_RNase_vs_control.png", width = 9, height = 7, units = "in", res = 400)
upset(df_upset_FL_or_RNase_vs_control, 
      sets = c("H1_endoderm_MicroC","H1_ES_MicroC","K562","HFF","RNase","FL","NH4OAc"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.5, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()



df_upset_FL_or_RNase_vs_control_no_ES_micro = df_upset_loops_10000_ALL_full[which(df_upset_loops_10000_ALL_full$Control == 0 &
                                                                        (df_upset_loops_10000_ALL_full$FL == 1 |
                                                                           df_upset_loops_10000_ALL_full$RNase == 1) &
                                                                          df_upset_loops_10000_ALL_full$H1_ES_MicroC == 0),]

png("/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/df_upset_FL_or_RNase_vs_control_no_ES_micro.png", width = 9, height = 7, units = "in", res = 400)
upset(df_upset_FL_or_RNase_vs_control_no_ES_micro, 
      sets = c("H1_endoderm_MicroC","K562","HFF","RNase","FL","NH4OAc"),
      sets.x.label = "Number of Loops",
      order.by = "freq",
      keep.order = T,
      number.angles = 45,
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.5, 0),
      # y-label
      # y-ticks
      # number of loops
      # loops axis
      # samples
      # bar labels
      point.size = 2, line.size = 0.5)
dev.off()



#######################################

extract_back_loops_from_upset <- function(index_loop,
                                          Gr_upset_1,
                                          Gr_upset_2,
                                          Gr_loop_1,
                                          Gr_loop_2)
{
  overlap1 = countOverlaps(Gr_loop_1[index_loop], Gr_upset_1)>0
  overlap2 = countOverlaps(Gr_loop_2[index_loop], Gr_upset_2)>0
  
  if (overlap1 & overlap2){
    return(index_loop)
  } else {
    return()
  }
}


emerging_loop_statistics <- function(loop_emerging,
                                     loop_not_emerging,
                                     loop_emerging_label,
                                     loop_not_emerging_label,
                                     output_path){
  
  out_stats_list = list()
  my_comparisons = list(c(loop_emerging_label, loop_not_emerging_label))
  label_size = 8
  
  ### CTCF
  out_stats_list[["CTCF"]] = t.test(loop_emerging$ctcf_both, loop_not_emerging$ctcf_both)
  
  df = data.frame(rbind(cbind(loop_emerging_label, loop_emerging$ctcf_both), cbind(loop_not_emerging_label, loop_not_emerging$ctcf_both)))
  colnames(df) = c("sample","value")
  df$value = as.numeric(as.character(df$value))
  df$sample = factor(df$sample, levels = c(loop_emerging_label,loop_not_emerging_label))
  
  p1 <- ggplot(df, aes(x=sample, y = value, fill=sample)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.2, size = jitter_size) +
    xlab('')+
    ylab("Number of CTCF peaks in loop anchors") +
    ylim(c(0, max(df$value)*1.2)) +
    stat_compare_means(comparisons = my_comparisons, method='t.test', label='p.signif', label.x = 1.5, label.y=max(df$value)*1.1, size=label_size-4) +
    theme_bw() +
    theme(text = element_text(size=label_size),
          plot.title = element_text(size=label_size),
          axis.title.x = element_text(size=label_size),
          axis.title.y = element_text(size=label_size),
          legend.position = "none")
  
  ### Density ratio HiC
  a = loop_emerging[which(!is.na(loop_emerging$density_ratio_hic) &
                            !is.infinite(loop_emerging$density_ratio_hic)),"density_ratio_hic"]
  b = loop_not_emerging[which(!is.na(loop_not_emerging$density_ratio_hic) &
                                !is.infinite(loop_not_emerging$density_ratio_hic)),"density_ratio_hic"]
  
  out_stats_list[["density_ratio_hic"]] = t.test(a, b)
  
  df = data.frame(rbind(cbind(loop_emerging_label, a), cbind(loop_not_emerging_label, b)))
  colnames(df) = c("sample","value")
  df$value = as.numeric(as.character(df$value))
  df$sample = factor(df$sample, levels = c(loop_emerging_label,loop_not_emerging_label))
  
  p2 <- ggplot(df, aes(x=sample, y = value, fill=sample)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.2, size = jitter_size) +
    xlab('')+
    ylab("HiC density ratio") +
    ylim(c(0, max(df$value)*1.2)) +
    stat_compare_means(comparisons = my_comparisons, method='t.test', label='p.signif', label.x = 1.5, label.y=max(df$value)*1.1, size=label_size-4) +
    theme_bw() +
    theme(text = element_text(size=label_size),
          plot.title = element_text(size=label_size),
          axis.title.x = element_text(size=label_size),
          axis.title.y = element_text(size=label_size),
          legend.position = "none")
  
  
  ### Density ratio iMARGI
  a = loop_emerging[which(!is.na(loop_emerging$density_ratio) &
                            !is.infinite(loop_emerging$density_ratio)),"density_ratio"]
  b = loop_not_emerging[which(!is.na(loop_not_emerging$density_ratio) &
                                !is.infinite(loop_not_emerging$density_ratio)),"density_ratio"]
  
  out_stats_list[["density_ratio_margi"]] = t.test(a, b)
  
  df = data.frame(rbind(cbind(loop_emerging_label, a), cbind(loop_not_emerging_label, b)))
  colnames(df) = c("sample","value")
  df$value = as.numeric(as.character(df$value))
  df$sample = factor(df$sample, levels = c(loop_emerging_label,loop_not_emerging_label))
  
  p3 <- ggplot(df, aes(x=sample, y = value, fill=sample)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    #geom_jitter(width = 0.2, size = jitter_size) +
    xlab('')+
    ylab("iMARGI density ratio") +
    ylim(c(0, max(df$value)*1.2)) +
    stat_compare_means(comparisons = my_comparisons, method='t.test', label='p.signif', label.x = 1.5, label.y=max(df$value)*1.1, size=label_size-4) +
    theme_bw() +
    theme(text = element_text(size=label_size),
          plot.title = element_text(size=label_size),
          axis.title.x = element_text(size=label_size),
          axis.title.y = element_text(size=label_size),
          legend.position = "none")
  
  ### Loop domain
  a = nrow(loop_emerging[which(loop_emerging$loop_domain == "loop_domain"),]) / nrow(loop_emerging) * 100
  b = nrow(loop_emerging[which(loop_emerging$loop_domain != "loop_domain"),]) / nrow(loop_emerging) * 100
  c = nrow(loop_not_emerging[which(loop_not_emerging$loop_domain == "loop_domain"),]) / nrow(loop_not_emerging) * 100
  d = nrow(loop_not_emerging[which(loop_not_emerging$loop_domain != "loop_domain"),]) / nrow(loop_not_emerging) * 100
  
  temp = as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T))
  chisq_temp = chisq.test(temp)
  chisq_temp$odds_ratio = (a*d) / (b*c)
  out_stats_list[["loop_domain"]] = chisq_temp
  
  df = data.frame(rbind(cbind(loop_emerging_label, "loop_domain", a),
                        cbind(loop_emerging_label, "not_loop_domain", b),
                        cbind(loop_not_emerging_label, "loop_domain", c),
                        cbind(loop_not_emerging_label, "not_loop_domain", d)))
  
  colnames(df) = c("sample","fill","value")
  df$value = as.numeric(as.character(df$value)) / 100
  df$sample = factor(df$sample, levels = c(loop_emerging_label,loop_not_emerging_label))
  df$fill = factor(df$fill, levels = c("loop_domain","not_loop_domain"))
  
  p4 <- ggplot(df, aes(x=sample, y = value, fill=fill)) +
    geom_bar(position="stack", stat="identity") +
    xlab('')+
    ylab("Number of loops / total loops") +
    theme_bw() +
    theme(text = element_text(size=label_size),
          plot.title = element_text(size=label_size),
          axis.title.x = element_text(size=label_size),
          axis.title.y = element_text(size=label_size),
          legend.position = "top",
          legend.title = element_blank())
  
  
  ### EP
  a = nrow(loop_emerging[which(!is.na(loop_emerging$EP)),]) / nrow(loop_emerging) * 100
  b = nrow(loop_emerging[which(is.na(loop_emerging$EP)),]) / nrow(loop_emerging) * 100
  c = nrow(loop_not_emerging[which(!is.na(loop_not_emerging$EP)),]) / nrow(loop_not_emerging) * 100
  d = nrow(loop_not_emerging[which(is.na(loop_not_emerging$EP)),]) / nrow(loop_not_emerging) * 100
  
  temp = as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T))
  chisq_temp = chisq.test(temp)
  chisq_temp$odds_ratio = (a*d) / (b*c)
  out_stats_list[["EP"]] = chisq_temp
  
  df = data.frame(rbind(cbind(loop_emerging_label, "EP", a),
                        cbind(loop_emerging_label, "not_EP", b),
                        cbind(loop_not_emerging_label, "EP", c),
                        cbind(loop_not_emerging_label, "not_EP", d)))
  
  colnames(df) = c("sample","fill","value")
  df$value = as.numeric(as.character(df$value)) / 100
  df$sample = factor(df$sample, levels = c(loop_emerging_label,loop_not_emerging_label))
  df$fill = factor(df$fill, levels = c("EP","not_EP"))
  
  p5 <- ggplot(df, aes(x=sample, y = value, fill=fill)) +
    geom_bar(position="stack", stat="identity") +
    xlab('')+
    ylab("Number of loops / total loops") +
    theme_bw() +
    theme(text = element_text(size=label_size),
          plot.title = element_text(size=label_size),
          axis.title.x = element_text(size=label_size),
          axis.title.y = element_text(size=label_size),
          legend.position = "top",
          legend.title = element_blank())
  
  p <- plot_grid(p1,p2,p3,p4,p5, nrow = 2)
  png(paste0(output_path, "/", loop_emerging_label, "_stats.png"), width = 9, height = 6, res = 200, units = "in")
  print(p)
  dev.off()
  
  return(out_stats_list)
}


###### Emerging loops in FL
temp=mclapply(1:length(Gr_loop_FL_1),
              extract_back_loops_from_upset,
              Gr_upset_1=Gr_loop_union_10000_1[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                             df_upset_H1_loops_10000_full$FL == 1),])],
              Gr_upset_2=Gr_loop_union_10000_2[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                             df_upset_H1_loops_10000_full$FL == 1),])],
              Gr_loop_1 = Gr_loop_FL_1,
              Gr_loop_2 = Gr_loop_FL_2,
              mc.cores=32)

loops_FL_no_control_index = do.call("rbind",temp)
loops_FL_no_control = loops_FL[loops_FL_no_control_index,1:6]
loops_FL_no_control$V1 = paste0("chr",loops_FL_no_control$V1)
loops_FL_no_control$V4 = paste0("chr",loops_FL_no_control$V4)
loops_FL_no_control = merge(loops_FL_no_control, loops_FL_full, by.x=paste0("V",c(1:6)), by.y=colnames(loops_FL_full)[1:6], all.x=T)
loops_FL_no_control = loops_FL_no_control[which(loops_FL_no_control$V1 != "chrY"),]

loops_FL_and_control = loops_FL_full[which(loops_FL_full$index %in% setdiff(loops_FL_full$index,loops_FL_no_control$index)),]

loops_FL_no_control_stats = emerging_loop_statistics(loop_emerging = loops_FL_no_control,
                         loop_not_emerging = loops_FL_and_control,
                         loop_emerging_label = "FL_no_control",
                         loop_not_emerging_label = "FL_and_control",
                         output_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/")


###### Emerging loops in RNase
temp=mclapply(1:length(Gr_loop_RNase_1),
              extract_back_loops_from_upset,
              Gr_upset_1=Gr_loop_union_10000_1[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                             df_upset_H1_loops_10000_full$RNase == 1),])],
              Gr_upset_2=Gr_loop_union_10000_2[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                             df_upset_H1_loops_10000_full$RNase == 1),])],
              Gr_loop_1 = Gr_loop_RNase_1,
              Gr_loop_2 = Gr_loop_RNase_2,
              mc.cores=32)

loops_RNase_no_control_index = do.call("rbind",temp)
loops_RNase_no_control = loops_RNase[loops_RNase_no_control_index,1:6]
loops_RNase_no_control$V1 = paste0("chr",loops_RNase_no_control$V1)
loops_RNase_no_control$V4 = paste0("chr",loops_RNase_no_control$V4)
loops_RNase_no_control = merge(loops_RNase_no_control, loops_RNase_full, by.x=paste0("V",c(1:6)), by.y=colnames(loops_RNase_full)[1:6], all.x=T)
loops_RNase_no_control = loops_RNase_no_control[which(loops_RNase_no_control$V1 != "chrY"),]

loops_RNase_and_control = loops_RNase_full[which(loops_RNase_full$index %in% setdiff(loops_RNase_full$index,loops_RNase_no_control$index)),]

loops_RNase_no_control_stats = emerging_loop_statistics(loop_emerging = loops_RNase_no_control,
                                                        loop_not_emerging = loops_RNase_and_control,
                                                        loop_emerging_label = "RNase_no_control",
                                                        loop_not_emerging_label = "RNase_and_control",
                                                        output_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/")


###### Emerging loops in FL and RNase
temp_FL=mclapply(1:length(Gr_loop_FL_1),
              extract_back_loops_from_upset,
              Gr_upset_1=Gr_loop_union_10000_1[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                             df_upset_H1_loops_10000_full$FL == 1 &
                                                                                             df_upset_H1_loops_10000_full$RNase == 1),])],
              Gr_upset_2=Gr_loop_union_10000_2[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                             df_upset_H1_loops_10000_full$FL == 1 &
                                                                                             df_upset_H1_loops_10000_full$RNase == 1),])],
              Gr_loop_1 = Gr_loop_FL_1,
              Gr_loop_2 = Gr_loop_FL_2,
              mc.cores=32)

temp_RNase=mclapply(1:length(Gr_loop_RNase_1),
                 extract_back_loops_from_upset,
                 Gr_upset_1=Gr_loop_union_10000_1[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                                df_upset_H1_loops_10000_full$FL == 1 &
                                                                                                df_upset_H1_loops_10000_full$RNase == 1),])],
                 Gr_upset_2=Gr_loop_union_10000_2[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                                df_upset_H1_loops_10000_full$FL == 1 &
                                                                                                df_upset_H1_loops_10000_full$RNase == 1),])],
                 Gr_loop_1 = Gr_loop_RNase_1,
                 Gr_loop_2 = Gr_loop_RNase_2,
                 mc.cores=32)

loops_FL_RNase_no_control_index_AND = do.call("rbind",temp_FL)
loops_FL_RNase_no_control = loops_FL[loops_FL_RNase_no_control_index_AND,1:6]
loops_FL_RNase_no_control$V1 = paste0("chr",loops_FL_RNase_no_control$V1)
loops_FL_RNase_no_control$V4 = paste0("chr",loops_FL_RNase_no_control$V4)
loops_FL_RNase_no_control = merge(loops_FL_RNase_no_control, loops_FL_full, by.x=paste0("V",c(1:6)), by.y=colnames(loops_FL_full)[1:6], all.x=T)
loops_FL_RNase_no_control = loops_FL_RNase_no_control[which(loops_FL_RNase_no_control$V1 != "chrY"),]

loops_RNase_FL_no_control_index_AND = do.call("rbind",temp_RNase)
loops_RNase_FL_no_control = loops_RNase[loops_RNase_FL_no_control_index_AND,1:6]
loops_RNase_FL_no_control$V1 = paste0("chr",loops_RNase_FL_no_control$V1)
loops_RNase_FL_no_control$V4 = paste0("chr",loops_RNase_FL_no_control$V4)
loops_RNase_FL_no_control = merge(loops_RNase_FL_no_control, loops_RNase_full, by.x=paste0("V",c(1:6)), by.y=colnames(loops_RNase_full)[1:6], all.x=T)
loops_RNase_FL_no_control = loops_RNase_FL_no_control[which(loops_RNase_FL_no_control$V1 != "chrY"),]

loops_FL_and_RNase_no_control = rbind(loops_FL_RNase_no_control, loops_RNase_FL_no_control)

loops_FL_RNase_and_control = loops_FL_full[which(loops_FL_full$index %in% setdiff(loops_FL_full$index,loops_FL_RNase_no_control$index)),]
loops_RNase_FL_and_control = loops_RNase_full[which(loops_RNase_full$index %in% setdiff(loops_RNase_full$index,loops_RNase_FL_no_control$index)),]
loops_FL_and_RNase_and_control = rbind(loops_FL_RNase_and_control, loops_RNase_FL_and_control)

loops_FL_and_RNase_no_control_stats = emerging_loop_statistics(loop_emerging = loops_FL_and_RNase_no_control,
                                                        loop_not_emerging = loops_FL_and_RNase_and_control,
                                                        loop_emerging_label = "FL_and_RNase_no_control",
                                                        loop_not_emerging_label = "FL_and_RNase_and_control",
                                                        output_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/")


###### Emerging loops in FL or RNase
temp_FL=mclapply(1:length(Gr_loop_FL_1),
                 extract_back_loops_from_upset,
                 Gr_upset_1=Gr_loop_union_10000_1[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                                (df_upset_H1_loops_10000_full$FL == 1 |
                                                                                                   df_upset_H1_loops_10000_full$RNase == 1)),])],
                 Gr_upset_2=Gr_loop_union_10000_2[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                                (df_upset_H1_loops_10000_full$FL == 1 |
                                                                                                   df_upset_H1_loops_10000_full$RNase == 1)),])],
                 Gr_loop_1 = Gr_loop_FL_1,
                 Gr_loop_2 = Gr_loop_FL_2,
                 mc.cores=32)


temp_RNase=mclapply(1:length(Gr_loop_RNase_1),
                    extract_back_loops_from_upset,
                    Gr_upset_1=Gr_loop_union_10000_1[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                                   (df_upset_H1_loops_10000_full$FL == 1 |
                                                                                                      df_upset_H1_loops_10000_full$RNase == 1)),])],
                    Gr_upset_2=Gr_loop_union_10000_2[rownames(df_upset_H1_loops_10000_full[which(df_upset_H1_loops_10000_full$Control == 0 & 
                                                                                                   (df_upset_H1_loops_10000_full$FL == 1 |
                                                                                                      df_upset_H1_loops_10000_full$RNase == 1)),])],
                    Gr_loop_1 = Gr_loop_RNase_1,
                    Gr_loop_2 = Gr_loop_RNase_2,
                    mc.cores=32)

loops_FL_RNase_no_control_index_OR = do.call("rbind",temp_FL)
loops_FL_RNase_no_control = loops_FL[loops_FL_RNase_no_control_index_OR,1:6]
loops_FL_RNase_no_control$V1 = paste0("chr",loops_FL_RNase_no_control$V1)
loops_FL_RNase_no_control$V4 = paste0("chr",loops_FL_RNase_no_control$V4)
loops_FL_RNase_no_control = merge(loops_FL_RNase_no_control, loops_FL_full, by.x=paste0("V",c(1:6)), by.y=colnames(loops_FL_full)[1:6], all.x=T)
loops_FL_RNase_no_control = loops_FL_RNase_no_control[which(loops_FL_RNase_no_control$V1 != "chrY"),]

loops_RNase_FL_no_control_index_OR = do.call("rbind",temp_RNase)
loops_RNase_FL_no_control = loops_RNase[loops_RNase_FL_no_control_index_OR,1:6]
loops_RNase_FL_no_control$V1 = paste0("chr",loops_RNase_FL_no_control$V1)
loops_RNase_FL_no_control$V4 = paste0("chr",loops_RNase_FL_no_control$V4)
loops_RNase_FL_no_control = merge(loops_RNase_FL_no_control, loops_RNase_full, by.x=paste0("V",c(1:6)), by.y=colnames(loops_RNase_full)[1:6], all.x=T)
loops_RNase_FL_no_control = loops_RNase_FL_no_control[which(loops_RNase_FL_no_control$V1 != "chrY"),]

loops_FL_or_RNase_no_control = rbind(loops_FL_RNase_no_control, loops_RNase_FL_no_control)

loops_FL_RNase_and_control = loops_FL_full[which(loops_FL_full$index %in% setdiff(loops_FL_full$index,loops_FL_RNase_no_control$index)),]
loops_RNase_FL_and_control = loops_RNase_full[which(loops_RNase_full$index %in% setdiff(loops_RNase_full$index,loops_RNase_FL_no_control$index)),]
loops_FL_or_RNase_and_control = rbind(loops_FL_RNase_and_control, loops_RNase_FL_and_control)

loops_FL_or_RNase_no_control_stats = emerging_loop_statistics(loop_emerging = loops_FL_or_RNase_no_control,
                                                               loop_not_emerging = loops_FL_or_RNase_and_control,
                                                               loop_emerging_label = "FL_or_RNase_no_control",
                                                               loop_not_emerging_label = "FL_or_RNase_and_control",
                                                               output_path = "/dataOS/rcalandrelli/phase_separation/HiC/loops/differential_loop_analysis/")

