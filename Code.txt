#PCoA
library(vegan)
library(ggplot2)
library(plyr)
fun_otu <- read.delim('otu_chouping.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
fun_otu <- data.frame(t(fun_otu))
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = F)
fun_distance <- vegdist(fun_otu, method = 'bray')
fun_pcoa <- cmdscale(fun_distance, k = (nrow(fun_otu) - 1), eig = T)
ordiplot(scores(fun_pcoa)[ ,c(1, 2)], type = 't')
fun_pcoa$eig
fun_point <- data.frame(fun_pcoa$point)
fun_species <- wascores(fun_pcoa$points[,1:2], fun_otu)
write.csv(as.matrix(fun_distance), 'fun_distance.csv', quote = F)
write.csv(fun_point, 'fun_pcoa.sample.csv')
write.csv(fun_species, 'fun_pcoa.otu.csv')
fun_pcoa_eig <- (fun_pcoa$eig)[1:2] / sum(fun_pcoa$eig)
fun_sample_site <- data.frame({fun_pcoa$point})[1:2]
fun_sample_site$names <- rownames(fun_sample_site)
names(fun_sample_site)[1:2] <- c('PCoA1', 'PCoA2')
fun_sample_site <- merge(fun_sample_site, group, by = 'names', all.x = T)
write.csv(fun_sample_site, 'fun_sample_site.csv', quote = F)
fun_sample_site = read.csv("fun_sample_site.csv", header=T)
treatment_border <- ddply(fun_sample_site, 'treatment', function(df) df[chull(df[[3]], df[[4]]), ])
 p2<-ggplot(fun_sample_site, aes(x=PCoA1, y=PCoA2))+
   geom_polygon(data = treatment_border, aes(fill = treatment)) + 
   geom_point(aes(color = treatment), size = 3, alpha = 1) + 
   scale_color_manual(values = c('#C673FF', '#73D5FF', '#49C35A')) + 
   scale_fill_manual(values = c('#C673FF2E', '#73D5FF2E', '#49C35A2E')) +  
 # geom_point(aes(color=treatment),size=3,alpha=1)+
#geom_point(aes(color = treatment, shape = Season), size = 1.5, alpha = 0.8) + 
  # scale_color_manual(values = c("DodgerBlue","Coral","DodgerBlue4","Coral3"))+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
      panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) + 
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  labs(x = paste('PCoA1: ', round(100 * fun_pcoa_eig[1], 2), '%'), 
       y = paste('PCoA2: ', round(100 * fun_pcoa_eig[2], 2), '%')) 
p2


##PERMANOVA 
fun_adonis_result_otu=adonis2(fun_otu~group$treatment, data=fun_otu,permutations = 999,method="bray")
fun_adonis_result_otu


#Stacked bar chart
library(reshape2)	
library(ggplot2)	
library(dplyr)
phylum_top10 <- read.delim('phylum_top10.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
phylum_top10$phylum<- factor(rownames(phylum_top10), levels = rev(rownames(phylum_top10)))
phylum_top10 <- melt(phylum_top10, id = 'phylum')
p2 <- ggplot(phylum_top10, aes(variable, value, fill = phylum)) +
  geom_col(position = 'stack', width = 0.6) +
  # facet_wrap(~group, scales = 'free_x', ncol = 3) +
  scale_fill_manual(values =  rev(c( "#DC0000B2","#00A087B2", "lightskyblue","#F39B7FB2","#91D1C2B2",  "#8491B4B2",  "darkgoldenrod1", "#3C5488B2", "khaki2",
                                     "darkorange1", "#CC9933","#996666","#FFFF99","#FF9900","#3399FF","#FFCCCC","#CC6633","gray"))) +
  labs(x = 'Cultivar', y = 'Relative Abundance(%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))
p2


#Heatmap 
library("ggplot2")
library("reshape2")
library("pheatmap") 
test<-read.csv("F-H.csv",row.names=1,header=TRUE)
p<-pheatmap(test, scale = "row",cluster_rows = F,cluster_cols = F)
p


#Bar 
library(ggpubr)
library(dplyr)
library(ggplot2)
library(MetBrewer)
k<-met.brewer("Hokusai3",n=17) #Hokusai3
k
#df_res$t = factor(df_res$t, levels=c('ck','Rs1','Rs2','up1','up2',
 #                                    'down1','down2','st1','st2','ax1','ax2' )) 
setwd("D:/PhD/LSS")
df_res<-read.csv("FOC.csv")
ggplot() +
  geom_bar(df_res, mapping = aes(t, N, fill = t, color = t), stat = "identity", width = 0.7) +
  #geom_jitter(df, mapping = aes(t, N), position = position_jitter(0.1), size = 2, color = "gray") +
  geom_errorbar(
    aes(x = t, ymin = N - NSD, ymax = N + NSD),
    data = df_res, position = position_dodge(0.2), width = 0.2
  ) +
  scale_fill_manual(values = k) +  
  scale_color_manual(values = k) +  
  theme_bw() +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  theme_classic() +
  labs(y='Abundance of pathogen (lg copies g-1 soil)')+
  xlab(NULL) +
  coord_cartesian(ylim=c(4,7))+
  theme(plot.margin = unit(c(2, 2, 2, 1), "cm"))


#Box
library(ggbeeswarm)
k<-met.brewer("Johnson",n=28) #Hokusai3
k
a<-read.csv("FOC-N-YZ.csv",header=T)
p1<-ggplot(a,aes(x=t, y=Y, fill=t,color=t))+
  scale_fill_manual(values=(k)) + #scale_fill_manual(values=rev(k)) 
  scale_color_manual(values=(k))+
  #scale_fill_manual(values=c('#0a5c31', '#778f2f', '#efce6b', '#e69c3b'))+
  #scale_color_manual(values=c('#0a5c31', '#778f2f', '#efce6b', '#e69c3b'))+
  geom_quasirandom(alpha =0.8, #0.7, 
                   width = 0.4, 
                   size = 1.5) +
  geom_boxplot(size = 0.45, 
               width = 0.6, 
               #fill = "lightskyblue", 
               color = "black", 
               alpha =0.2, #0.2, 
               outlier.color = "black", 
               outlier.fill = "black", 
               outlier.shape = 19, 
               outlier.size = 3, 
               outlier.stroke = 0.5, 
               outlier.alpha = 0.4, 
               notch = F, 
               notchwidth = 0.3) +

  stat_boxplot(geom = "errorbar", 
               width = 0.2, 
               size = 0.9, 
               color = "black") +  #geom_smooth(method = "loess", aes(group=g))+
  theme_bw() +
  theme(
    legend.position='none',
)+
  theme_classic()+
  #ylim(0,100)+
  guides (fill=FALSE)+
  labs(y='F_all')+
  theme (axis.title.x = element_blank ()) +
  #xlab(NULL)
  #geom_hline(yintercept =0)+
  #theme_classic()+
  #coord_cartesian(ylim=c(-1,1))+
  theme(plot.margin = unit(c(1,1,2,1),"cm"))
p1
ggsave('FOC-N-YZ.pdf', p1, width = 10, height = 4.38)


# Line chart
library(MetBrewer)
library(ggpubr)
k<-met.brewer("Hokusai3",n=4)
k
data<-read.csv("supernatent.csv")
p1 <-ggplot(data, aes(x=time, y=Y, color=t)) +
  scale_color_manual(values = c(k)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  #scale_y_continuous(limits = c(0,0.6),breaks = seq(0,0.6,0.15),expand = c(0,0)) +
  #scale_x_continuous(limits = c(0,7.5),breaks = seq(0,7,1),expand = c(0,0)) +
  geom_errorbar(aes(ymin=Y-SD, ymax=Y+SD), color = "grey",width=0.2) +
  labs(title='',
       x='Time (h)',
       y='absorbance (OD600)') +
  theme_classic()+
  theme(plot.margin = unit(c(2,4,2,3),"cm"))
p1
ggsave('supernatent.pdf', width = 6, height = 4)


