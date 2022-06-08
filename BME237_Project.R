setwd('/Users/wenyuliang/Documents/BME237/counts/')
library(tximportData)
library(tximport)
library(DESeq2)
library(tidyverse)
library(janitor)
library(ggthemes)
library(pheatmap)
#-----------------------------------------mouse--------------------------------------------#
mouse = c(paste0("SRR678110",6:9),paste0("SRR67811",10:11))
mouse <- file.path('.', mouse, "abundance.h5")
mouseSample <- tximport(mouse, type = "kallisto", txOut = TRUE)
#--------------------------------------filter----------------------------------------------#
mouseFilter <- data.frame(j= rowSums(mouseSample$counts)>5)$j
mouseSample$counts <- mouseSample$counts[mouseFilter,]
mouseSample$abundance <- mouseSample$abundance[mouseFilter,]
mouseSample$length <- mouseSample$length[mouseFilter,]
#--------------------------------------filter----------------------------------------------#
colnames(mouseSample$counts) <- c(paste0('CTRL',1:3),paste0('KI',1:3))
mouseTable <- data.frame(condition = factor(rep(c("CTRL", "KI"), each = 3)))
rownames(mouseTable) <- colnames(mouseSample$counts)
ddsmouse <- DESeqDataSetFromTximport(mouseSample, mouseTable, ~condition)

ddsmouse <- DESeq(ddsmouse)
resmouse <- results(ddsmouse)
write.table(mouseSample$counts, 'mouseCount',sep = '\t')
mouseRes <- resmouse %>% as.data.frame() %>% filter(padj < 0.05) %>% rownames_to_column('transcript') 
mouseMart <- read.table('mouseTrans.txt', fill = TRUE, sep = '\t',header = TRUE) 
mouseRes <- merge(mouseRes, mouseMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') 
mouseD <- c('Foxg1', 'Id4','Fezf2', 'Sox3', 'Six3')
mouseU <- c('Cntn2', 'Nefl', 'Gap43', 'Sox10')


# mouseAll <- resmouse %>% as.data.frame()%>% rownames_to_column('transcript')
# mouseAll <- merge(mouseAll, mouseMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') 
# mouseAll %>% filter(Gene.name %in% mouseD) %>% view()
#-----------------------------------------mouse--------------------------------------------#


#-----------------------------------------human--------------------------------------------#
human = paste0("SRR678110",0:5)
human <- file.path('.', human, "abundance.h5")
humanSample <- tximport(human, type = "kallisto", txOut = TRUE)
humanFilter <- data.frame(j= rowSums(humanSample$counts)>5)$j
humanSample$counts <- humanSample$counts[humanFilter,]
humanSample$abundance <- humanSample$abundance[humanFilter,]
humanSample$length <- humanSample$length[humanFilter,]
write.table(humanSample$counts, 'humanCount',sep = '\t')
#humanCounts <- humanSample$counts
colnames(humanSample$counts) <- c(paste0('WT',1:3),paste0('KO',1:3))
humanTable <- data.frame(condition = factor(rep(c("WT", "KO"), each = 3)))
rownames(humanTable) <- colnames(humanSample$counts)
ddshuman <- DESeqDataSetFromTximport(humanSample, humanTable, ~condition)
ddshuman <- DESeq(ddshuman)
reshuman <- results(ddshuman)
humanRes <- reshuman %>% as.data.frame() %>% filter(padj < 0.05) %>% rownames_to_column('transcript') 
humanMart <- read.table('humanTrans.txt', fill = TRUE, sep = '\t', header = TRUE) 
humanRes <- merge(humanRes, humanMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') 
#-----------------------------------------human--------------------------------------------#
humanCounts <-humanSample$counts %>% as.data.frame() %>% mutate(sum = apply(humanSample$counts, 1, sum)) %>% view()
NOTCH2NL_human <- humanCounts %>% as.data.frame()  %>% rownames_to_column('transcript') %>% merge(humanMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') %>% filter(Gene.name%in%c('NOTCH2NLA','NOTCH2NLB','NOTCH2NLC'))
NOTCH2NL_human <- NOTCH2NL_human %>% group_by(Gene.name) %>% summarise_if(is.numeric, sum) %>% select(-sum) %>% t() %>% 
  janitor::row_to_names(1) %>% as.data.frame() %>% rownames_to_column('sample') %>% mutate(across(NOTCH2NLA:NOTCH2NLC, as.numeric))
NOTCH2NL_human$Condition <- gsub('[123]','',NOTCH2NL_human$sample)
ggplot(data = NOTCH2NL_human, aes(x=sample, y=NOTCH2NLA, color = Condition))+geom_point(size = 3)+ theme_light() + scale_color_stata()+
  ggtitle('Human')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('HumanA.png',dpi = 300)
dev.off()

ggplot(data = NOTCH2NL_human, aes(x=sample, y=NOTCH2NLB, color = Condition))+geom_point(size = 3)+ theme_light() + scale_color_stata()+
  ggtitle('Human')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('HumanB.png',dpi = 300)
dev.off()

ggplot(data = NOTCH2NL_human, aes(x=sample, y=NOTCH2NLC, color = Condition))+geom_point(size = 3)+ theme_light() + scale_color_stata()+
  ggtitle('Human')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('HumanC.png',dpi = 300)
dev.off()

# write.table(humanCounts %>% filter(sum>0) %>% select(-sum),'humanCount',  sep = '\t')
# NOTCH2NL <- humanCounts %>% as.data.frame()  %>% rownames_to_column('transcript') %>% merge(humanMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') %>% filter(Gene.name%in%c('NOTCH2NLA','NOTCH2NLB','NOTCH2NLC'))
# NOTCH2NL %>% group_by(Gene.name) %>% summarise_if(is.numeric, sum) %>% select(-sum) 




#-----------------------------------------benchmark--------------------------------------------#
benchmark = c(paste0("SRR678110",6:9),paste0("SRR67811",10:11))
benchmark <- file.path('./benchmark/benchmark/', benchmark, "abundance.h5")
benchmarkSample <- tximport(benchmark, type = "kallisto", txOut = TRUE)
colnames(benchmarkSample$counts) <- c(paste0('CTRL',1:3),paste0('KI',1:3))
benchmarkTable <- data.frame(condition = factor(rep(c("CTRL", "KI"), each = 3)))
rownames(benchmarkTable) <- colnames(benchmarkSample$counts)
ddsbenchmark <- DESeqDataSetFromTximport(benchmarkSample, benchmarkTable, ~condition)
ddsbenchmark <- DESeq(ddsbenchmark)
resbenchmark <- results(ddsbenchmark)
benchmarkRes <- resbenchmark %>% as.data.frame()  %>% rownames_to_column('transcript')
benchmarkMart <- read.table('humanTrans.txt', fill = TRUE, sep = '\t', header = TRUE) 
benchmarkRes <- merge(benchmarkRes, benchmarkMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version')
benchmarkRes %>% filter(Gene.name%in%c('NOTCH2NLA','NOTCH2NLB','NOTCH2NLC')) %>% view()
benchmarkCounts <-benchmarkSample$counts %>% as.data.frame() %>% mutate(sum = apply(benchmarkSample$counts, 1, sum))
NOTCH2NL_benchmark <- benchmarkCounts %>% as.data.frame()  %>% rownames_to_column('transcript') %>% merge(benchmarkMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') %>% filter(Gene.name%in%c('NOTCH2NLA','NOTCH2NLB','NOTCH2NLC'))
NOTCH2NL_benchmark <- NOTCH2NL_benchmark %>% group_by(Gene.name) %>% summarise_if(is.numeric, sum) %>% select(-sum) %>% t() %>% 
  janitor::row_to_names(1) %>% as.data.frame() %>% rownames_to_column('sample') %>% mutate(across(NOTCH2NLA:NOTCH2NLC, as.numeric))
NOTCH2NL_benchmark$Condition <- gsub('[123]','',NOTCH2NL_benchmark$sample)
ggplot(data = NOTCH2NL_benchmark, aes(x=sample, y=NOTCH2NLA, color = Condition))+geom_point(size = 3)+ theme_light() + scale_color_stata()+
  ggtitle('Mouse')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('MouseA.png',dpi = 300)
dev.off()

ggplot(data = NOTCH2NL_benchmark, aes(x=sample, y=NOTCH2NLB, color = Condition))+geom_point(size = 3)+ theme_light() + scale_color_stata()+
  ggtitle('Mouse')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('MouseB.png',dpi = 300)
dev.off()

ggplot(data = NOTCH2NL_benchmark, aes(x=sample, y=NOTCH2NLC, color = Condition))+geom_point(size = 3)+ theme_light() + scale_color_stata()+
  ggtitle('Mouse')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('MouseC.png',dpi = 300)
dev.off()


#-----------------------------------------benchmark--------------------------------------------#


k_human = kmeans(t(humanSample$counts), centers = 2)
k_mouse = kmeans(t(mouseSample$counts), centers = 2, iter.max = 50)

library(FactoMineR)
library(factoextra)
humanCounts <- t(humanSample$counts) %>% as.data.frame() %>% rownames_to_column('group')
humanCounts$group <- c(rep(c("KO", "WT"), each = 3))
pcaHuman <- PCA(humanCounts[,2:ncol(humanCounts)], graph = FALSE)
fviz_pca_ind(pcaHuman,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = humanCounts$group, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


#-----------------------------------------QC--------------------------------------------#
library(DEGreport)
mouseQC <- counts(ddsmouse, normalized = TRUE)
mouseDesign <- as.data.frame(colData(ddsmouse))
degCheckFactors(mouseQC) + ggtitle('Mouse')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('MouseLib.png', dpi = 300)
dev.off()
humanQC <- counts(ddshuman, normalized = TRUE)
humanDesign <- as.data.frame(colData(ddshuman))
degCheckFactors(humanQC) + ggtitle('human')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('HumanLib.png', dpi = 300)
dev.off()

png('MouseMV.png',width = 4, height = 4, units = 'in',res = 300)
degQC(mouseQC, mouseDesign[['condition']], pvalue = resmouse[['pvalue']]) 
dev.off()

png('HumanMV.png',width = 4, height = 4, units = 'in',res = 300)
degQC(humanQC, humanDesign[['condition']], pvalue = reshuman[['pvalue']]) 
dev.off()
#-----------------------------------------QC--------------------------------------------#


#-----------------------------------------GO--------------------------------------------#
humanRes %>% filter(log2FoldChange>0 & Gene.name !='') %>% arrange(padj) %>% top_n(250) %>% 
  select(Gene.name) %>% write.table('humanUp',sep='\t',row.names = FALSE,col.names = FALSE, quote = FALSE)
humanRes %>% filter(log2FoldChange<0 & Gene.name !='') %>% arrange(padj) %>% top_n(250) %>% 
  select(Gene.name) %>% write.table('humanDown',sep='\t',row.names = FALSE,col.names = FALSE, quote = FALSE)
humanRes %>% filter(Gene.name !='') %>% arrange(padj) %>% top_n(500) %>% 
  select(Gene.name) %>% write.table('human500',sep='\t',row.names = FALSE,col.names = FALSE, quote = FALSE)

write.table(unique(humanRes$Gene.name), 'humanID',sep='\t',row.names = FALSE,col.names = FALSE, quote = FALSE)
write.table(unique(mouseRes$Gene.name), 'mouseID',sep='\t',row.names = FALSE,col.names = FALSE, quote = FALSE)
mouseRes %>% filter(log2FoldChange<0)%>% 
  select(Gene.name) %>% write.table('mouseUp',sep='\t',row.names = FALSE,col.names = FALSE, quote = FALSE)

#------------------------------enrichment-------------------------------------------#
library(RColorBrewer)
mouseGo <- read.table('mouseNeuronGene.txt',fill = TRUE,sep = '\t')
mouseGene <- mouseGo$V2
#mouseGene <- c('Nnat','Ptprz1','Foxg1','Wnt5a','Fezf2', 'Dmrt3', 'Wnt7b')
mousePaint <- mouseQC %>% as.data.frame() %>% rownames_to_column('transcript') %>%
  merge(mouseMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') %>% 
  filter(Gene.name%in%mouseGene) %>% group_by(Gene.name) %>% summarise_if(is.numeric,sum) %>% column_to_rownames('Gene.name')
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
png('mouseHeatmap.png',width = 4, height = 4, units = 'in',res = 300)
pheatmap(mousePaint,scale = "row",cluster_rows = FALSE,
         cluster_cols = FALSE,color = coul,border_color = NA, main = 'ND Associated Genes')
dev.off()
#heatmap(as.matrix(mousePaint),Rowv = NA,Colv = NA,col = coul)
# mouseGGplot <- mousePaint %>% rownames_to_column('Gene') %>% pivot_longer(!Gene, names_to = 'Sample', values_to = 'count')
# library(viridis)
# ggplot(mouseGGplot, aes(Sample,Gene,fill = count)) + 
#   geom_tile() +
#   scale_fill_viridis(discrete=FALSE)

allMouse <- read.table('allMouseNeuron.txt',sep = '\t', fill = TRUE)
humanUD <- c('BCL11B', 'CTIP2', 'DLX1', 'SEMA3A', 'UNC5D', 'FGFR2', 'WNT7B', 'EPHA4', 'SLC1A3', 'TNC', 'UNC5D','DLX1', 'SLITR', 'ISL1', 'PBX3')
humanGo <- read.table('HumanNeuronGene.txt',fill = TRUE,sep = '\t')
humanUp <- humanRes %>% filter(Gene.name%in%humanGo$V2 & log2FoldChange<0) 
humanD <- humanRes %>% filter(Gene.name%in%humanGo$V2 & log2FoldChange>0)  
humanPaintUp <- humanQC %>% as.data.frame() %>% rownames_to_column('transcript') %>%
  merge(humanMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') %>% 
  filter(Gene.name%in%humanUp$Gene.name) %>% group_by(Gene.name) %>% summarise_if(is.numeric,sum) %>% column_to_rownames('Gene.name')
humanPaintD <- humanQC %>% as.data.frame() %>% rownames_to_column('transcript') %>%
  merge(humanMart, by.x = 'transcript', by.y = 'Transcript.stable.ID.version') %>% 
  filter(Gene.name%in%humanD$Gene.name) %>% group_by(Gene.name) %>% summarise_if(is.numeric,sum) %>% column_to_rownames('Gene.name')
png('humanUp.png',width = 4, height = 4, units = 'in',res = 300)
pheatmap(humanPaintUp,scale = "row",cluster_rows = FALSE,
         cluster_cols = FALSE,color = coul,border_color = NA, main = 'Up regulated Genes')
dev.off()
png('humanD.png',width = 4, height = 4, units = 'in',res = 300)
pheatmap(humanPaintD,scale = "row",cluster_rows = FALSE,
         cluster_cols = FALSE,color = coul,border_color = NA, main = 'Down regulated Genes')
dev.off()
#------------------------------enrichment-------------------------------------------#
