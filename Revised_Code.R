#Part1=================================Main_celltype_ann========================
library(Seurat)
library(ggplot2)
library(dplyr)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
library(stringi)
library(stringr)
#===============================================================================
scRNA <- CreateSeuratObject(count, meta.data = scRNA@meta.data)
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
dim(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA)
load(file="cycle.rda")
scRNA <- CellCycleScoring(scRNA,g2m.features = g2m_genes,s.features = s_genes)
DimPlot(scRNA,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
#===============================================================================
scRNA <-ScaleData(scRNA, vars.to.regress = "Phase")
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
print(scRNA[["pca"]], dims = 1:5, nfeatures = 5)
#
DimHeatmap(scRNA, dims = 1:25, cells = 500, balanced = TRUE)
ElbowPlot(scRNA)
# 
scRNA <- FindNeighbors(object = scRNA, 
                       dims = 1:15)
#                                
scRNA <- FindClusters(object = scRNA,
                      resolution = c(0.2, 0.3, 0.4, 0.5, 0.6,0.7,0.8))
clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
scRNA <- RunUMAP(scRNA,reduction ="pca",dims = 1:15)
#
B_cell=c(6)
Plasma_cells=c(3)
TNKIL_cells=c(0)
Myeloids_cells=c(2,18)
Epi=c(1,7,9,14,5,4,17,16,19,8,15)
Endo=c(11)
Peri=c(12)
Fibro=c(10) 
DimPlot(scRNA,
        reduction = "umap",
        label = TRUE,
        label.size = 5,group.by ="celltype")
#
length(table(scRNA@meta.data$patient))
DimPlot(scRNA,
        reduction = "tsne",
        label = TRUE,
        label.size = 5,group.by ="patient")
library(RColorBrewer)
display.brewer.all(type="all")
display.brewer.pal(9, "RdYlBu")
a <- brewer.pal(11, "RdYlBu")
display.brewer.pal(11, "PRGn")
b <- brewer.pal(11, "PRGn")[c(1,2,3,4,5,7,8,9,10,11)]
display.brewer.pal(11, "PiYG")
c <- brewer.pal(11, "PiYG")[c(1,2,3,4,5)]
d <- c(a,b,c)
library(scales)
show_col(d,labels=F,ncol=1)
#===============================================================================
TSNE_coord <- scRNA@reductions$tsne@cell.embeddings
TSNE_coord <- as.data.frame(TSNE_coord)
metadata=scRNA@meta.data
metadata$CB <- row.names(metadata)
TSNE_coord$CB <- row.names(TSNE_coord)
plot_data <- merge(TSNE_coord,metadata,by="CB")
library(dplyr)
library(ggthemes)
library(ggsci)
library(ggplot2)
colnames(plot_data)
ggplot(data=plot_data, aes(x=tSNE_1, y=tSNE_2,color=patient)) + 
  geom_point(shape=16,alpha=0.8,size=0.5)+
  scale_color_manual(values=d)+
  #scale_color_manual(values=c("#D7301F","#B8E186","#0570B0"))+
  #scale_color_npg()+
  theme_few()
#
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}
library(Seurat)

table(scRNA@meta.data$celltype)
metadata=scRNA@meta.data

colnames(metadata)
table(metadata$Group)
meta.tb=metadata[,c("patient","celltype","Group")]
colnames(meta.tb) <- c("patient","meta.cluster","loc")
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             out.prefix=sprintf("%s.STARTRAC.dist",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)
#===============================================================================
#===============================================================================
#Part2 Malignant epithelial cell identification
library(copykat)
counts = as.matrix(Epi@assays$RNA@counts)
sc_cnv = copykat(rawmat = counts,ngene.chr = 5,sam.name = 'Malignant_epithelium',n.cores = 20)
#
import scanpy
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
from cnmf import cNMF
#

cnmf_obj = cNMF(output_dir=output_directory, name=run_name)

cnmf_obj.prepare(counts_fn=countfn, components=np.arange(2,15), n_iter=10, seed=123, num_highvar_genes=2000)

cnmf_obj.factorize(worker_i=0, total_workers=1)
cnmf_obj.combine()
cnmf_obj.k_selection_plot()
#
selected_K = 10
density_threshold = 2.00
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
#
density_threshold2 = 0.01
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold2, show_clustering=True, close_clustergram_fig=False)
usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold2)
#
#===============================================================================
#===============================================================================
#Part_3
library(AUCell)
library(clusterProfiler)
#
cells_rankings <- AUCell_buildRankings(Malignant_Epithelium@assays$RNA@data) 
Hallmarker <- read.gmt("h.all.v7.5.1.symbols.gmt") 
geneSets <- lapply(unique(Hallmarker$term), function(x){print(x);Hallmarker$gene[Hallmarker$term == x]})
names(geneSets) <- unique(Hallmarker$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
geneSet <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
Malignant_Epithelium$AUC  <- aucs

library(ggraph)
ggplot(data.frame(Malignant_Epithelium@meta.data, Malignant_Epithelium@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "TNFA_SIGNALING_VIA_NFKB")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))
#===============================================================================
#===============================================================================
#Part04
#
scRNA <- CreateSeuratObject(counts = Myeloids_cells@assays$RNA@counts, 
                            meta.data = scRNA@meta.data,
                            min.cells = 3, min.features = 200)

scRNA <- NormalizeData(scRNA)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
scRNA<- RunPCA(scRNA,verbose = F)
ElbowPlot(scRNA,ndims = 50)
DimHeatmap(scRNA,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:40
scRNA <- scRNA %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
DimPlot(scRNA,group.by = "RNA_snn_res.0.2",reduction = "umap",label = T)
DimPlot(scRNA,group.by = "clMidwayPr",reduction = "umap",label = T)
#===============================================================================
#===============================================================================
Monocyte <- subset(scRNA,clMidwayPr=="Mono")
Monocyte <- CreateSeuratObject(counts = Monocyte@assays$RNA@counts, 
                               meta.data = Monocyte@meta.data,
                               min.cells = 3, min.features = 200)

Monocyte <- NormalizeData(Monocyte)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
Monocyte<- RunPCA(Monocyte,verbose = F)
ElbowPlot(Monocyte,ndims = 50)
DimHeatmap(Monocyte,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:50
Monocyte <- Monocyte %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
Monocyte <- FindNeighbors(Monocyte, dims=pc.num)%>% FindClusters(resolution=c(0.1,0.2,0.3,0.4,0.5))
DimPlot(Monocyte,group.by = "RNA_snn_res.0.2",reduction = "umap",label = T)
#===============================================================================
Mono_CD14 <- c(1,2)
Mono_CD14_CD16 <- c(0,5)
Mono_CD16 <- c(6)
Mono_CD3E <- c(3)
Undefined <- c(4)
current.cluster.ids <- c(Mono_CD14,
                         Mono_CD14_CD16,
                         Mono_CD16,
                         Mono_CD3E,Undefined)

new.cluster.ids <- c(rep("Mono_CD14",length(Mono_CD14)),
                     rep("Mono_CD14_CD16",length(Mono_CD14_CD16)),
                     rep("Mono_CD16",length(Mono_CD16)),
                     rep("Mono_CD3E",length(Mono_CD3E)),
                     rep("Undefined",length(Undefined)))

Monocyte@meta.data$LV_03 <- plyr::mapvalues(x = as.character(Monocyte@meta.data$RNA_snn_res.0.2), 
                                            from = current.cluster.ids, to = new.cluster.ids)
#===============================================================================
table(scRNA@meta.data$clMidwayPr)
Macrophage <- subset(scRNA,clMidwayPr=="Macro")
Macrophage <- CreateSeuratObject(counts = Macrophage@assays$RNA@counts, 
                                 meta.data = Macrophage@meta.data,
                                 min.cells = 3, min.features = 200)

Macrophage <- NormalizeData(Macrophage)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
Macrophage<- RunPCA(Macrophage,verbose = F)
ElbowPlot(Macrophage,ndims = 50)
DimHeatmap(Macrophage,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:20
Macrophage <- Macrophage %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
Macrophage <- FindNeighbors(Macrophage, dims=pc.num)%>% FindClusters(resolution=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))
DimPlot(Macrophage,group.by = "RNA_snn_res.0.8",reduction = "umap",label = T)
DimPlot(Macrophage,group.by = "RNA_snn_res.0.2",reduction = "umap",label = T)
#===============================================================================
#标识基因热图帮助注释
Macro_FCN1 <- c("S100A9", "S100A8", "THBS1", "VCAN", "FCN1", "HLA-DQA1", "HLA-DQB1", "HSPA1A", "HSPA1B", "CXCR4", "FCGR2A", "FCGR3A")
Macro_NLPR3 <- c("NLRP3", "CCL20", "PTGS2", "EREG", "AQP9", "METRNL", "GPR132", "ANPEP", "PLAUR", "EREG", "VEGFA", "PTGS2")
Macro_PLTP <- c("KLF2", "KLF4", "EGR1", "CREM", "HIF1A","PLTP", "LYVE1", "IL10", "STAB1", "SEPP1", "SNX6", "CTSL", "THBS1", "MARCO", "CXCL3", "CALR", "CD163", "CD36")
Macro_IL1B <- c("HES1", "REL", "BTG2","IL1B", "C1QA", "C1QC", "C1QB", "CCL3L3", "HLA-DQA1", "HLA-DQA2","HLA-DQB1", "GPR183", "CD83", "CLEC10A", "MRC1")
Macro_C1QC <- c("MAFB", "ATF3", "MAF", "MEF2A", "PA2G4","C1QA","C1QB","C1QC", "MERTK", "FPR3", "TREM2", "MS4A4A", "SLCO2B1", "NRP1", "SLAMF8", "FCGR1A", "ENG", "SIRPA","APOE", "TREM2", "MS4A4A")
Macro_SPP1 <- c("CEBPB","SPP1", "IL1RN", "OLR1", "CXCL2", "VEGFA", "EREG", "C15ORF48", "GPNMB", "PHLDA1", "AQP9", "TNS3", "NDRG1", "FN1", "CLEC5A", "OLR1", "IL1RN")
Macro_Mki67 <- c("MKI67","PCNA","UBE2C")
Macro_CD3E <- c("CD3D","CD3E")
Annother <- c("APOE","APOC1","IL10","HLA-DQA2","HLA-DRA","HLA-DRB1","HLA-DRB5")
genelist <- c(Macro_FCN1,Macro_NLPR3,Macro_PLTP,Macro_IL1B,Macro_C1QC,Macro_SPP1,Macro_Mki67,Macro_CD3E,Annother)
genelist <- genelist[!duplicated(genelist)]
data=FetchData(Macrophage,genelist)
#===============================================================================
metadata=Macrophage@meta.data
colnames(metadata)
metadata=metadata[,c("cell_id","RNA_snn_res.0.5")]
data$cell_id <- row.names(data)
plot_data=merge(metadata,data,by="cell_id")
colnames(plot_data)
plot_data=plot_data[,-1]
plot_data$RNA_snn_res.0.5<- paste("Cluster",plot_data$RNA_snn_res.0.5)
#===============================================================================
plot_data=aggregate(plot_data,by=list(plot_data$RNA_snn_res.0.8),mean)
plot_data=plot_data[,-2]
row.names(plot_data) <- plot_data$Group.1
plot_data=plot_data[,-1]
plot_data=t(plot_data)
#
library(pheatmap)
pheatmap(plot_data,cluster_rows = F,scale = "column",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#===============================================================================
Idents(Macrophage) <- "RNA_snn_res.0.5"
Markers <- FindAllMarkers(Macrophage,min.pct = 0.5,min.diff.pct = 0.25,logfc.threshold = 0.25,only.pos = T)
diff.markers <- Markers%>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
#===============================================================================
Macro_SPP1 <- c(2,3,1,5,6)
Macro_Mki67 <- c(9)
Macro_CD3E <- c(8,14)
Macro_IL1B <- c(13)
Macro_C1QC <- c(4,7)
Macro_PLTP <- c(0,10,12,16)
Macro_FCN1 <- c(15)
Macro_NLPR3 <- c(11)
#+==============================================================================
current.cluster.ids <- c(Macro_Mki67,
                         Macro_FCN1,
                         Macro_SPP1,
                         Macro_NLPR3,
                         Macro_IL1B,Macro_C1QC,Macro_PLTP,Macro_CD3E)

new.cluster.ids <- c(rep("Macro_Mki67",length(Macro_Mki67)),
                     rep("Macro_FCN1",length(Macro_FCN1)),
                     rep("Macro_SPP1",length(Macro_SPP1)),
                     rep("Macro_NLPR3",length(Macro_NLPR3)),
                     rep("Macro_IL1B",length(Macro_IL1B)),
                     rep("Macro_C1QC",length(Macro_C1QC)),
                     rep("Macro_PLTP",length(Macro_PLTP)),
                     rep("Macro_CD3E",length(Macro_CD3E)))

Macrophage@meta.data$LV_03 <- plyr::mapvalues(x = as.character(Macrophage@meta.data$RNA_snn_res.0.8), 
                                              from = current.cluster.ids, to = new.cluster.ids)
Macrophage@meta.data
DimPlot(Macrophage,group.by = "LV_03",reduction = "umap",label = T)
#
tab <- table(Macrophage$LV_03, Macrophage$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#==============================================================================
table(scRNA@meta.data$clMidwayPr)
DC <- subset(scRNA,clMidwayPr=="DC")
DC <- CreateSeuratObject(counts = DC@assays$RNA@counts, 
                         meta.data = DC@meta.data,
                         min.cells = 3, min.features = 200)

DC <- NormalizeData(DC)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
DC<- RunPCA(DC,verbose = F)
ElbowPlot(DC,ndims = 50)
DimHeatmap(DC,cells=500,balanced = T,dims = 1:20)
pc.num=1:10
DC <- DC %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
DC <- FindNeighbors(DC, dims=pc.num)%>% FindClusters(resolution=c(0.1,0.2,0.3,0.4,0.5))
colnames(DC@meta.data)
#===============================================================================
pDC_LILRA4 <- c("cM07","cM08")
cDC2_C1QC <-c("cM05")
cDC2_CD1C <-c("cM04")
cDC1_CLEC9A <- c("cM03")
cDC_LAMP3 <- c("cM09")
cDC_IL22RA2 <- c("cM06")
#===============================================================================
current.cluster.ids <- c(pDC_LILRA4,
                         cDC2_C1QC,
                         cDC2_CD1C,
                         cDC1_CLEC9A,
                         cDC_LAMP3,
                         cDC_IL22RA2)
new.cluster.ids <- c(rep("pDC_LILRA4",length(pDC_LILRA4)),
                     rep("cDC2_C1QC",length(cDC2_C1QC)),
                     rep("cDC2_CD1C",length(cDC2_CD1C)),
                     rep("cDC1_CLEC9A",length(cDC1_CLEC9A)),
                     rep("cDC_LAMP3",length(cDC_LAMP3)),
                     rep("cDC_IL22RA2",length(cDC_IL22RA2)))
DC@meta.data$LV_03 <- plyr::mapvalues(x = as.character(DC@meta.data$cl295v11SubShort), 
                                      from = current.cluster.ids, to = new.cluster.ids)
DimPlot(DC,group.by = "LV_03",reduction = "umap",label = T)
#
library(ggpie)
ggnestedpie(data = Macro, group_key = c("LV_03", "Group"), count_type = "full",
            inner_label_info = "all",
            inner_label_split = NULL,
            inner_label_threshold = 3,# 设置内层环形的阈值
            inner_label_size = 2,
            outer_label_type = "circle", # 设置外层环形
            outer_label_pos = "in",
            outer_label_info = "all",r0 = 1,r2 = 2)
#===============================================================================
#===============================================================================
#Part_05
scRNA <- NormalizeData(TNKIL_cells)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
scRNA<- RunPCA(scRNA,verbose = F)
ElbowPlot(scRNA,ndims = 50)
DimHeatmap(scRNA,cells=500,balanced = T,dims = 1:20)

#
pc.num=1:40
scRNA <- scRNA %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
scRNA@meta.data
DimPlot(scRNA,group.by = "clMidwayPr",reduction = "umap",label = T)
table(scRNA@meta.data$clMidwayPr)
tab <- table(scRNA$clMidwayPr, scRNA$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#===============================================================================
CD8T=subset(scRNA,clMidwayPr=="TCD8")
CD8T <- NormalizeData(CD8T)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
CD8T<- RunPCA(CD8T,verbose = F)
ElbowPlot(CD8T,ndims = 50)
DimHeatmap(CD8T,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:30
CD8T <- CD8T %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
CD8T <- FindNeighbors(CD8T, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
table(CD8T@meta.data$cl295v11SubFull)
DimPlot(CD8T,group.by = "RNA_snn_res.0.8",reduction = "umap",label = T)
VlnPlot(CD8T,features = "TCF7",group.by = "RNA_snn_res.0.6")
#
tab <- table(CD8T$RNA_snn_res.0.8, CD8T$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#==============================================================================-
#===============================================================================
#CD8T细胞注释
CD8_LEF1 <- c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1")
CD8_GPR183 <- c("CCR7","SELL","IL7R","CD27","CD28","PRF1","GZMA","CCL5","GPR183","S1PR1")
CD8_CX3CR1 <- c("KLRG1","CX3CR1","FCGR3A","FGFBP2","PRF1","GZMH","TBX21","EOMES","S1PR1","S1PR5")
CD8_GZMK <- c("GZMK","CXCR4","CXCR3","CD44")
CD8_CD6 <- c("CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE")
CD8_CD160 <- c("CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3","NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE")
CD8_LAYN <- c("HAVCR2","CXCL13","PDCD1","LAYN","TOX",'IFNG',"GZMB","MIR155HG","TNFRSF9","ITGAE")
CD8_SLC4A10 <- c("SLC4A10","KLRB1","ZBTB16","NCR3","RORC","RORA")
#
genelist <- c(CD8_LEF1,
              CD8_GPR183,
              CD8_CX3CR1,
              CD8_GZMK,
              CD8_CD6,
              CD8_CD160,
              CD8_LAYN,
              CD8_SLC4A10)
genelist <- genelist[!duplicated(genelist)]
data=FetchData(CD8T,genelist)
#===============================================================================
metadata=CD8T@meta.data
colnames(metadata)
metadata=metadata[,c("cells","RNA_snn_res.0.8")]
data$cells <- row.names(data)
plot_data=merge(metadata,data,by="cells")
colnames(plot_data)
plot_data=plot_data[,-1]
plot_data$RNA_snn_res.0.8<- paste("Cluster",plot_data$RNA_snn_res.0.8)
#===============================================================================
plot_data=aggregate(plot_data,by=list(plot_data$RNA_snn_res.0.8),mean)
plot_data=plot_data[,-2]
row.names(plot_data) <- plot_data$Group.1
plot_data=plot_data[,-1]
plot_data=t(plot_data)
#
library(pheatmap)
pheatmap(plot_data,cluster_rows = F,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
VlnPlot(CD8T,features ="MKI67",group.by = "RNA_snn_res.0.")
tab <- table(CD8T$RNA_snn_res.0.4, CD8T$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#
DimPlot(CD8T,group.by = "RNA_snn_res.0.8",reduction = "umap",label = T)
DimPlot(CD8T,group.by = "Location",reduction = "umap",label = T)
CD8_LEF1 <- c(9,3)
CD8_GPR183 <-c(4,0) 
CD8_CX3CR1 <- c(14)
CD8_GZMK <- c(12)
CD8_CD6 <- c(1)
CD8_CD160 <- c(6,7,11)
CD8_LAYN <- c(5,8,10,2,13)

#
current.cluster.ids <- c(CD8_LEF1,
                         CD8_GPR183,
                         CD8_CX3CR1,
                         CD8_GZMK,
                         CD8_CD6,
                         CD8_CD160,
                         CD8_LAYN)

new.cluster.ids <- c(rep("CD8_LEF1",length(CD8_LEF1)),
                     rep("CD8_GPR183",length(CD8_GPR183)),
                     rep("CD8_CX3CR1",length(CD8_CX3CR1)),
                     rep("CD8_GZMK",length(CD8_GZMK)),
                     rep("CD8_CD6",length(CD8_CD6)),
                     rep("CD8_CD160",length(CD8_CD160)),
                     rep("CD8_LAYN",length(CD8_LAYN)))

CD8T@meta.data$LV_03 <- plyr::mapvalues(x = as.character(CD8T@meta.data$RNA_snn_res.0.8), 
                                        from = current.cluster.ids, to = new.cluster.ids)
DimPlot(CD8T,group.by = "Group",reduction = "umap",label = T)
DimPlot(CD8T,group.by = "LV_03",reduction = "umap",label = T)
tab <- table(CD8T$LV_03, CD8T$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#===============================================================================
#===============================================================================
TCD4=subset(scRNA,clMidwayPr=="TCD4")
TCD4 <- NormalizeData(TCD4)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
TCD4<- RunPCA(TCD4,verbose = F)
ElbowPlot(TCD4,ndims = 50)
DimHeatmap(TCD4,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:30
TCD4 <- TCD4 %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
TCD4 <- FindNeighbors(TCD4, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
table(TCD4@meta.data$cl295v11SubFull)
DimPlot(TCD4,group.by = "cl295v11SubFull",reduction = "umap",label = T)
VlnPlot(TCD4,features = "CCR7",group.by = "cl295v11SubFull")
tab <- table(TCD4$cl295v11SubFull, TCD4$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#
CD4_GZMK <- "cTNI04 (CD4+ IL7R+CCL5+)"
CD4_CXCL13 <- "cTNI07 (CD4+ CXCL13+)"
CD4_CCR7 <-c("cTNI01 (CD4+ IL7R+)","cTNI02 (CD4+ IL7R+SELL+)","cTNI03 (CD4+ IL7R+HSP+)") 
CD4_IL23R <- "cTNI05 (CD4+ IL17+)"
CD4_CXCR5 <- "cTNI06 (CD4+ TFH)"
CD4_FOXP3 <- "cTNI08 (CD4+ Treg)"
CD4_CTLA4 <- "cTNI09 (CD4+ Treg prolif)"
#
current.cluster.ids <- c(CD4_GZMK,
                         CD4_CXCL13,
                         CD4_CCR7,
                         CD4_IL23R,
                         CD4_CXCR5,
                         CD4_FOXP3,
                         CD4_CTLA4)

new.cluster.ids <- c(rep("CD4_GZMK",length(CD4_GZMK)),
                     rep("CD4_CXCL13",length(CD4_CXCL13)),
                     rep("CD4_CCR7",length(CD4_CCR7)),
                     rep("CD4_IL23R",length(CD4_IL23R)),
                     rep("CD4_CXCR5",length(CD4_CXCR5)),
                     rep("CD4_FOXP3",length(CD4_FOXP3)),
                     rep("CD4_CTLA4",length(CD4_CTLA4)))

TCD4@meta.data$LV_03 <- plyr::mapvalues(x = as.character(TCD4@meta.data$cl295v11SubFull), 
                                        from = current.cluster.ids, to = new.cluster.ids)
#===============================================================================
#===============================================================================
table(scRNA@meta.data$clMidwayPr)
NK=subset(scRNA,clMidwayPr=="NK")
NK <- NormalizeData(NK)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
NK<- RunPCA(NK,verbose = F)
ElbowPlot(NK,ndims = 50)
DimHeatmap(NK,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:15
NK <- NK %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
NK <- FindNeighbors(NK, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
table(NK@meta.data$cl295v11SubFull)
DimPlot(NK,group.by = "cl295v11SubFull",reduction = "umap",label = T)
VlnPlot(NK,features = "CD160",group.by = "cl295v11SubFull")
tab <- table(NK$cl295v11SubFull, NK$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#
table(NK@meta.data$cl295v11SubFull)
NK_CD160 <- "cTNI24 (NK GZMK+)"
NK_CD16 <- c("cTNI23 (NK CD16A+)","cTNI22 (cTNI22)")
NK_XCL1 <- "cTNI25 (NK XCL1+)"
#
current.cluster.ids <- c(NK_CD160,
                         NK_CD16,
                         NK_XCL1)

new.cluster.ids <- c(rep("NK_CD160",length(NK_CD160)),
                     rep("NK_CD16",length(NK_CD16)),
                     rep("NK_XCL1",length(NK_XCL1)))

NK@meta.data$LV_03 <- plyr::mapvalues(x = as.character(NK@meta.data$cl295v11SubFull), 
                                      from = current.cluster.ids, to = new.cluster.ids)
#===============================================================================
Tgd=subset(scRNA,clMidwayPr=="Tgd")
Tgd <- NormalizeData(Tgd)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
Tgd<- RunPCA(Tgd,verbose = F)
ElbowPlot(Tgd,ndims = 50)
DimHeatmap(Tgd,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:15
Tgd <- Tgd %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
Tgd <- FindNeighbors(Tgd, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
table(Tgd@meta.data$cl295v11SubFull)
DimPlot(Tgd,group.by = "cl295v11SubFull",reduction = "umap",label = T)
VlnPlot(Tgd,features = "CD160",group.by = "cl295v11SubFull")
tab <- table(Tgd$cl295v11SubFull, Tgd$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#
gd_like_T <- "cTNI17 (gd-like T)"
gd_like_T_PDCD1 <- "cTNI18 (gd-like T PDCD1+)"
gd_like_T_prolif <- "cTNI19 (gd-like T prolif)"
#
current.cluster.ids <- c(gd_like_T,
                         gd_like_T_PDCD1,
                         gd_like_T_prolif)

new.cluster.ids <- c(rep("gd_like_T",length(gd_like_T)),
                     rep("gd_like_T_PDCD1",length(gd_like_T_PDCD1)),
                     rep("gd_like_T_prolif",length(gd_like_T_prolif)))

Tgd@meta.data$LV_03 <- plyr::mapvalues(x = as.character(Tgd@meta.data$cl295v11SubFull), 
                                       from = current.cluster.ids, to = new.cluster.ids)

#
TZBTB16=subset(scRNA,clMidwayPr=="TZBTB16")
TZBTB16 <- NormalizeData(TZBTB16)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
TZBTB16<- RunPCA(TZBTB16,verbose = F)
ElbowPlot(TZBTB16,ndims = 50)
DimHeatmap(TZBTB16,cells=500,balanced = T,dims = 1:20)
#
pc.num=1:15
TZBTB16 <- TZBTB16 %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
TZBTB16 <- FindNeighbors(TZBTB16, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
table(TZBTB16@meta.data$cl295v11SubFull)
DimPlot(TZBTB16,group.by = "cl295v11SubFull",reduction = "umap",label = T)
VlnPlot(TZBTB16,features = "CD160",group.by = "cl295v11SubFull")
tab <- table(TZBTB16$cl295v11SubFull, TZBTB16$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
#
PLZF_T <- "cTNI20 (PLZF+ T)"
PLZF_T_pro <- "cTNI21 (PLZF+ T prolif)"
#
current.cluster.ids <- c(PLZF_T,
                         PLZF_T_pro)

new.cluster.ids <- c(rep("PLZF_T",length(PLZF_T)),
                     rep("PLZF_T_pro",length(PLZF_T_pro)))

TZBTB16@meta.data$LV_03 <- plyr::mapvalues(x = as.character(TZBTB16@meta.data$cl295v11SubFull), 
                                           from = current.cluster.ids, to = new.cluster.ids)
#=
ILC=subset(scRNA,clMidwayPr=="ILC")
ILC@meta.data$LV_03 <- "ILC"
table(TZBTB16@meta.data$LV_03)
#===============================================================================

T_cell=merge(ILC,TZBTB16)
T_cell=merge(T_cell,Tgd)
T_cell=merge(T_cell,NK)
T_cell=merge(T_cell,TCD4)
T_cell=merge(T_cell,CD8T)

table(T_cell@meta.data$LV_03)
T_cell <- NormalizeData(T_cell)%>% FindVariableFeatures(nfeatures = 3000) %>%ScaleData()
T_cell<- RunPCA(T_cell,verbose = F)
ElbowPlot(T_cell,ndims = 50)
DimHeatmap(T_cell,cells=500,balanced = T,dims = 1:20)

#
pc.num=1:20
T_cell <- T_cell %>%RunTSNE(dims=pc.num)%>% RunUMAP(dims=pc.num)
T_cell <- FindNeighbors(T_cell, dims=pc.num)%>% FindClusters(resolution=c(0.2,0.4,0.6,0.8,1.0,1.2))
T_cell@meta.data
DimPlot(T_cell,group.by = "LV_03",reduction = "umap",label = T)
tab <- table(T_cell$clMidwayPr, T_cell$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
color_plot <- c(ILC,ZBTB16,Tgd,NK,CD4,CD8)
Idents(T_cell) <- "LV_03"
names(color_plot) <- unique(Idents(T_cell))
do_DimPlot(T_cell,legend.position = "none", colors.use = color_plot,reduction = "umap",pt.size=0.8,plot_cell_borders = F)

tab <- table(T_cell$LV_03, T_cell$Group)
chi <- chisq.test(tab)
Roe <- chi$observed / chi$expected
Roe <- as.data.frame(Roe)
#===============================================================================
#===============================================================================
#Part_06
cco.list <- list(low=low,high=high)
cellchat <- mergeCellChat(cco.list, add.names= names(cco.list),cell.prefix = TRUE)
#通讯数量与强度对比
compareInteractions(cellchat,show.legend = F, group= c(1,2), measure = "count")
compareInteractions(cellchat,show.legend = F, group= c(1,2), measure = "weight")
#
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#
#保守和特异性信号通路的识别与可视化
#
rankNet(cellchat, mode ="comparison", stacked = T,do.stat = TRUE)
rankNet(cellchat, mode = "comparison",stacked = F, do.stat = TRUE)
#
pathway.union <-union(cco.list[[1]]@netP$pathways,cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern= "all", signaling =pathway.union, title = names(cco.list)[1],
                                        width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern= "all", signaling= pathway.union,title = names(cco.list)[2],
                                        width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#
pathways.show <- c("WNT")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute= pathways.show)
par(mfrow = c(1,2),xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]],signaling = pathways.show,layout = "circle", 
                      edge.weight.max =weight.max[1],edge.width.max = 10,signaling.name = paste(pathways.show,names(cco.list)[i]))
}
#
netVisual_bubble(cellchat, sources.use = c(1),targets.use = c(4), comparison = c(1, 2), max.dataset = 2,title.name ="Increased signaling in TIL",
                 angle.x = 45, remove.isolate = T)
#
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
# 查看潜在的配体
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)
# 运行分析
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

# 条形图展示配体活性的直方图
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
# top配体靶基因处理
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
# 可视化之前的处理
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
nrow(active_ligand_target_links_df)

# 热图可视化
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Mph-ligands","PS genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
p_ligand_target_network
#===============================================================================
#===============================================================================
Load10X_Spatial_change <- function(data.dir, filename = "filtered_feature_bc_matrix.h5", 
                                   image.name = "tissue_lowres_image.png",
                                   assay = "Spatial", slice = "slice1", filter.matrix = TRUE, 
                                   to.upper = FALSE, image = NULL, ...) 
{
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'", 
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  if(grepl("h5", filename )){
    data <- Read10X_h5(filename = file.path(data.dir, filename), 
                       ...)
  } else {
    #data = Read10X(filtered_dir)
    data = Read10X(paste(data.dir, filename,sep = "/") )
  }
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  if (is.null(x = image)) {
    image <- Read10X_Image(image.dir = file.path(data.dir, 
                                                 "spatial"), 
                           image.name = image.name,
                           filter.matrix = filter.matrix)
  }
  else {
    if (!inherits(x = image, what = "VisiumV1")) 
      stop("Image must be an object of class 'VisiumV1'.")
  }
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  if(image.name == "tissue_lowres_image.png") {
    object = object
  }else {
    object@images[[1]]@scale.factors$lowres = object@images[[1]]@scale.factors$hires
    #object@images$slice1@scale.factors$lowres = object@images$slice1@scale.factors$hires
    
  }
  return(object)
}
#
library(Seurat)
data.dir = ""

T1 <- Load10X_Spatial_change(data.dir = data.dir,
                             filename = "filtered_feature_bc_matrix")
for (i in colnames((T1@images$slice1@coordinates))) {
  T1@images$slice1@coordinates[[i]] <- as.integer(T1@images$slice1@coordinates[[i]])
}
p1 <- SpatialDimPlot(T1,alpha = 0)
p2 <- SpatialFeaturePlot(T1, features = "nFeature_Spatial",
                         pt.size.factor = 3)
p1 + p2
scRNA=T1
scRNA <- SCTransform(scRNA, assay = "Spatial", verbose = FALSE)
scRNA <- RunPCA(scRNA, assay = "SCT", verbose = FALSE,npcs = 50)
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:50)
scRNA <- FindClusters(scRNA, verbose = FALSE,,resolution=c(1.0,1.2))
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:50)
DimPlot(scRNA, reduction = "umap", label = TRUE)
scRNA@meta.data
SpatialDimPlot(scRNA, label = TRUE, label.size = 5,group.by = "SCT_snn_res.1.2")
#重新定义细胞类群
Epi <- c(13,11,0,10,12,9,7,14,8)
Immune <-c(3,1,4) 
Stromal <- c(6,5,2)
#
current.cluster.ids <- c(Epi,
                         Immune,
                         Stromal)

new.cluster.ids <- c(rep("Epi",length(Epi)),
                     rep("Immune",length(Immune)),
                     rep("Stromal",length(Stromal)))

scRNA@meta.data$LV_01 <- plyr::mapvalues(x = as.character(scRNA@meta.data$SCT_snn_res.1.2), 
                                         from = current.cluster.ids, to = new.cluster.ids)
SpatialDimPlot(scRNA, label = TRUE, label.size = 5,group.by = "LV_01",pt.size.factor = 1)
#基因集评分
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(msigdbr)
homo_KEGG = msigdbr(species = "Homo sapiens",
                    category = "H") %>% dplyr::select(gs_name,gene_symbol)
#
homo_KEGG_gene = homo_KEGG %>% split(x =.$gene_symbol, f =.$gs_name)

#选择其中一条通路（我这里选择的是氨基酸和核苷酸糖代谢），将其也转为list
features <- list(homo_KEGG_gene$HALLMARK_G2M_CHECKPOINT)
scRNA<- AddModuleScore(scRNA,
                       features = features,
                       ctrl = 100,
                       name = "HALLMARK_G2M_CHECKPOINT")
scRNA@meta.data
SpatialFeaturePlot(object = scRNA, features = "", alpha = c(0.1, 1),pt.size.factor = 1)
#===============================================================================
#===============================================================================