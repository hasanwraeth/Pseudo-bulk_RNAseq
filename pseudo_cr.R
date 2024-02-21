#Library---------------
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(metap)
library(presto)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(velocyto.R)
library(pagoda2)
library(SeuratWrappers)
library(monocle3)
library(CellChat)
library(viridis)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

#Pseudobulk RNAseq-------------------

counts=H@assays$RNA@counts
metadata=H@meta.data
metadata$cluster_id=factor(H@active.ident)
sce=SingleCellExperiment(assays=list(counts=counts),
                         colData=metadata)

saveRDS(sce, file = "./HLA_sce.rds")
assays(sce)
dim(counts(sce))
counts(sce)[1:6, 1:6]

dim(colData(sce))
head(colData(sce))

cluster_names=levels(colData(sce)$cluster_id)
cluster_names
length(cluster_names)

groups <- colData(sce)[, c("cluster_id", "seurat_clusters")]
head(groups)

aggr_counts=aggregate.Matrix(t(counts(sce)),
                             groupings=groups, fun='sum')

class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

aggr_counts=t(aggr_counts)
aggr_counts[1:6, 1:6]

tstrsplit(colnames(aggr_counts),"_") %>% str()
head(colnames(aggr_counts),n=10)
head(tstrsplit(colnames(aggr_counts),"_")[[1]],n=10)

chol_idx=which(head(tstrsplit(colnames(aggr_counts),"_")[[1]]=="Cholangiocyte"))
chol_idx
colnames(aggr_counts)[chol_idx]
aggr_counts[1:10, chol_idx]

cluster_names
counts_ls=list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

str(counts_ls)
counts_ls=counts_ls[-2]

head(colData(sce))
metadata=colData(sce) %>% as.data.frame() %>%
  dplyr::select(seurat_clusters, cluster_id)
dim(metadata)
head(metadata)
metadata=metadata[!duplicated(metadata),]
dim(metadata)
head(metadata)

rownames(metadata)=metadata$seurat_clusters
head(metadata)

t=table(colData(sce)$seurat_clusters,
        colData(sce)$cluster_id)
t[1:6,1:4]

metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$seurat_clusters  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$seurat_cluster, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  ###Not usedful in this case
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  #df <- plyr::join(df, metadata, 
  #                 by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
}

# Explore the different components of the list
str(metadata_ls)

cluster_names
all(names(counts_ls)==names(metadata_ls))
idx=which(names(counts_ls)=="Endothelium")
cluster_counts=counts_ls[[idx]]
cluster_counts=cbind(cluster_counts,counts_ls[[idx]])
cluster_metadata=metadata_ls[[idx]]
cluster_metadata=rbind(cluster_metadata,metadata_ls[[idx]])
cluster_counts[1:6,1:2]
head(cluster_metadata)
all(colnames(cluster_counts)==rownames(cluster_metadata))

dds=DESeqDataSetFromMatrix(cluster_counts,
                           colData=cluster_metadata,
                           design=~cluster_id)
rld=rlog(dds,blind=T)
DESeq2::plotPCA(rld,ntop=500,intgroup="cluster_id")
rld_mat=assay(rld)
rld_cor=cor(rld_mat)
pheatmap(rld_cor,annotation = cluster_metadata[,c("cluster_id"),drop=F])


dds=DESeq(dds)
plotDispEsts(dds)
resultsNames(dds)
res=results(dds, 
            name="cluster_id_Endothelium_vs_Cholangiocyte",
            alpha=0.05)

res=lfcShrink(dds, res=res, type="apeglm",
            coef="cluster_id_Endothelium_vs_Cholangiocyte")

res_tbl=res %>% data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% arrange(padj)
res_tbl

write.csv(res_tbl, "Endo_vs_Chol.csv", row.names = F)

sig_res=dplyr::filter(res_tbl, padj < 0.05) %>%
  dplyr::arrange(padj)
sig_res1=sig_res[abs(sig_res$log2FoldChange)>1,]
norm_counts=counts(dds, normalized=T)
sig_counts=norm_counts[rownames(norm_counts) %in% 
                         sig_res1$gene,]

## Set a color-blind friendly palette
heat_colors <- rev(brewer.pal(11, "PuOr"))

## Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_m, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 

q_mine=c('KRT19','PLAT','EPCAM','SEMA5A','TACSTD2',
         'AQP1','MMP7','SOX9','ATP1A1','FGFR2',
         'TUBB1','KCNC1','PRSS48','GPRC5D','GPR34',
         'CRAMP1','DNASE1','CRHR2','AKAP10','AKAP13')
sig_m=subset(sig_counts, rownames(sig_counts) %in% q_mine)






