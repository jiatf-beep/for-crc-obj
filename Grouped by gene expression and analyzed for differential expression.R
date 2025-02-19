library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)
vdjcd4 <- readRDS( 'F:/rds-liu/8-26/vdjcd4.rds')
CXCL13 <- vdjcd4@assays$RNA
gene_CXCL13 <- CXCL13['CXCL13',]
gene_CXCL13 <- as.data.frame(gene_CXCL13)
gene_CXCL13 <-as.data.frame(t(gene_CXCL13))
summary(gene_CXCL13$CXCL13)
gene_CXCL13$group1 <- ifelse(gene_CXCL13$CXCL13 > 0,'CXCL13+','CXCL13-')
unique(gene_CXCL13$group1)
vdjcd4@meta.data$exCXCL13 <- gene_CXCL13$group1
vdjcd4$exCXCL13 <- factor(vdjcd4$exCXCL13,levels = c("CXCL13+", "CXCL13-"))
TNFRSF18 <- vdjcd4@assays$RNA
gene_TNFRSF18 <- TNFRSF18['TNFRSF18',]
gene_TNFRSF18 <- as.data.frame(gene_TNFRSF18)
gene_TNFRSF18 <-as.data.frame(t(gene_TNFRSF18))
summary(gene_TNFRSF18$TNFRSF18)
gene_TNFRSF18$group1 <- ifelse(gene_TNFRSF18$TNFRSF18 > 0,'TNFRSF18+','TNFRSF18-')
unique(gene_TNFRSF18$group1)
vdjcd4@meta.data$exTNFRSF18 <- gene_TNFRSF18$group1
vdjcd4$exTNFRSF18 <- factor(vdjcd4$exTNFRSF18,levels = c("TNFRSF18+", "TNFRSF18-"))
unique(vdjcd4$exCXCL13)
dir.create('F:/rds-liu/2025.1.17/cd4t_cxcl13_deg')
deg <-FindMarkers(vdjcd4, ident.1 = 'CXCL13+',ident.2 = 'CXCL13-',
                  group.by = vdjcd4$exCXCL13,logfc.threshold =0, min.pct = 0.1,pseudocount.use = 0.01)
# use the function -- add_regulate to add a regulate column 
# to the DEG result Pdeg. 
library('ggVolcano')
deg <- add_regulate(deg, log2FC_name = "avg_log2FC",
                    fdr_name = "p_val",log2FC = 1, fdr = 0.05)
deg$geneName <- rownames(deg)
p1 <- ggvolcano(deg, x = 'log2FoldChange', y = 'padj',
                label = "row", label_number = 10,  log2FC_cut = 1,output = FALSE,FDR_cut = 0.05,
                fills = c("blue","#b4b4d8","#e94234"),
                colors = c("blue","#b4b4d8","#e94234"))

p1

ggsave('F:/rds-liu/2025.1.17/cd4t_cxcl13_deg/CXCL13+vsCXCL13-.png', p1, width = 5, height = 5)

##GSEA
Markers_genelist <- deg$log2FoldChange
names(Markers_genelist)= rownames(deg)
head(Markers_genelist)
Markers_genelist <- sort(Markers_genelist, decreasing = T)
# 导入MSigDB
m_df = msigdbr(species = 'Homo sapiens' , category = "C2")
mf_df = m_df %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df)<-c("term","gene")
gsea.results <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea.results)
gseaplot(gsea.results, gsea.results@result$ID[1])
library('GseaVis')
dotplotGsea(data = gsea.results,topn = 10,
            order.by = 'NES',
            add.seg = T)
write.csv(gsea.results@result,file = 'F:/rds-liu/2025.1.17/cd4t_cxcl13_deg/vdjcd4CXCL13_gsea.csv')

dir.create('F:/rds-liu/2025.1.17/cd4t_TNFRSF18_deg')
deg <-FindMarkers(vdjcd4, ident.1 = 'TNFRSF18+',ident.2 = 'TNFRSF18-',
                  group.by = vdjcd4$exTNFRSF18,logfc.threshold =0, min.pct = 0.1,pseudocount.use = 0.01)
# use the function -- add_regulate to add a regulate column 
# to the DEG result Pdeg. 
library('ggVolcano')
deg <- add_regulate(deg, log2FC_name = "avg_log2FC",
                    fdr_name = "p_val",log2FC = 1, fdr = 0.05)
deg$geneName <- rownames(deg)
p1 <- ggvolcano(deg, x = 'log2FoldChange', y = 'padj',
                label = "row", label_number = 10,  log2FC_cut = 1,output = FALSE,FDR_cut = 0.05,
                fills = c("blue","#b4b4d8","#e94234"),
                colors = c("blue","#b4b4d8","#e94234"))

p1

ggsave('F:/rds-liu/2025.1.17/cd4t_TNFRSF18_deg/TNFRSF18+vsTNFRSF18-.png', p1, width = 5, height = 5)

##GSEA富集分析
Markers_genelist <- deg$log2FoldChange
names(Markers_genelist)= rownames(deg)
head(Markers_genelist)
Markers_genelist <- sort(Markers_genelist, decreasing = T)
# 导入MSigDB
m_df = msigdbr(species = 'Homo sapiens' , category = "C2")
mf_df = m_df %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df)<-c("term","gene")
gsea.results <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea.results)
gseaplot(gsea.results, gsea.results@result$ID[1])
library('GseaVis')
dotplotGsea(data = gsea.results,topn = 10,
            order.by = 'NES',
            add.seg = T)
write.csv(gsea.results@result,file = 'F:/rds-liu/2025.1.17/cd4t_TNFRSF18_deg/vdjcd4TNFRSF18_gsea.csv')
saveRDS(vdjcd4,file = 'F:/rds-liu/8-26/vdjcd4.rds')
