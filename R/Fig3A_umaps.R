# from the run_gsva.R
library(SingleCellExperiment)
library(scuttle)
library(scater)
library(ggpubr)
library(cowplot)
library(scran)
library(clusterProfiler)
library(org.Hs.eg.db)
context = 'COAD'
countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancersea/scRNASeq/COAD/GSE81861_CRC_tumor_all_cells_FPKM.csv"))
dfmat = countdf[,1:ncol(countdf)]
rownames(dfmat) = sapply(countdf$X, function(x) paste0('ENSG', strsplit(x, "ENSG")[[1]][2]))
rownames(dfmat) = sapply(rownames(dfmat), function(x) strsplit(x, "[.]")[[1]][1])
dfmat = dfmat[rownames(dfmat) %in% annotationdf$ENSG,]
subdfmat = dfmat[,sapply(colnames(dfmat), function(x) grepl('Epithelial|Fibroblast|Endothelial', x))]
alldfmat = dfmat[,sapply(colnames(dfmat), function(x) grepl('Bcell|Macrophage|MastCell|Tcell|Epithelial|Fibroblast|Endothelial', x))]
nondfmat = dfmat[,sapply(colnames(dfmat), function(x) grepl('Bcell|Macrophage|MastCell|Tcell', x))]

current_mat = subdfmat
# current_mat = nondfmat
# original
gsva.es = gsva.es_list[[context]]

# context specific
gsva.es_context = gsva.es_context_list[[context]]

# discarded
gsva.es_discard = gsva.es_discard_list[[context]]
gsva.es_retain = gsva.es_retain_list[[context]]
gsva.es_hpa = gsva.es_hpa_list[[context]]


set.seed(1234)
sce = SingleCellExperiment(list(counts=current_mat))
sce = logNormCounts(sce)
sce = runUMAP(sce)
sce$UMAP1 = reducedDims(sce)$UMAP[,1]
sce$UMAP2 = reducedDims(sce)$UMAP[,2]
sce$Context_agnostic = gsva.es['Quiescence',]
sce$COAD_specific = gsva.es_context['Quiescence',]
sce$discarded = gsva.es_discard['Quiescence',]
sce$retained = gsva.es_retain['Quiescence',]
sce$hpa = gsva.es_hpa['Quiescence',]
sce$scaled_rank1 = scales::rescale(rank(sce$Context_agnostic), range(sce$Context_agnostic))
sce$scaled_rank2 = scales::rescale(rank(sce$COAD_specific), range(sce$COAD_specific))
sce$scaled_rank3 = scales::rescale(rank(sce$discarded), range(sce$discarded))
sce$scaled_rank4 = scales::rescale(rank(sce$retained), range(sce$retained))
sce$scaled_rank5 = scales::rescale(rank(sce$hpa), range(sce$hpa))

cpallete = 'Spectral'
b = c(-0.5,-0.2, 0,0.2, 0.5)
limits = c(min(sce$scaled_rank1, sce$scaled_rank2, sce$scaled_rank3, sce$scaled_rank4, sce$scaled_rank5), 
           max(sce$scaled_rank1, sce$scaled_rank2, sce$scaled_rank3, sce$scaled_rank4, sce$scaled_rank5))
p1 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank1'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='Context-agnostic')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p1

p2 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank2'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='COAD-specific')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p2

p3 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank3'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='COAD-discarded')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p3

p4 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank4'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='COAD-retained')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p4

p5 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank5'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='COAD-FrostRefined')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p5

cowplot::plot_grid(p3, p5, labels = c('', ''))
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3a_coad_umap_sup.png'), 
        width = 10, height = 6, dpi=700) 

#markers for each group
cell_to_num = list(Epithelial=1, Fibroblast=2, Endothelial=3)
sce$celltype = unlist(sapply(rownames(colData(sce)), function(x) 
  strsplit(x,'__')[[1]][2]))

sce$cluster = as.factor(as.numeric(sce$UMAP1>0.0 & sce$UMAP2>0.0))
p3 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='cluster'))+
  geom_point(size=1)+scale_colour_manual(values=brewer.pal(6, "Dark2"))+
  labs(col='COAD-specific')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p3

plotColData(sce, x = "UMAP1", y="UMAP2", colour_by="cluster") 

out <- findMarkers(sce, groups=sce$cluster)
high = rownames(out[['1']][out[['1']]$logFC.0>0.2& out[['1']]$FDR<0.1,])
low = rownames(out[['0']][out[['0']]$logFC.1>2& out[['0']]$FDR<0.001,])
yy <- enrichGO(high, 
               'org.Hs.eg.db',keyType = "ENSEMBL", 
               ont="BP", pvalueCutoff=0.05); dim(yy)
head(yy$Description)

yy <- enrichGO(low, 
               'org.Hs.eg.db',keyType = "ENSEMBL", 
               ont="BP", pvalueCutoff=0.05); dim(yy)
head(yy$Description)

ck <- compareCluster(geneClusters = list(Active=high, Inactive=low), OrgDb = org.Hs.eg.db,
                     fun = 'enrichGO', pvalueCutoff=0.05, keyType="ENSEMBL")
head(ck) 
dotplot(ck)

b = c(0,2.5,5,7.5,10)
limits = c(min(umapdf$EPCAM, umapdf$CEACAM5), max(umapdf$EPCAM, umapdf$CEACAM5))

p3 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='ENSG00000119888'))+
  geom_point(size=1)+scale_color_distiller(palette='Spectral')+
  labs(col='EPCAM')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p3

p4 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='ENSG00000105388'))+
  geom_point(size=1)+scale_color_distiller(palette='Spectral')+
  labs(col='CEACAM5')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p4

# legacy
umapdf = as.data.frame(reducedDims(sce)$UMAP); head(umapdf)
umapdf$`Context_agnostic` = gsva.es['Quiescence',]
umapdf$`COAD_specific` = gsva.es_context['Quiescence',]
umapdf$celltype = 'cancer'
if(context=='COAD'){
  umapdf$celltype[sapply(rownames(umapdf), function(x) grepl('Bcell|Macrophage|MastCell|Tcell', x))] = 'non-cancer'
}
logdf = assays(sce)$logcounts
umapdf$EPCAM = as.numeric(logdf['ENSG00000119888',]); summary(umapdf$EPCAM)
umapdf$CEACAM5 = as.numeric(logdf['ENSG00000105388',]); summary(umapdf$CEACAM5)


p3 = umapdf %>% ggplot(aes(x=V1, y=V2, color=EPCAM))+
  scale_color_gradient2(low='blue', high='red', breaks=b, labels=format(b), limits=limits, midpoint=5)+ # 
  geom_point()+labs(x='UMAP1', y='UMAP2')+theme_bw()+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
  theme(legend.position = "top"); p3

p4 = umapdf %>% ggplot(aes(x=V1, y=V2, color=CEACAM5))+
  scale_color_gradient2(low='blue', high='red', breaks=b, labels=format(b), limits=limits, midpoint=5)+ # 
  geom_point()+labs(x='UMAP1', y='UMAP2')+theme_bw()+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
  theme(legend.position = "top"); p4

cowplot::plot_grid(p3, p4, labels = c('A', 'B'))

umapdf$active = umapdf$COAD_specific>0.2
umapdf %>% ggplot(aes(x=active, y=EPCAM, color=EPCAM))+
  geom_boxplot()+labs(x='COAD_specific', y='CEACAM5')+theme_bw()+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
  theme(legend.position = "top"); 



 #---------------- gbm

context = 'GBM'
countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/GBM/GSE84465_GBM_All_data.csv"), sep=" ")
mat = countdf
mat = mat[rownames(mat) %in% annotationdf$gene_name,]
rownames(mat) = sapply(rownames(mat), function(x) annotationdf$ENSG[annotationdf$gene_name==x])



# current_mat = alldfmat
current_mat = as.matrix(mat)
reshape_mat <- sapply(current_mat, as.numeric)
reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                     dimnames=list(rownames(current_mat), colnames(current_mat)))

# current_mat = nondfmat
# original
gs = lapply(signatures, read_genes)
names(gs) = signatures
gsvapar = gsvaParam(reshape_mat, gs)
gsva.es = gsva(gsvapar)
#gsva.es = gsva(as.matrix(current_mat), gs, verbose=FALSE)

# context specific
gs_context = lapply(signatures, function(x) read_genes(x, context))
names(gs_context) = signatures
gsvapar = gsvaParam(reshape_mat, gs_context)
gsva.es_context = gsva(gsvapar)
#gsva.es_context = gsva(as.matrix(current_mat), gs_context, verbose=FALSE)


sce = SingleCellExperiment(list(counts=current_mat))
sce = logNormCounts(sce)
sce = runUMAP(sce)
sce$UMAP1 = reducedDims(sce)$UMAP[,1]
sce$UMAP2 = reducedDims(sce)$UMAP[,2]
sce$Context_agnostic = gsva.es['Angiogenesis',]
sce$GBM_specific = gsva.es_context['Angiogenesis',]
sce$discarded = gsva.es_discard['Angiogenesis',]
sce$scaled_rank1 = scales::rescale(rank(sce$Context_agnostic), range(sce$Context_agnostic))
sce$scaled_rank2 = scales::rescale(rank(sce$GBM_specific), range(sce$GBM_specific))
sce$scaled_rank3 = scales::rescale(rank(sce$discarded), range(sce$discarded))

#cpallete = 'Spectral'
b = c(-0.5,-0.2, 0,0.2, 0.5)
limits = c(min(sce$scaled_rank1, sce$scaled_rank2), max(sce$scaled_rank1, sce$scaled_rank2))

p1 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank1'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral', breaks=b, labels=format(b), limits=limits)+
  labs(col='Context-agnostic')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p1

p2 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank2'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral', breaks=b, labels=format(b), limits=limits)+
  labs(col='GBM-specific')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p2

p3 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank3'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral', breaks=b, labels=format(b), limits=limits)+
  labs(col='GBM-discarded')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p3

cowplot::plot_grid(p1, p3, labels = c('', ''))
ggsave2('/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_3a_gbm_umap.png', 
        width = 10, height = 6,  dpi=700)

p3 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='ENSG00000106991'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral')+
  labs(col='Endoglin')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p3

p4 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='ENSG00000011422'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral')+
  labs(col='PLAUR')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p4

cowplot::plot_grid(p3, p4, labels = c('', ''))
ggsave2('/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_3a_gbm_marker_umap.png', 
        width = 10, height = 6, dpi=700)


cell_to_matkers = list(Myeloid=c('CD45', 'HLA'), Neoplastic=c('EGFR', 'SOX9'), 
                       Oligodendrocyte=c('MOG', 'MBP', 'OPALIN'), Vascular=c('DCN'),
     OPC=c('GPR17'), Endothelial=c('HIF1A'), Neurons=c('STMN2', 'L1CAM'), 
     Astrocytes=c('AGXT2L1', 'ALDH1L1', 'WIF1', 'NTSR2'))
#'VEGFA', 
symbol_to_ensg = list(CD45='ENSG00000081237', HLA='ENSG00000196126',
                      EGFR='ENSG00000146648', SOX9='ENSG00000125398',
                      MOG='ENSG00000204655', MBP='ENSG00000197971', OPALIN='ENSG00000197430',
                      DCN='ENSG00000011465', GPR17='ENSG00000144230', HIF1A='ENSG00000100644',
                      STMN2='ENSG00000104435', L1CAM='ENSG00000198910',
                      ALDH1L1='ENSG00000144908', 
                      WIF1='ENSG00000156076', NTSR2='ENSG00000169006')
sapply(names(symbol_to_ensg), function(gene){
  print(gene)
  p5 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col=symbol_to_ensg[[gene]]))+
    geom_point(size=0.4)+scale_color_distiller(palette='Spectral')+
    labs(col=gene)+theme_pubr(base_size = 18)+ labs_pubr()+
    theme(legend.position = "top"); print(p5)
})

p5 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='ENSG00000081237'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral')+
  labs(col='CD45')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p5

p6 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='ENSG00000196126'))+
  geom_point(size=0.4)+scale_color_distiller(palette='Spectral')+
  labs(col='HLA-DRB1')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p6

cowplot::plot_grid(p5, p6, labels = c('', ''))
ggsave2('/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_3a_gbm_myeloid_umap.png', 
        width = 10, height = 6, dpi=700)


# box/violin for - Endoglin, Plaur, cd45, hla-drb1, boxplots for cells > 0 of umap1
genes = c('ENSG00000106991', 'ENSG00000011422', 'ENSG00000081237', 'ENSG00000196126')
expdf = melt(logcounts(sce[genes,]))
expdf$cluster = 'negative'
expdf$cluster[expdf$Var2 %in% colnames(sce)[sce$UMAP1>0]] = 'positive'
symbols <- c('ENSG00000106991'="Endoglin", 'ENSG00000011422'="PLAUR", 'ENSG00000081237'="CD45", 'ENSG00000196126'="HLA-DRB1")
expdf$symbol <- sapply(expdf$Var1, function(x) symbols[x])

ggplot(expdf, aes(x = factor(symbol, level=c('Endoglin', 'PLAUR', 'CD45', 'HLA-DRB1')), y = value)) + 
  geom_boxplot(aes(fill = cluster), position = position_dodge(0.9)) + theme_bw()+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+xlab('Gene')+ylab('Logcounts')+
  theme(legend.position = "top"); 
ggsave2('/Users/sharminm2/OneDrive - National Institutes of Health/context/plots/hPPIN_DE2/Fig_3b_gbm_boxplots.png', 
        width = 10, height = 6, dpi=700)

#cluster and find biomarkers
library(bluster)
sce$cluster <- clusterCells(sce, use.dimred="UMAP", BLUSPARAM=NNGraphParam(cluster.fun="louvain"))

p5 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='cluster'))+
  geom_point(size=0.4)+#scale_color_distiller(palette='Spectral')+
  labs(col='cluster')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p5

out <- findMarkers(sce, groups=sce$cluster, pval.type="all")



# legacy
b = c(0,2.5,5,7.5,10)
limits = c(min(umapdf$EGFR, umapdf$SOX9), max(umapdf$EGFR, umapdf$SOX9))

umapdf = as.data.frame(reducedDims(sce)$UMAP); head(umapdf)
umapdf$Context_agnostic = gsva.es['Angiogenesis',]
umapdf$GBM_specific = gsva.es_context['Angiogenesis',]
umapdf$celltype = 'cancer'
logdf = assays(sce)$logcounts
umapdf$EGFR = as.numeric(logdf['ENSG00000146648',]); summary(umapdf$EGFR)
umapdf$SOX9 = as.numeric(logdf['ENSG00000125398',]); summary(umapdf$SOX9)
umapdf$Endoglin = as.numeric(logdf['ENSG00000106991',]); summary(umapdf$Endoglin)
umapdf$PLAUR = as.numeric(logdf['ENSG00000011422',]); summary(umapdf$PLAUR)

p5 = umapdf %>% ggplot(aes(x=V1, y=V2, color=EGFR))+
  scale_color_distiller(palette='Spectral', breaks=b, labels=format(b), limits=limits)+ # 
  geom_point(size=0.4)+labs(x='UMAP1', y='UMAP2')+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
  theme(legend.position = "top"); p5

p6 = umapdf %>% ggplot(aes(x=V1, y=V2, color=SOX9))+
  scale_color_distiller(palette='Spectral', breaks=b, labels=format(b), limits=limits)+ # 
  geom_point(size=0.4)+labs(x='UMAP1', y='UMAP2')+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
  theme(legend.position = "top"); p6

umapdf$active = ifelse(umapdf$GBM_specific>0.1, 'Active', 'Inactive')
umapdf %>% ggplot(aes(x=active, y=FMOD, color=FMOD))+
  geom_boxplot(size=0.4)+labs(x='GBM_specific', y='FMOD')+
  theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
  theme(legend.position = "top")


## BRCA
context = 'BRCA'
cdat = read.csv(glue("/Users/sharminm2/Dropbox/nih/context/data/cancerSEA/scRNASeq/Breast/GSE75688_final_sample_information.txt"), sep="\t")
countdf = read.csv(glue("/Users/sharminm2/Dropbox/nih/context/data/cancerSEA/scRNASeq/Breast/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"), sep="\t")
#dfmat = countdf[,cdat[cdat$type=='SC' & cdat$index3=='Tumor',]$sample]
dfmat = countdf[,cdat[cdat$type=='SC',]$sample] # all single cells
#dfmat = countdf[,4:ncol(countdf)]
#rownames(dfmat) = rownames(countdf)
rownames(dfmat) = sapply(countdf$gene_id, function(x) strsplit(x, "[.]")[[1]][1])
dfmat = dfmat[countdf$gene_type=='protein_coding',]

current_mat = dfmat
# current_mat = nondfmat
# original
gs = lapply(signatures, read_genes)
names(gs) = signatures
gsva.es = gsva(as.matrix(current_mat), gs, verbose=FALSE)

# context specific
gs_context = lapply(signatures, function(x) read_genes(x, context))
names(gs_context) = signatures
gsva.es_context = gsva(as.matrix(current_mat), gs_context, verbose=FALSE)

set.seed(1234)
sce = SingleCellExperiment(list(counts=current_mat))
sce = logNormCounts(sce)
sce = runUMAP(sce)
sce$UMAP1 = reducedDims(sce)$UMAP[,1]
sce$UMAP2 = reducedDims(sce)$UMAP[,2]
sce$Context_agnostic = gsva.es['Quiescence',]
sce$BRCA_specific = gsva.es_context['Quiescence',]

sce$scaled_rank1 = scales::rescale(rank(sce$Context_agnostic), range(sce$Context_agnostic))
sce$scaled_rank2 = scales::rescale(rank(sce$BRCA_specific), range(sce$BRCA_specific))


cpallete = 'Spectral'
b = c(-0.5,-0.2, 0,0.2, 0.5)
limits = c(min(sce$scaled_rank1, sce$scaled_rank2), max(sce$scaled_rank1, sce$scaled_rank2))
p1 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank1'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='Context-agnostic')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p1

p2 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='scaled_rank2'))+
  geom_point(size=1)+scale_color_distiller(palette=cpallete, breaks=b, labels=format(b), limits=limits)+
  labs(col='BRCA-specific')+theme_pubr(base_size = 25)+ labs_pubr()+
  theme(legend.position = "top"); p2
cowplot::plot_grid(p1, p2, labels = c('', ''))
ggsave2('/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_3c_brca_umap.png', 
        width = 10, height = 6, dpi=700) 

sce$subtype = sapply(rownames(colData(sce)), function(x) strsplit(x, '_')[[1]][1])
sce$index = sapply(rownames(colData(sce)), function(x) cdat$index[cdat$sample==x])
sce$index2 = sapply(rownames(colData(sce)), function(x) cdat$index2[cdat$sample==x])
sce$index3 = sapply(rownames(colData(sce)), function(x) cdat$index3[cdat$sample==x])

p3 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='index2'))+
  geom_point(size=0.4)+scale_colour_manual(values=brewer.pal(6, "Dark2"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(col='cluster')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p3

p4 = scater::ggcells(sce, aes_string(x='UMAP1', y='UMAP2', col='subtype'))+
  geom_point(size=0.4)+scale_colour_discrete(drop=TRUE,limits = levels(sce$subtype))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(col='cluster')+theme_pubr(base_size = 18)+ labs_pubr()+
  theme(legend.position = "top"); p4

cowplot::plot_grid(p3, p4, labels = c('', ''))
ggsave2('/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_3c_brca_umap_cells.png', 
        width = 10, height = 6, dpi=700) 

