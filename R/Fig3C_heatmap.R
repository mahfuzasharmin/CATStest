library(GSVA)
library(glue)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(org.Hs.eg.db)
library(cowplot)
library(ggpubr)
library(SingleCellExperiment)
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  #library(rafalib)
  library(umap)
})
#library(dplyr)
#library(tidyverse)

signatures = c("Angiogenesis", "Apoptosis", "Cell_Cycle", "Differentiation",
               "DNA_damage", "DNA_repair", "EMT", "Hypoxia", "Inflammation",
               "Invasion", "Metastasis", "Proliferation", "Quiescence", "Stemness") #"SenMayo", "Neuroendocrine"
cancerTypes = c('BRCA', 'COAD', 'ESCA', 'GBM', 'PAAD', 'PRAD') #'LUAD', 
cancerType_to_tissue = list('BRCA'='breast', 'COAD'='colon', 'ESCA'='Esophagus', 
                            'GBM'='brain', #'LUAD'='lung', 
                            'PAAD'='pancreas', 'PRAD'='Prostate')

cancerType_to_hpa_tissue = list('ACC'='adrenal.gland', 'BRCA'='breast', 'CESC'='cervix', 
                                'COAD'='colon', 'ESCA'='esophagus', 'KICH'='kidney', 
                                'KIRC'='kidney', 'KIRP'='kidney',  'LIHC'='liver', 'LUAD'='lung', 
                                'LUSC'='lung', 'OV'='ovary', 'PAAD'='pancreas', 'PRAD'='prostate', 
                                'READ'='rectum', 'SKCM'='skin' , 'STAD'='stomach', 'THCA'='thyroid.gland')

USERNAME = "sharminm2"
DB_location = glue("/Users/{USERNAME}/OneDrive - National Institutes of Health/db")
PROJECT_location = glue("/Users/{USERNAME}/OneDrive - National Institutes of Health/context")
GENE_NETWORK = 'hPPIN'
DE_K = '2'
CUR_DST = glue("{GENE_NETWORK}_DE{DE_K}")
extra = ''
if(GENE_NETWORK=='STRING'){
  extra = 'String_'
}

get_annotation = function(){
  genedf = read.csv(glue("{DB_location}/gencode/gencode.v26.GRCh38.genes.gtf"),
                    header=F, sep="\t", comment.char="#")
  head(genedf)
  pcgenedf = genedf[genedf$V3=='gene' & grepl("protein_coding", genedf$V9),]; dim(pcgenedf) #
  annotationdf = data.frame(Reduce(rbind, strsplit(pcgenedf$V9, ";")))[c('X1', 'X4')]
  rownames(annotationdf) = 1:nrow(annotationdf)
  colnames(annotationdf) = c('gene_id', 'gene_name')
  annotationdf$gene_id = sapply(annotationdf$gene_id, function(x) strsplit(x, " ")[[1]][2])
  annotationdf$gene_name = sapply(annotationdf$gene_name, function(x) strsplit(x, " ")[[1]][3])
  annotationdf = cbind(pcgenedf[c('V1', 'V4', 'V5', 'V7')], annotationdf)
  colnames(annotationdf) = c('chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name')
  rownames(annotationdf) = 1:nrow(annotationdf)
  annotationdf$ENSG = unlist(lapply(strsplit(annotationdf$gene_id, "[.]"), function(x) x[[1]][1]))
  annotationdf = annotationdf[!duplicated(annotationdf$gene_name),]
  rownames(annotationdf) = annotationdf$ENSG
  return(annotationdf)
}

annotationdf = get_annotation()

read_genes <- function(signature, cancerType=NULL, setname=NULL){
  if (is.null(cancerType)){
    genefile = glue('{PROJECT_location}/data/cancersea/hallmarks/{signature}.txt')
    genetab = read.table(genefile, header = T)
    
  }else{
    if(cancerType %in% cancerTypes){
      if(is.null(setname)){
        genefile = glue('{PROJECT_location}/data/hPPIN/Step3-Prune_DE2/{cancerType}/{cancerType}_{signature}_DE2')
        #print(genefile)
        genetab = read.table(genefile, header = T)
        genetab = unique(annotationdf[annotationdf$gene_name  %in% genetab$x,][,c('ENSG', 'gene_name')])
      }else{
        if(setname %in% c('discarded', 'retained', 'frost')){
          orifile = glue('{PROJECT_location}/data/cancersea/hallmarks/{signature}.txt')
          oritab = read.table(orifile, header = T)
          nsize = nrow(oritab)
          
          if(setname=='discarded'){
            genetab = oritab[!(oritab$ENSG %in% genetab$ENSG),] # dicarded
          }else if(setname=='retained'){
            genetab = oritab[(oritab$ENSG %in% genetab$ENSG),] # retained
          }else{
            tissueType = cancerType_to_hpa_tissue[[cancerType]]
            genetab  = read.table(glue('{PROJECT_location}/data/HPA/TissueSpecificWeightGeneration/computed_weights_{signature}.tsv'), sep="\t")
            
            genes = unique(rownames(genetab[genetab[[tissueType]] > 0.6927579, ])); 
            # genes = unique(rownames(genetab[order(genetab[[tissueType]], decreasing=T), ]))[1:nsize]
            genetab = unique(annotationdf[annotationdf$ENSG  %in% genes,][,c('ENSG', 'gene_name')])
          }
        }
      }
      
    }else if(cancerType=='random'){
      #set.seed(872436)           # Set seed
      #c_rand = sample(70:1200, 1) #large set
      c_rand = sample(70:200, 1) # small set
      x_rand <- sample(nrow(annotationdf))[1:c_rand]
      genetab = unique(annotationdf[x_rand,][,c('ENSG', 'gene_name')])
    }else{
      tissueType = cancerType
      genefile = glue("{PROJECT_location}/data/DE genes/{tissueType}_de2.txt")
      genetab =  read.table(genefile, sep="\t", header=TRUE)
      genetab = unique(annotationdf[annotationdf$gene_name  %in% genetab$x,][,c('ENSG', 'gene_name')])
    }
  }
  genetab$ENSG
}

#context_name = 'PRAD'
myargs = commandArgs(trailingOnly=TRUE)
context_name = myargs[1]

get_data_matrix <- function(context_name){
  if(context_name=='BRCA'){
    cdat = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/Breast/GSE75688_final_sample_information.txt"), sep="\t")
    #row.names(cdat) = cdat$X
    #cdat = cdat[,2:ncol(cdat)]
    countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/Breast/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"), sep="\t")
    mat = countdf[,cdat[cdat$type=='SC' & cdat$index3=='Tumor',]$sample] # only single cell tumors
    mat = countdf[,4:ncol(countdf)] # take bulk and single cell
    mat = countdf[,cdat[cdat$type=='SC',]$sample] # all single cells
    rownames(mat) = sapply(countdf$gene_id, function(x) strsplit(x, "[.]")[[1]][1])
    mat = mat[countdf$gene_type=='protein_coding',]
    
  }
  if(context_name=='COAD'){
    countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/COAD/GSE81861_CRC_tumor_all_cells_FPKM.csv"))
    mat = countdf[,1:ncol(countdf)]
    
    rownames(mat) = sapply(countdf$X, function(x) paste0('ENSG', strsplit(x, "ENSG")[[1]][2]))
    rownames(mat) = sapply(rownames(mat), function(x) strsplit(x, "[.]")[[1]][1])
    mat = mat[rownames(mat) %in% annotationdf$ENSG,]
    mat = mat[,sapply(colnames(mat), function(x) grepl('Epithelial|Fibroblast|Endothelial', x))]
  }
  if(context_name=='ESCA'){
    countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/ESCA/GSE81812_Normalized_counts.txt"))
    mat = countdf[,2:ncol(countdf)]
    
    rownames(mat) = sapply(countdf$X, function(x) strsplit(x, "[.]")[[1]][1])
    mat = mat[rownames(mat) %in% annotationdf$ENSG,]
    
  }
  if(context_name=='GBM'){
    countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/GBM/GSE84465_GBM_All_data.csv"), sep=" ")
    mat = countdf
    
    mat = mat[rownames(mat) %in% annotationdf$gene_name,]
    rownames(mat) = sapply(rownames(mat), function(x) annotationdf$ENSG[annotationdf$gene_name==x])
  }
  if(context_name=='LUAD'){
    # countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/LUAD/GSE69405_PROCESSED_GENE_TPM_ALL.txt"), sep='\t')
    # mat = countdf[,2:ncol(countdf)]
    # 
    # rownames(mat) = sapply(countdf$gene_id, function(x) paste0('ENSG', strsplit(x, "ENSG")[[1]][2]))
    # rownames(mat) = sapply(rownames(mat), function(x) strsplit(x, "[.]")[[1]][1])
    # mat = mat[rownames(mat) %in% annotationdf$ENSG,]
    # mat = mat[,sapply(colnames(mat), function(x) grepl('RHL', x))]
    
  }
  if(context_name=='PAAD'){
    countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/PAAD/GSE99305_Fishel.gene.counts.matrix.txt"), sep="\t")
    mat = countdf[,2:ncol(countdf)]
    
    rownames(mat) = countdf$gene
    mat = mat[rownames(mat) %in% annotationdf$ENSG,]
    mat = mat[,sapply(colnames(mat), function(x) grepl('tophat', x))]
    
  }
  if(context_name=='PRAD'){
    countdf = read.csv(glue("/Users/sharminm2/OneDrive - National Institutes of Health/context/data/cancerSEA/scRNASeq/PRAD/GSE99795_rpkm.txt"), sep="\t")
    mat = countdf[,9:ncol(countdf)]
    rownames(mat) = countdf$Transcript.RepeatID..cmd.analyzeRepeats.pl..data.resource.Homo_sapiens.Homo_sapiens.UCSC.hg19.Annotation.Archives.archive.2014.06.02.13.47.56.Genes.genes.gtf.none..condenseGenes..strand.both..count.exons..d.Sample_1_10.homer.Sample_1_11.homer.Sample_1_12.homer.Sample_1_13.homer.Sample_1_14.homer.Sample_1_15.homer.Sample_1_16.homer.Sample_1_17.homer.Sample_1_18.homer.Sample_1_19.homer.Sample_1_1.homer.Sample_1_20.homer.Sample_1_21.homer.Sample_1_22.homer.Sample_1_23.homer.Sample_1_24.homer.Sample_1_25.homer.Sample_1_26.homer.Sample_1_27.homer.Sample_1_28.homer.Sample_1_29.homer.Sample_1_2.homer.Sample_1_30.homer.Sample_1_31.homer.Sample_1_32.homer.Sample_1_33.homer.Sample_1_34.homer.Sample_1_35.homer.Sample_1_36.homer.Sample_1_37.homer.Sample_1_38.homer.Sample_1_39.homer.Sample_1_3.homer.Sample_1_40.homer.Sample_1_41.homer.Sample_1_42.homer.Sample_1_43.homer.Sample_1_44.homer.Sample_1_45.homer.Sample_1_46.homer.Sample_1_47.homer.Sample_1_48.homer.Sample_1_4.homer.Sample_1_5.homer.Sample_1_6.homer.Sample_1_7.homer.Sample_1_8.homer.Sample_1_9.homer.Sample_2_10.homer.Sample_2_11.homer.Sample_2_12.homer.Sample_2_13.homer.Sample_2_14.homer.Sample_2_15.homer.Sample_2_16.homer.Sample_2_17.homer.Sample_2_18.homer.Sample_2_19.homer.Sample_2_1.homer.Sample_2_20.homer.Sample_2_21.homer.Sample_2_22.homer.Sample_2_23.homer.Sample_2_24.homer.Sample_2_25.homer.Sample_2_26.homer.Sample_2_27.homer.Sample_2_28.homer.Sample_2_29.homer.Sample_2_2.homer.Sample_2_30.homer.Sample_2_31.homer.Sample_2_32.homer.Sample_2_33.homer.Sample_2_34.homer.Sample_2_35.homer.Sample_2_36.homer.Sample_2_37.homer.Sample_2_38.homer.Sample_2_39.homer.Sample_2_3.homer.Sample_2_40.homer.Sample_2_41.homer.Sample_2_42.homer.Sample_2_43.homer.Sample_2_44.homer.Sample_2_45.homer.Sample_2_46.homer.Sample_2_47.homer.Sample_2_48.homer.Sample_2_4.homer.Sample_2_5.homer.Sample_2_6.homer.Sample_2_7.homer.Sample_2_8.homer.Sample_2_9.homer.Sample_3_10.homer.Sample_3_11.homer.Sample_3_12.homer.Sample_3_13.homer.Sample_3_14.homer.Sample_3_15.homer.Sample_3_16.homer.Sample_3_17.homer.Sample_3_18.homer.Sample_3_19.homer.Sample_3_1.homer.Sample_3_20.homer.Sample_3_21.homer.Sample_3_22.homer.Sample_3_23.homer.Sample_3_24.homer.Sample_3_25.homer.Sample_3_26.homer.Sample_3_27.homer.Sample_3_28.homer.Sample_3_29.homer.Sample_3_2.homer.Sample_3_30.homer.Sample_3_31.homer.Sample_3_32.homer.Sample_3_33.homer.Sample_3_34.homer.Sample_3_35.homer.Sample_3_36.homer.Sample_3_37.homer.Sample_3_38.homer.Sample_3_39.homer.Sample_3_3.homer.Sample_3_40.homer.Sample_3_41.homer.Sample_3_42.homer.Sample_3_43.homer.Sample_3_44.homer.Sample_3_45.homer.Sample_3_46.homer.Sample_3_47.homer.Sample_3_48.homer.Sample_3_4.homer.Sample_3_5.homer.Sample_3_6.homer.Sample_3_7.homer.Sample_3_8.homer.Sample_3_9.homer.Sample_D_1.homer.Sample_D_2.homer.Sample_D_3.homer..rpkm.
    mat = mat[sapply(rownames(mat), function(x) grepl('NM_', x)), ]
    cols = c("SYMBOL", "ENSEMBL")
    genemap = select(org.Hs.eg.db, keys=rownames(mat), columns=cols, keytype="REFSEQ")
    genemap = genemap[genemap$ENSEMBL %in% annotationdf$ENSG,]
    
    genes = unlist(sapply(rownames(mat), function(x) {
      genemap[genemap$REFSEQ==x,]$ENSEMBL
    }))
    
    mat = mat[rownames(mat) %in% names(genes),]
    genes = genes[names(genes) %in% rownames(mat)]
    
    mat = mat[match(names(genes), rownames(mat)),]
    #mat = mat[order(names(genes)),]
    uniques = !duplicated(genes)
    genes = genes[uniques]
    mat = mat[uniques,]
    rownames(mat) = sapply(rownames(mat), function(x) genes[x])
    
  }
  
  # sce
  return(mat)
}


# check for activity of random set of genes: larger geneset size
gs_random_smaller = lapply(1:1000, function(x) read_genes(x, 'random'))
names(gs_random_smaller) = 1:1000
gsva.es_random_smaller <- gsva(as.matrix(mat), gs_random_smaller, verbose=FALSE)
summary(rowMaxs(gsva.es_random_smaller))
summary(rowSums(gsva.es_random_smaller>0.1)/ncol(gsva.es_random_smaller))
summary(unlist(lapply(gs_random_smaller, length)))

# do we see any significant score distribution between two es of random where one set has more: smaller geneset size
gs_random_larger = lapply(1:1000, function(x) read_genes(x, 'random'))
names(gs_random_larger) = 1:1000
gsva.es_random_larger <- gsva(as.matrix(mat), gs_random_larger, verbose=FALSE)
summary(rowMaxs(gsva.es_random_larger))
summary(rowSums(gsva.es_random_larger>0.1)/ncol(gsva.es_random_larger))
summary(unlist(lapply(gs_random_larger, length)))


summary(rowMaxs(gsva.es_random_smaller, useNames = TRUE))
summary(rowMaxs(gsva.es_random_larger, useNames = TRUE))

summary(rowMeans(gsva.es_random_smaller))
summary(rowMeans(gsva.es_random_larger))

#significantly different
summary(rowMedians(gsva.es_random_smaller, useNames = TRUE))
summary(rowMedians(gsva.es_random_larger, useNames = TRUE))

#significantly different
summary(rowSums(gsva.es_random_smaller>0.1)/ncol(gsva.es_random_smaller))
summary(rowSums(gsva.es_random_larger>0.1)/ncol(gsva.es_random_larger))

summary(unlist(lapply(gs_random_smaller, length)))
summary(unlist(lapply(gs_random_larger, length)))

matlist = list()
for(cancerType in cancerTypes){
  print(cancerType)
  matlist[[cancerType]] = get_data_matrix(cancerType)
}
lapply(matlist, dim)


# original
gsva.es_list = lapply(cancerTypes, function(x){
  current_mat = matlist[[x]]
  reshape_mat <- sapply(current_mat, as.numeric)
  reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                       dimnames=list(rownames(current_mat), colnames(current_mat)))
  
  gs = lapply(signatures, read_genes)
  names(gs) = signatures
  gsvapar = gsvaParam(reshape_mat, gs, kcdf="Poisson")
  gsva.es = gsva(gsvapar)
  #gsva.es = gsva(as.matrix(current_mat), gs, verbose=FALSE)
  return(gsva.es)
})
names(gsva.es_list) = cancerTypes
  
# context specific
gsva.es_context_list = lapply(cancerTypes, function(x){
  current_mat = matlist[[x]]
  reshape_mat <- sapply(current_mat, as.numeric)
  reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                       dimnames=list(rownames(current_mat), colnames(current_mat)))
  
  gs_context = lapply(signatures, function(y) read_genes(y, x))
  names(gs_context) = signatures
  gsvapar = gsvaParam(reshape_mat, gs_context, kcdf="Poisson")
  gsva.es = gsva(gsvapar)
  #gsva.es = gsva(as.matrix(current_mat), gs_context, verbose=FALSE)
  return(gsva.es)
})
names(gsva.es_context_list) = cancerTypes

# de
gsva.es_de_list = lapply(cancerTypes, function(x){
  current_mat = matlist[[x]]
  reshape_mat <- sapply(current_mat, as.numeric)
  reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                       dimnames=list(rownames(current_mat), colnames(current_mat)))
  
  tissueType = cancerType_to_tissue[[x]]
  gs_de = list()
  gs_de[['DE']] = read_genes(signatures[1], tissueType)
  gsvapar = gsvaParam(reshape_mat, gs_de, kcdf="Poisson")
  gsva.es = gsva(gsvapar)
  #gsva.es_de <- gsva(as.matrix(current_mat), gs_de, verbose=FALSE)
  return(gsva.es)
})
names(gsva.es_de_list) = cancerTypes

#discard
gsva.es_discard_list = lapply(cancerTypes, function(x){
  current_mat = matlist[[x]]
  reshape_mat <- sapply(current_mat, as.numeric)
  reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                       dimnames=list(rownames(current_mat), colnames(current_mat)))
  
  gs_discard = lapply(signatures, function(y) read_genes(y, x, 'discarded'))
  names(gs_discard) = signatures
  gsvapar = gsvaParam(reshape_mat, gs_discard, kcdf="Poisson")
  gsva.es = gsva(gsvapar)
  #gsva.es = gsva(as.matrix(current_mat), gs_context, verbose=FALSE)
  return(gsva.es)
})
names(gsva.es_discard_list) = cancerTypes

#retained
gsva.es_retain_list = lapply(cancerTypes, function(x){
  current_mat = matlist[[x]]
  reshape_mat <- sapply(current_mat, as.numeric)
  reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                       dimnames=list(rownames(current_mat), colnames(current_mat)))
  
  gs_retain = lapply(signatures, function(y) read_genes(y, x, 'retained'))
  names(gs_retain) = signatures
  gsvapar = gsvaParam(reshape_mat, gs_retain, kcdf="Poisson")
  gsva.es = gsva(gsvapar)
  #gsva.es = gsva(as.matrix(current_mat), gs_context, verbose=FALSE)
  return(gsva.es)
})
names(gsva.es_retain_list) = cancerTypes

# frost
gsva.es_frost_list = lapply(intersect(names(cancerType_to_hpa_tissue), cancerTypes), function(x){
  print(x)
  current_mat = matlist[[x]]
  reshape_mat <- sapply(current_mat, as.numeric)
  reshape_mat = matrix(reshape_mat, nrow=nrow(current_mat), ncol=ncol(current_mat), 
                       dimnames=list(rownames(current_mat), colnames(current_mat)))
  
  gs_frost = lapply(signatures, function(y) read_frost_genes(y, x))
  names(gs_frost) = signatures
  gsvapar = gsvaParam(reshape_mat, gs_frost, kcdf="Poisson")
  gsva.es = gsva(gsvapar)
  #gsva.es = gsva(as.matrix(current_mat), gs_context, verbose=FALSE)
  return(gsva.es)
})
names(gsva.es_frost_list) = intersect(names(cancerType_to_hpa_tissue), cancerTypes)

saveRDS(list(reference=gsva.es_list, context=gsva.es_context_list, de=gsva.es_de_list, 
             discard=gsva.es_discard_list, retain=gsva.es_retain_list, frost=gsva.es_frost_list), 
        file='/Users/sharminm2/OneDrive - National Institutes of Health/context/data/results/hPPIN_DE2/gsva/poission.RDS')
mydata = readRDS('/Users/sharminm2/OneDrive - National Institutes of Health/context/data/results/hPPIN_DE2/gsva/poission.RDS')
gsva.es_list = mydata[[1]] 
gsva.es_context_list=mydata[[2]]
gsva.es_de_list=mydata[[3]]
gsva.es_discard_list=mydata[[4]]
gsva.es_retain_list=mydata[[5]]
gsva.es_frost_list=mydata[[6]]


# counts plot
count_list = lapply(cancerTypes, function(x){
  gsva.es = gsva.es_list[[x]]
  gsva.es_context = gsva.es_context_list[[x]]
  k = 0.1
  agnostic = rowSums(gsva.es>k)/ncol(gsva.es)
  projected = rowSums(gsva.es_context>k)/ncol(gsva.es)
  df = data.frame(agnostic, projected)
  df['hallmark'] = rownames(df)
  df = melt(df)
  colnames(df) = c('hallmark', 'Context', 'value')
  #df$Context[df$Context == 'projected'] <- x
  levels(df$Context) = c('agnostic', x)
  p = ggplot(data=df, aes(x=hallmark, y=value, fill=Context)) +
    geom_bar(stat="identity", position=position_dodge())+
    theme_pubr(base_size = 13.5, legend='right')+ labs_pubr()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))+
    ylab("Fraction of cells\nwith\nfunctional activity > 0.1")+
    xlab('') #print(p)
  
  return(p)
})
names(count_list) = cancerTypes


plot_grid(count_list[['COAD']], count_list[['GBM']], rel_widths = c(3,3), 
              labels = c('', ''), label_size = 12, ncol = 1)
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3c_coad_gbm.png'), height=10, width=8)


outfile = glue('{PROJECT_location}/Manuscript/Supplementary/Fig_S3 - ssgsva_count.pdf')
pdf(outfile, width=10, height=5, onefile = TRUE)
for (i in names(count_list)) {
  if(i=='COAD'|i=='GBM'){next}
  print(i)
  print(count_list[[i]])
}
dev.off()

# # context/reference 
logcount_list = lapply(cancerTypes, function(x){
  gsva.es = gsva.es_list[[x]]
  gsva.es_context = gsva.es_context_list[[x]]
  k = 0.1
  agnostic = rowSums(gsva.es>k)/ncol(gsva.es)
  projected = rowSums(gsva.es_context>k)/ncol(gsva.es)
  df = data.frame(agnostic, projected)
  df['hallmark'] = rownames(df)
  df['logFC'] = log2(df$projected/df$agnostic)
  df['cancerType'] = x
  return(df)
})
df = do.call(rbind, logcount_list)
b  = c(2,1, 0, -1,-2)
limits = c(min(df$logFC), max(df$logFC))
ggplot(df, aes(x=cancerType, y=hallmark, fill=logFC)) + 
  geom_tile(show.legend=T)+ 
  scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
  labs_pubr()+#theme_pubr(base_size=5) +
  theme(text = element_text(size=18))#+
  #ggtitle('LOG(#active cells by CATS-geneset / #active cells by REF-geneset)')
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3c_logcounts.png'))

# context/de
logcount_list = lapply(cancerTypes, function(x){
  gsva.es = gsva.es_list[[x]]
  gsva.es_context = gsva.es_context_list[[x]]
  gsva.es_de = gsva.es_de_list[[x]]
  k = 0.1
  agnostic = rowSums(gsva.es>k)/ncol(gsva.es)
  projected = rowSums(gsva.es_context>k)/ncol(gsva.es)
  dexp = rowSums(gsva.es_de>k)/ncol(gsva.es)
  df = data.frame(projected, dexp)
  df['hallmark'] = rownames(df)
  df['logFC'] = log2(df$projected/df$dexp)
  df['cancerType'] = x
  return(df)
})
df = do.call(rbind, logcount_list)
b  = c(2,1, 0, -1,-2)
limits = c(min(df$logFC), max(df$logFC))
ggplot(df, aes(x=cancerType, y=hallmark, fill=logFC)) + 
  geom_tile(show.legend=T)+ 
  scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
  labs_pubr()+#theme_pubr(base_size=5) +
  theme(text = element_text(size=18))#+
  #ggtitle('LOG(#active cells by CATS-geneset / #active cells by DE-geneset)')
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3c_logcounts2.png'))

# retain/discard, cancerTypes
logcount_list = lapply(cancerTypes, function(x){
  print(x)
  gsva.es = gsva.es_list[[x]]
  gsva.es_discard = gsva.es_discard_list[[x]]
  gsva.es_retain = gsva.es_retain_list[[x]]
  k = 0.1
  retain = rowSums(gsva.es_retain>k)/ncol(gsva.es)
  discard = rowSums(gsva.es_discard>k)/ncol(gsva.es)
  s = intersect(names(retain), names(discard))
  df = data.frame(discard[s], retain[s])
  df['hallmark'] = rownames(df)
  df['logFC'] = log2(df$retain/df$discard)
  df['cancerType'] = x
  
  return(df)
})
df = do.call(rbind, logcount_list)
b  = c(2,1, 0, -1,-2)
limits = c(min(df$logFC), max(df$logFC))
ggplot(df, aes(x=cancerType, y=hallmark, fill=logFC)) + 
  geom_tile(show.legend=T)+ 
  scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
  labs_pubr()+#theme_pubr(base_size=5) +
  theme(text = element_text(size=18))#+
#ggtitle('LOG(#active cells by CATS-geneset / #active cells by DE-geneset)')
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3c_logcounts3.png'))


# context/frost
logcount_list = lapply(intersect(names(cancerType_to_hpa_tissue), cancerTypes), function(x){
  gsva.es = gsva.es_list[[x]]
  gsva.es_context = gsva.es_context_list[[x]]
  gsva.es_hpa = gsva.es_frost_list[[x]]
  k = 0.1
  hpa = rowSums(gsva.es_hpa>k)/ncol(gsva.es)
  projected = rowSums(gsva.es_context>k)/ncol(gsva.es)
  df = data.frame(hpa, projected)
  df['hallmark'] = rownames(df)
  df['logFC'] = log2(df$projected/df$hpa)
  df['cancerType'] = x
  return(df)
})
df = do.call(rbind, logcount_list)
b  = c(2,1, 0, -1,-2)
limits = c(min(df$logFC), max(df$logFC))
ggplot(df, aes(x=cancerType, y=hallmark, fill=logFC)) + 
  geom_tile(show.legend=T)+ 
  scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
  labs_pubr()+#theme_pubr(base_size=5) +
  theme(text = element_text(size=18))#+
  #ggtitle('LOG(#active cells by CATS-geneset / #active cells by HPA-geneset)')
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3c_logcounts4.png'))


# heatplot
heat_list = lapply(cancerTypes, function(x){
  gsva.es = gsva.es_list[[x]]
  gsva.es_context = gsva.es_context_list[[x]]
  gsva.es_discard = gsva.es_discard_list[[x]]
  gsva.es_de = gsva.es_de_list[[x]]
  cell_order = order(gsva.es_de, decreasing = TRUE)
  tissueType = cancerType_to_tissue[[x]]
  b  = c(3,2,1,0,-1,-2,-3)
  limits = c(min(gsva.es_de, gsva.es, gsva.es_context), max(gsva.es_de, gsva.es, gsva.es_context))
  # original
  df = melt(gsva.es[,cell_order], na.rm = FALSE, value.name='value')
  colnames(df) = c('Signatures', 'Cells', 'activity')
  p1 = ggplot(df, aes(x=Cells, y=Signatures, fill=activity)) + 
    geom_tile(show.legend=F)+
    scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
    labs_pubr()+#theme_pubr(base_size=5) +
    theme(text = element_text(size=18),
          axis.title.x=element_blank(), axis.text.x=element_blank(), 
          axis.ticks.x=element_blank())+ylab('')+
    ggtitle('Context-agnostic') 
  
  # context
  df = melt(gsva.es_context[,cell_order], na.rm = FALSE, value.name='value')
  colnames(df) = c('Signatures', 'Cells', 'activity')
  p2 = ggplot(df, aes(x=Cells, y=Signatures, fill=activity)) + 
    geom_tile(show.legend=F)+ 
    scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
    labs_pubr()+#theme_pubr(base_size=5) +
    theme(text = element_text(size=18),
          axis.title.x=element_blank(), axis.text.x=element_blank(), 
          axis.ticks.x=element_blank())+ylab('')+
    ggtitle(glue('{x}-specific'))
  
  # discard
  df = melt(gsva.es_discard[,cell_order], na.rm = FALSE, value.name='value')
  colnames(df) = c('Signatures', 'Cells', 'activity')
  p2 = ggplot(df, aes(x=Cells, y=Signatures, fill=activity)) + 
    geom_tile(show.legend=F)+ 
    scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
    labs_pubr()+#theme_pubr(base_size=5) +
    theme(text = element_text(size=18),
          axis.title.x=element_blank(), axis.text.x=element_blank(), 
          axis.ticks.x=element_blank())+ylab('')+
    ggtitle(glue('{x}-specific'))
  
  # non context: de
  df = melt(t(gsva.es_de[,cell_order]), na.rm = FALSE, value.name='value')
  colnames(df) = c('DE genes', 'Cells', 'activity')
  #df['DE genes'] = '                  '
  df['DE genes'] = "                       "; head(df) #equivalent to angiogenesis
  p3 = ggplot(df, aes(x=Cells, y=`DE genes`, fill=activity)) + 
    geom_tile()+ 
    scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
    labs_pubr()+#theme_pubr(base_size=5) +
    theme(legend.position="bottom", text = element_text(size=18),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          #axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
    )+ylab('')+ggtitle(paste(" ", tissueType, 'DE')) 
  
  
  #top_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12, ncol = 1)
  p = plot_grid(p1, p2, p3, rel_widths = c(3,3,0.01), 
                labels = c('', '', ''), label_size = 12, ncol = 1); print(p)
  #figure <- ggarrange(p1, p2, p3, nrow = 3, widths = c(5, 5, 0.1),
  #                    align = "v", common.legend = TRUE)
  #p = annotate_figure(figure, fig.lab = x)
  return(p)
})
names(heat_list) = cancerTypes

heat_list[['BRCA']]
ggsave2(glue('{PROJECT_location}/plots/{CUR_DST}/Fig_3b_brca.png'))

outfile = glue('{PROJECT_location}/Manuscript/Supplementary/Fig_S2 - ssgsva_heatmap.pdf')
pdf(outfile, width=10, height=8, onefile = TRUE)
for (i in names(heat_list)[2:length(cancerTypes)]) {
  print(i)
  print(heat_list[[i]])
}
dev.off()

# gsva
run_gsva = function(context_name, plottype){
  
  
  
  if(plottype=='correlation'){
    print('correlation')
    m = cor(t(gsva.es))
    p1 = ggplot(melt(m), aes(Var1, Var2, fill=value)) + 
      geom_tile(show.legend=F) +
      scale_fill_distiller(palette = "RdPu", limits=c(-1,1)) +
      #theme_ipsum()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
    
    m = cor(t(gsva.es_context))
    p2 = ggplot(melt(m), aes(Var1, Var2, fill=value)) + 
      geom_tile() +
      scale_fill_distiller(palette = "RdPu", limits=c(-1,1)) +
      #theme_ipsum()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
    
    p = plot_grid(p1, p2, labels=c('A (Original)', 'B (Projected)'), rel_widths=c(4,5), label_size=12)+
      ggtitle(context_name)
    return(p)
  }
  
  # non context: de
  tissueType = cancerType_to_tissue[[context_name]]
  if(plottype=='boxplot'){
    gs_de = lapply(signatures, function(x) read_genes(x, tissueType))
    names(gs_de) = signatures
  }else{
    gs_de = list()
    gs_de[['DE']] = read_genes(signatures[1], tissueType)
  }
  
  gsva.es_de <- gsva(as.matrix(mat), gs_de, verbose=FALSE)
  cell_order = order(gsva.es_de, decreasing = TRUE)
  
  # non context: other cancers
  gsva.es_noncontext = list()
  # for(cancerType in cancerTypes){
  #   if(cancerType==context_name){
  #     next
  #   }
  #   gs_noncontext = lapply(signatures, function(x) {
  #     g = read_genes(x, cancerType)
  #     g = setdiff(g, gs_context[x])
  #     g
  #   })
  #   names(gs_noncontext) = signatures
  #   
  #   gsva.es_noncontext[[cancerType]] <- gsva(as.matrix(mat), gs_noncontext, verbose=FALSE)
  # }
  
  if(plottype=='boxplot'){
    print('boxplot')
    df_original = melt(gsva.es, na.rm = FALSE, value.name='value')[,c('Var1', 'value')]
    df_original$L1 = 'Original'
    df_context = melt(gsva.es_context, na.rm = FALSE, value.name='value')[,c('Var1', 'value')]
    df_context$L1 = 'Projected'
    df_de = melt(gsva.es_de, na.rm = FALSE, value.name='value')[,c('Var1', 'value')]
    df_de$L1 = 'DE'
    df_noncontext = melt(gsva.es_noncontext, na.rm = FALSE, value.name='value')[,c('Var1', 'value', 'L1')]
    df = do.call('rbind', list(df_original, df_context, df_de, df_noncontext))
    df$L1 <- factor(df$L1, 
                    levels = c('Original','Projected', 'DE', sort(cancerTypes)),
                    ordered = TRUE)
    colnames(df) = c('Signatures', 'Functional_activity', 'Genes')
    
    p = ggplot(df, aes(x=Signatures, y= Functional_activity, fill=Genes)) + 
      geom_boxplot()+ theme_bw()+facet_wrap(~Signatures, scale="free")+
      scale_fill_brewer(palette="Set2")+
    ggtitle(context_name)
    return(p)
  }else if(plottype=='heatmap'){
    print('heatmap')
    p_list = list()
    b  = c(3,2,1,0,-1,-2,-3)
    limits = c(min(gsva.es_de, gsva.es, gsva.es_context), max(gsva.es_de, gsva.es, gsva.es_context))
    # original
    df = melt(gsva.es[,cell_order], na.rm = FALSE, value.name='value')
    colnames(df) = c('Signatures', 'Cells', 'activity')
    p1 = ggplot(df, aes(x=Cells, y=Signatures, fill=activity)) + 
      geom_tile(show.legend=F)+
      scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
      theme(#legend.position="bottom",
        axis.title.x=element_blank(), axis.text.x=element_blank(), 
            axis.ticks.x=element_blank())+
      ggtitle('Original') 
    
    # context
    df = melt(gsva.es_context[,cell_order], na.rm = FALSE, value.name='value')
    colnames(df) = c('Signatures', 'Cells', 'activity')
    p2 = ggplot(df, aes(x=Cells, y=Signatures, fill=activity)) + 
      geom_tile(show.legend=F)+ 
      scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
      theme(#legend.position="bottom",
            axis.title.x=element_blank(), axis.text.x=element_blank(), 
            axis.ticks.x=element_blank())+
      ggtitle('Context')
    
    # non context: de
    df = melt(t(gsva.es_de[,cell_order]), na.rm = FALSE, value.name='value')
    colnames(df) = c('DE genes', 'Cells', 'activity')
    #df['DE genes'] = '                     '
    p3 = ggplot(df, aes(x=Cells, y=`DE genes`, fill=activity)) + 
      geom_tile()+ 
      scale_fill_gradient2(high='red', low='blue', breaks=b, labels=format(b), midpoint=0, limits=limits)+
      theme(legend.position="bottom",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
        )+
      ggtitle(paste(" ", tissueType, 'DE')) 
    
    
    #top_row <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12, ncol = 1)
    p = plot_grid(p1, p2, p3, rel_widths = c(3,3,0.1),
                  labels = c('A', 'B', 'C'), label_size = 12, ncol = 1)
    #figure <- ggarrange(p1, p2, p3, nrow = 3, widths = c(5, 5, 1),
    #                    align = "v", common.legend = TRUE)
    #p = annotate_figure(figure, fig.lab = context_name)
    return(p)
    
  }
  
  
}

p_list = lapply(cancerTypes, function(x) run_gsva(x, 'counts'))
names(p_list) = cancerTypes

p_list = lapply(cancerTypes, function(x) run_gsva(x, 'heatmap'))
names(p_list) = cancerTypes

outfile = glue('{PROJECT_location}/plots/{CUR_DST}/ssgsva_heatmap.pdf')
pdf(outfile, width=10, height=5, onefile = TRUE)
for (i in names(p_list)) {
  print(i)
  print(p_list[[i]])
}
dev.off()





