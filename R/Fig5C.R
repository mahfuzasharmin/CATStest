library(GEOquery)
library(glue)
#library(xlsx) # the installation went wrong and cannot be used now
library(limma)
library(gplots)
library("RColorBrewer")
library("viridis")
library(ggplot2)
#library(cowplot) # could not use
#library(ggpubr) # could not use
#library(ComplexHeatmap) # could not use
library(gridGraphics) # actually used
library(grid)
library(gridExtra)


library(ggplot2)
library(reshape2)


USERNAME = "sharminm2"
DB_location = glue("/Users/{USERNAME}/Dropbox/nih/db")
RESDIR = glue("Dropbox/nih/context/data/results/hPPIN_DE2/responder_data/")
PROJECT_location = glue("/Users/{USERNAME}/Dropbox/nih/context")
GENE_NETWORK = 'hPPIN'
DE_K = '2'
CUR_DST = glue("{GENE_NETWORK}_DE{DE_K}")
extra = ''
if(GENE_NETWORK=='STRING'){
  extra = 'String_'
}
annotationdf = get_annotation()


#--------------------GSE19697----------
gseno = 'GSE19697'
gse <- getGEO(gseno,GSEMatrix=TRUE)

p = pData(phenoData(gse$GSE19697_series_matrix.txt.gz))
p = p[,c('characteristics_ch1.9', 'recur:ch1')]
table(p$`characteristics_ch1.9`, p$`recur:ch1`)

e = exprs(gse$GSE19697_series_matrix.txt.gz)
f = pData(featureData(gse$GSE19697_series_matrix.txt.gz))

protein_coding = f[f$`Gene Symbol` %in% annotationdf$gene_name,c('ID', 'Gene Symbol')]
protein_e = e[protein_coding$ID,]
rownames(protein_e) = protein_coding[rownames(protein_e), 'Gene Symbol']

reduced_e = lapply(unique(rownames(protein_e)), function(x){
  sube = protein_e[rownames(protein_e)==x,]
  if(is.matrix(sube)){
    colMeans(sube)
  }else{
    sube
  }
})
names(reduced_e) = unique(rownames(protein_e))
reduced_e = t(data.frame(reduced_e))
y <- data.frame(normalizeQuantiles(log2(reduced_e)))
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', rownames(p)[p$`recur:ch1`=="0"])]
nonresponder = y[, c('Genes', rownames(p)[p$`recur:ch1`=="1"])]

write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------GSE123728--------
gseno = 'GSE123728'
gse <- getGEO(gseno,GSEMatrix=TRUE)

p = pData(phenoData(gse$GSE123728_series_matrix.txt.gz))
p = p[,c('characteristics_ch1.1', 'recurrence:ch1')]
table(p$`characteristics_ch1.1`, p$`recurrence:ch1`)

e = exprs(gse$GSE123728_series_matrix.txt.gz)
f = pData(featureData(gse$GSE123728_series_matrix.txt.gz))

protein_coding = f[f$`ORF` %in% annotationdf$gene_name,c('ID', 'ORF')]
protein_e = e[protein_coding$ID,]
rownames(protein_e) = protein_coding[rownames(protein_e), 'ORF']

reduced_e = lapply(unique(rownames(protein_e)), function(x){
  sube = protein_e[rownames(protein_e)==x,]
  if(is.matrix(sube)){
    colMeans(sube)
  }else{
    sube
  }
})
names(reduced_e) = unique(rownames(protein_e))
reduced_e = t(data.frame(reduced_e))
y <- data.frame(normalizeQuantiles(log2(reduced_e+0.001)))
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', rownames(p)[p$`recurrence:ch1`=="0"])]
nonresponder = y[, c('Genes', rownames(p)[p$`recurrence:ch1`=="1"])]

write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------GSE66999------------
gseno = 'GSE66999'
gse <- getGEO(gseno,GSEMatrix=TRUE)

p = pData(phenoData(gse$GSE66999_series_matrix.txt.gz)); head(p)
p = p[,c('characteristics_ch1.3', 'clinical response:ch1')]
table(p$characteristics_ch1.3, p$`clinical response:ch1`)

e = exprs(gse$GSE66999_series_matrix.txt.gz)
f = pData(featureData(gse$GSE66999_series_matrix.txt.gz))

protein_coding = f[f$`GENE_SYMBOL` %in% annotationdf$gene_name,c('ID', 'GENE_SYMBOL')]
protein_e = e[protein_coding$ID,]
rownames(protein_e) = protein_coding[rownames(protein_e), 'GENE_SYMBOL']

reduced_e = lapply(unique(rownames(protein_e)), function(x){
  sube = protein_e[rownames(protein_e)==x,]
  if(is.matrix(sube)){
    colMeans(sube)
  }else{
    sube
  }
})
names(reduced_e) = unique(rownames(protein_e))
reduced_e = t(data.frame(reduced_e))
y <- data.frame(normalizeQuantiles(log2(reduced_e+0.001)))
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', rownames(p)[p$`clinical response:ch1` %in% 
                                         c("Complete Response (CR)", "Partial Response (PR)")])]
nonresponder = y[, c('Genes', rownames(p)[p$`clinical response:ch1` %in% 
                                            c("Stable Disease (SD)", "Progressive Disease (PD)")])]

write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------GSE25066------------
gseno = 'GSE25066'
gse <- getGEO(gseno,GSEMatrix=TRUE)

p = pData(phenoData(gse$GSE25066_series_matrix.txt.gz)); head(p)
p = p[,c('chemosensitivity_prediction:ch1', 'pam50_class:ch1')]
table(p$`chemosensitivity_prediction:ch1`, p$`pam50_class:ch1`)

e = exprs(gse$GSE25066_series_matrix.txt.gz)
f = pData(featureData(gse$GSE25066_series_matrix.txt.gz))

protein_coding = f[f$`Gene Symbol` %in% annotationdf$gene_name,c('ID', 'Gene Symbol')]
protein_e = e[protein_coding$ID,]
rownames(protein_e) = protein_coding[rownames(protein_e), 'Gene Symbol']

reduced_e = lapply(unique(rownames(protein_e)), function(x){
  sube = protein_e[rownames(protein_e)==x,]
  if(is.matrix(sube)){
    colMeans(sube)
  }else{
    sube
  }
})
names(reduced_e) = unique(rownames(protein_e))
reduced_e = t(data.frame(reduced_e))
y <- data.frame(normalizeQuantiles(log2(reduced_e+0.001)))
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', rownames(p)[p$`chemosensitivity_prediction:ch1`=="Rx Sensitive"])]
nonresponder = y[, c('Genes', rownames(p)[p$`chemosensitivity_prediction:ch1`=="Rx Insensitive"])]

write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

for(i in c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')){
  responder = y[, c('Genes', rownames(p)[p$`chemosensitivity_prediction:ch1`=="Rx Sensitive" & p$`pam50_class:ch1`==i])]
  nonresponder = y[, c('Genes', rownames(p)[p$`chemosensitivity_prediction:ch1`=="Rx Insensitive" & p$`pam50_class:ch1`==i])]
  write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder_{i}'), sep="\t", row.names=F, quote = F)
  write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder_{i}'), sep="\t", row.names=F, quote = F)
  
}

#--------------------GSE91061------------
gseno = 'GSE91061'
clinical = read.xlsx('{RESDIR}/cinicalData.MasterFile.xlsx', sheetName = gseno)
sampledat = read.xlsx('{RESDIR}/sample.data.mapping.MasterFile.xlsx', sheetName = gseno)
sampledat = sampledat[sampledat$Treatment=='PRE',]
p = merge(clinical, sampledat, by = 'Patient')


e = read.csv(glue('{RESDIR}/salmon.tpm.{gseno}.txt'), sep='\t')
f = rownames(e)

protein_coding = f[f %in% annotationdf$gene_name]
protein_e = e[protein_coding,]

y <- protein_e
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PR', 'CR')], colnames(y)))]
nonresponder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PD', 'SD')], colnames(y)))]

write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------GSE78220------------
gseno = 'GSE78220'
clinical = read.xlsx('{RESDIR}/cinicalData.MasterFile.xlsx', sheetName = gseno)
clinical = na.omit(clinical[,c('Patient', 'Cohort', 'RECIST', 'Status', 'OS.days.')])
sampledat = read.xlsx('{RESDIR}/sample.data.mapping.MasterFile.xlsx', sheetName = gseno)
sampledat = sampledat[sampledat$Treatment=='PRE',]
p = merge(clinical, sampledat, by = 'Patient')

e = read.csv(glue('{RESDIR}/salmon.tpm.{gseno}.txt'), sep='\t')
f = rownames(e)

protein_coding = f[f %in% annotationdf$gene_name]
protein_e = e[protein_coding,]

y <- protein_e
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PR', 'CR')], colnames(y)))]
nonresponder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PD', 'SD')], colnames(y)))]

write.table(responder, glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------PRJEB23709, log2qn------------
gseno = 'PRJEB23709'
clinical = read.table(glue('{RESDIR}/cinicalData.MasterFile.{gseno}.txt'), sep='\t', header=T)
sampledat = read.table(glue('{RESDIR}/sample.data.mapping.MasterFile.{gseno}.txt'), sep='\t', header=T)
sampledat = sampledat[sampledat$Treatment=='PRE',]
p = merge(clinical, sampledat, by = 'Patient')
p = p[p$Drug=='ipiPD1',] # PD1 or ipiPD1
e = read.csv(glue('{RESDIR}/salmon.counts.{gseno}.txt'), sep='\t')
f = rownames(e)

protein_coding = f[f %in% annotationdf$gene_name]
protein_e = e[protein_coding,]

y <- normalizeQuantiles(log2(protein_e+0.001))
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PR', 'CR')], colnames(y)))]
nonresponder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PD', 'SD')], colnames(y)))]

write.table(responder, glue('{RESDIR}/{gseno}_ipiPD1_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_ipiPD1_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------PRJEB23709, tpm------------
gseno = 'PRJEB23709'
clinical = read.xlsx(glue('{RESDIR}/cinicalData.MasterFile.xlsx'), sheetName = gseno)
sampledat = read.xlsx(glue('{RESDIR}/sample.data.mapping.MasterFile.xlsx'), sheetName = gseno)
sampledat = sampledat[sampledat$Treatment=='PRE',]
p = merge(clinical, sampledat, by = 'Patient')
p = p[p$Drug=='ipiPD1',] # PD1 or ipiPD1
e = read.csv(glue('{RESDIR}/salmon.tpm.{gseno}.txt'), sep='\t')
f = rownames(e)

protein_coding = f[f %in% annotationdf$gene_name]
protein_e = e[protein_coding,]

y <- protein_e
samples = colnames(y)
y['Genes'] = rownames(y)
y = y[, c('Genes', samples)]
responder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PR', 'CR')], colnames(y)))]
nonresponder = y[, c('Genes', intersect(p$Sample[p$RECIST %in% c('PD', 'SD')], colnames(y)))]

write.table(responder, glue('{RESDIR}/{gseno}_tpm_ipiPD1_QN_Responder'), sep="\t", row.names=F, quote = F)
write.table(nonresponder, glue('{RESDIR}/{gseno}_tpm_ipiPD1_QN_Non-Responder'), sep="\t", row.names=F, quote = F)

#--------------------effectsize of PRJEB23709----------
gseno = 'PRJEB23709'
responder = read.table(glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", header = TRUE)
rownames(responder) = responder$Genes
responder = responder[,2:ncol(responder)]
responder[1:5,1:5]

non_responder = read.table(glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", header = TRUE)
rownames(non_responder) = non_responder$Genes
non_responder = non_responder[,2:ncol(non_responder)]
non_responder[1:5,1:5]

#--------------------effectsize of GSE25066----------
gseno = 'GSE25066'
responder = read.table(glue('{RESDIR}/{gseno}_QN_Responder'), sep="\t", header = TRUE)
rownames(responder) = responder$Genes
responder = responder[,2:ncol(responder)]
responder[1:5,1:5]

non_responder = read.table(glue('{RESDIR}/{gseno}_QN_Non-Responder'), sep="\t", header = TRUE)
rownames(non_responder) = non_responder$Genes
non_responder = non_responder[,2:ncol(non_responder)]
non_responder[1:5,1:5]

#--------------- now do the plot --------------
labs = c(rep('responder', ncol(responder)), rep('non_responder', ncol(non_responder)))
df = cbind(responder, non_responder)
#scaling
df2 = t(apply(df, 1, scale))
dim(df)

load_genelist <- function(cancerType, signature){
  if(cancerType=='original'){
    genefile = glue("{PROJECT_location}/data/cancersea/hallmarks/{signature}.txt")
    genetab =  read.table(genefile, sep="\t", header=TRUE)
    geneList = genetab$Symbol
  }else{
    genefile = glue("{PROJECT_location}/data/{GENE_NETWORK}/Step3-Prune_DE{DE_K}/{cancerType}/{cancerType}_{signature}_{extra}DE{DE_K}")
    
    genetab =  read.table(genefile, sep="\t", header=TRUE)
    geneList = genetab$x
  }
  return(geneList)
}

context = 'original'
original_genes = load_genelist('original', 'Quiescence')
original_genes = original_genes[original_genes %in% rownames(df)]

# use only for significance
original_genes = original_genes[sapply(original_genes, function(x){
  res = wilcox.test(as.numeric(responder[x,]), as.numeric(non_responder[x,]))
  if(res$p.value<0.05){
    return(TRUE)
  }else{
    return(FALSE)
  }
})]

context = 'context'
context_genes = load_genelist('SKCM', 'Quiescence')
context_genes = context_genes[context_genes %in% rownames(df)]

# use only for significance
context_genes = context_genes[sapply(context_genes, function(x){
  res = wilcox.test(as.numeric(responder[x,]), as.numeric(non_responder[x,]))
  if(res$p.value<0.05){
    return(TRUE)
  }else{
    return(FALSE)
  }
})]

pallete = brewer.pal(n=8, name='Dark2')
labs = c(rep(pallete[1], ncol(responder)), rep(pallete[3], ncol(non_responder)))


subdf = df2[original_genes,]
gplots::heatmap.2(t(as.matrix(subdf)), trace="none", density.info="none", Rowv=T,Colv=F,
                  #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
                  lhei=c(1.5,4),
                  RowSideColors=labs, col=brewer.pal(11,"RdBu"))
png(filename='/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_8b_effectsize_agnostic.png', 
    width= 12000, height = 10000, units = "px", pointsize=18, bg="white",  res=700)
gplots::heatmap.2(t(as.matrix(subdf)), trace="none", density.info="none", Rowv=T,Colv=F,
                  #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
                  keysize=1, key.xlab='Expression', key.title='', margins=c(8,8),
                  lhei=c(1.5,4), cexCol=1.5, labRow = "",
                  RowSideColors=labs, col=brewer.pal(11,"RdBu"))
dev.off()

subdf = df2[context_genes,]
gplots::heatmap.2(t(as.matrix(subdf)), trace="none", density.info="none", Rowv=T,Colv=F,
                  #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
                  lhei=c(1.5, 4 ), #cexCol=0.5,
                  RowSideColors=labs, col=brewer.pal(11,"RdBu"))

subdf = df2[context_genes,]
png(filename='/Users/sharminm2/Dropbox/nih/context/plots/hPPIN_DE2/Fig_8b_effectsize_context.png', 
    width= 12000, height = 10000, units = "px", pointsize=18, bg="white",  res=700)
gplots::heatmap.2(t(as.matrix(subdf)), trace="none", density.info="none", Rowv=T,Colv=F,
                  #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
                  keysize=1, key.xlab='Expression', key.title='', margins=c(8,8),
                  lhei=c(1.5,4), cexCol=1.5, labRow = "",
                  RowSideColors=labs, col=brewer.pal(11,"RdBu"))
dev.off()
