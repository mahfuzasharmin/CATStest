library(glue)

PROJECT_location = 'CATSevalData'
source("{PROJECT_location}/FROST/TSGST.R")
setwd(glue('{PROJECT_location}/FROST/TissueSpecificWeightGeneration/'))
gene.names = c("TP53", "MYC")
gene.ids = c('ENSG00000141510', 'ENSG00000136997')
tissue.types = c("skin", "liver", "kidney")
weights = getTissueSpecificGeneWeights(gene.names=gene.names,
                                      tissue.types=tissue.types,
                                      evidence.type="RNA")


cancerType_to_hpa_tissue = list('ACC'='adrenal.gland', 'BRCA'='breast', 'CESC'='cervix', 'COAD'='colon', 
                            'ESCA'='esophagus', 'KICH'='kidney', 'KIRC'='kidney', 'KIRP'='kidney', 
                            'LIHC'='liver', 'LUAD'='lung',  'LUSC'='lung', 'OV'='ovary', 
                            'PAAD'='pancreas', 'PRAD'='prostate', 'READ'='rectum', 
                            'SKCM'='skin', 'STAD'='stomach', 'THCA'='thyroid.gland')
tissues = c('adrenal.gland', 'breast', 'cervix', 'colon', 'esophagus', 'kidney', 
            'liver', 'lung', 'ovary', 'pancreas', 'prostate', 'rectum', 'skin', 
            'stomach', 'thyroid.gland')

# common
gene.ids = unique(unlist(lapply(signatures, read_genes)))

weights = getTissueSpecificGeneWeights(gene.ids=gene.ids,
                                       tissue.types=tissues,
                                       evidence.type="RNA")
write.table(weights, file='computed_weights.tsv', quote = F, sep='\t')

df =read.table(glue('{PROJECT_location}/FROST/TissueSpecificWeightGeneration/combined_hpa_data.csv'), sep=',', header=T); head(df)
df =read.table(glue('{PROJECT_location}/FROST/TissueSpecificWeightGeneration/computed_weights.tsv'), sep='\t', header=T); head(df)


# not common - get the signature; run frost across tissues

for(signature in signatures){
  gene.ids = read_genes(signature) 
  weights = getTissueSpecificGeneWeights(gene.ids=gene.ids, tissue.types=tissues, evidence.type="RNA")
  write.table(weights, file=glue('{PROJECT_location}/FROST/TissueSpecificWeightGeneration/computed_weights_{signature}.tsv'), quote = F, sep='\t')
}


read_frost_genes <- function(signature, cancerType){
  tissueType = cancerType_to_hpa_tissue[[cancerType]]
  genetab  = read.table(glue('{PROJECT_location}/FROST/TissueSpecificWeightGeneration/computed_weights_{signature}.tsv'), sep="\t")
  
  genes = unique(rownames(genetab[genetab[[tissueType]] > 0.6927579, ])); 
  # genes = unique(rownames(genetab[order(genetab[[tissueType]], decreasing=T), ]))[1:nsize]
  genetab = unique(annotationdf[annotationdf$ENSG  %in% genes,][,c('ENSG', 'gene_name')])
  genetab$ENSG
}

osizes = list()
fsizes = list()
for(signature in signatures){
  oritab = read_genes(signature); m = length(oritab)
  for(cancerType in names(cancerType_to_hpa_tissue)){
      gene.ids = read_frost_genes(signature, cancerType) 
      n = length(gene.ids)
      #print(glue('{signature} - {cancerType}: {m} vs {n}'))
      osizes[[glue('{signature} - {cancerType}')]] = m
      fsizes[[glue('{signature} - {cancerType}')]] = n
  }
}
summary(unlist(osizes))
summary(unlist(fsizes))

values = list()
for(signature in signatures){
      tissueType = cancerType_to_hpa_tissue[[cancerType]]
      genetab  = read.table(glue('{PROJECT_location}/FROST/TissueSpecificWeightGeneration/computed_weights_{signature}.tsv'), sep="\t")
      values[[glue('{signature}-{cancerType}')]] = genetab
}
df = do.call(rbind, values); dim(df); #1574
quantile(unlist(df), seq(0, 1, 0.05))
mean(unlist(df))
k = 0.6927579

#genetab = unique(annotationdf[annotationdf$gene_name  %in% genes,][,c('ENSG', 'gene_name')])