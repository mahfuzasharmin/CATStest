import numpy as np
import pandas as pd
from collections import Counter
from tqdm import tqdm
import pickle
import glob
import os

from itertools import chain
from scipy.stats import ranksums
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr, pearsonr
from scipy import stats

import multiprocessing as mp

from sklearn.preprocessing import scale
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

import seaborn as sns
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf as pdfbackend

from annot import *
from utils import read_genes

def compute_jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

def generate_jaccard_matrix(signature, annotationdf):
    
    js = {}
    for subType_1 in CANCER_TYPES:
        js[subType_1] = []
        for subType_2 in CANCER_TYPES:
            # if subType_1==subType_2:
            #     continue
            # else:
            gset_1 = read_genes("context", signature, subType_1, annotationdf)
            gset_2 = read_genes("context", signature, subType_2, annotationdf)
            js[subType_1].append(compute_jaccard(gset_1['ENSG'], gset_2['ENSG']))

    df = pd.DataFrame(js, index = BREAST_CANCER_SUBTYPES)
    
    return df

def prepare_depmap_scores(annotationdf, ensmbl_df, subType_option='.'):
    depdf =  pd.read_csv(f"{PROJECT_location}/depmap/CRISPR_gene_dependency.csv")
    depdf.set_index('DepMap_ID', inplace=True)
    entrez_ids  = [i.split(')')[0].split('(')[1] for i in depdf.columns]
    depdf.rename(columns=dict(zip(depdf.columns, entrez_ids)), inplace=True)
    
    metadata = pd.read_csv(f"{PROJECT_location}/depmap/metadata.csv", sep="\t")
    tcga_to_depmap_cancerType = TCGA_TO_DepMap_CANCERTYPE
    print(tcga_to_depmap_cancerType)
    
    cancerType_to_tissue = CANCERTYPE_TO_TISSUE
    cancerTypes = CANCER_TYPES
    subtype_info = None
    
    signatures = SIGNATURES 
    
    depmap_ids = np.unique(metadata['depmap_id'].values)
    for cancerType in cancerTypes:
        if cancerType not in tcga_to_depmap_cancerType.keys():
            continue
        disease = tcga_to_depmap_cancerType[cancerType]
        primary_disease = list(disease.keys())[0]
        subtype_disease = disease[primary_disease]
        
        disease_depmap_ids = metadata[(metadata['primary_disease']==primary_disease) & 
                                  (metadata['subtype_disease'].isin(subtype_disease))]['depmap_id'].values
            
        disease_depdf = depdf[depdf.index.isin(disease_depmap_ids)]

        nondisease_depmap_ids = np.unique(metadata[~(metadata['primary_disease']==primary_disease)]['depmap_id'].values)
        nondisease_depmap_ids = np.random.choice(nondisease_depmap_ids, size=disease_depdf.shape[0]*3, replace=False)
        nondisease_depdf = depdf[depdf.index.isin(nondisease_depmap_ids)].head(disease_depdf.shape[0])

        random_depmap_ids = np.random.choice(depmap_ids, size=disease_depdf.shape[0]*3, replace=False)
        random_depdf = depdf[depdf.index.isin(random_depmap_ids)].head(disease_depdf.shape[0])

        print(primary_disease, subtype_disease, disease_depdf.shape[0], nondisease_depdf.shape[0], random_depdf.shape[0])
        for signature in signatures:
            # signature genes
            genetab = read_genes("original", signature, cancerType, annotationdf); original_df = genetab
            sign_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # context genes
            genetab = read_genes("context", signature, cancerType, annotationdf); pruned = genetab
            context_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # discarded genes
            genetab = original_df.loc[~original_df['ENSG'].isin(pruned['ENSG']),:]
            discard_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # retained genes
            genetab = original_df.loc[original_df['ENSG'].isin(pruned['ENSG']),:]
            retain_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # hpa genes
            if cancerType not in CANCERTYPE_TO_HPA_TISSUE.keys():
                hpa_entrez_ids = np.array([])
            else:
                genetab = read_genes("frost", signature, cancerType, annotationdf); 
                hpa_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # ppi genes
            # genetab = read_genes("ppi", signature, cancerType, annotationdf)
            # ppi_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # de genes
            tissueType = cancerType_to_tissue[cancerType]
            genetab  = pd.read_csv(f'{PROJECT_location}/data/DE genes/{tissueType}_de2.txt', sep="\t")
            genetab = annotationdf[annotationdf['Symbol'].isin(genetab['x'].values)][['ENSG', 'Symbol']].drop_duplicates()
            de_entrez_ids = ensmbl_df[ensmbl_df['ensembl_gene_id'].isin(genetab['ENSG'].values)]['entrezgene_id']

            # sign_entrez_ids = list(set(sign_entrez_ids).difference(set(context_entrez_ids)))
            # context_entrez_ids = list(set(context_entrez_ids).difference(set(sign_entrez_ids)))
            # de_entrez_ids = list(set(de_entrez_ids).difference(set(sign_entrez_ids)))
            # ppi_entrez_ids = list(set(ppi_entrez_ids).difference(set(sign_entrez_ids)))
            # ppi_entrez_ids2 = list(set(ppi_entrez_ids).difference(set(context_entrez_ids)))

            if len(sign_entrez_ids)==0:
                print(cancerType, signature, ' -  no leftover signature genes')

            if len(context_entrez_ids)==0:
                print(cancerType, signature, ' -  no leftover context genes')

            if len(discard_entrez_ids)==0:
                print(cancerType, signature, ' -  no leftover discard genes')

            if len(retain_entrez_ids)==0:
                print(cancerType, signature, ' -  no leftover retain genes')
            
            if len(hpa_entrez_ids)==0:
                print(cancerType, signature, ' -  no leftover hpa genes')

            if len(de_entrez_ids)==0:
                print(cancerType, signature, ' -  no leftover de genes')

            # extract scores for signature genes
            signdf = disease_depdf.loc[:, disease_depdf.columns.isin(sign_entrez_ids)]
            signdf = pd.melt(signdf.reset_index(), id_vars=['DepMap_ID'], value_vars=signdf.columns)
            signdf['context'] = 'Signatures in relevant Cells'

            # extract scores for context genes
            contextdf = disease_depdf.loc[:, disease_depdf.columns.isin(context_entrez_ids)]
            contextdf = pd.melt(contextdf.reset_index(), id_vars=['DepMap_ID'], value_vars=contextdf.columns)
            contextdf['context'] = 'Context in relevant Cells'

            # extract scores for discarded genes
            discarddf = disease_depdf.loc[:, disease_depdf.columns.isin(discard_entrez_ids)]
            discarddf = pd.melt(discarddf.reset_index(), id_vars=['DepMap_ID'], value_vars=discarddf.columns)
            discarddf['context'] = 'Discard in relevant Cells'

             # extract scores for retained genes
            if len(retain_entrez_ids)==0:
                retaindf = None
            else:
                retaindf = disease_depdf.loc[:, disease_depdf.columns.isin(retain_entrez_ids)]
                retaindf = pd.melt(retaindf.reset_index(), id_vars=['DepMap_ID'], value_vars=retaindf.columns)
                retaindf['context'] = 'Retain in relevant Cells'

             # extract scores for hpa/frost genes
            if len(hpa_entrez_ids)==0 or cancerType not in CANCERTYPE_TO_HPA_TISSUE.keys():
                hpadf = None
            else:
                hpadf = disease_depdf.loc[:, disease_depdf.columns.isin(hpa_entrez_ids)]
                hpadf = pd.melt(hpadf.reset_index(), id_vars=['DepMap_ID'], value_vars=hpadf.columns)
                hpadf['context'] = 'HPA in relevant Cells'

            # extract scores for de genes
            dedf = disease_depdf.loc[:, disease_depdf.columns.isin(de_entrez_ids)]
            dedf = pd.melt(dedf.reset_index(), id_vars=['DepMap_ID'], value_vars=dedf.columns)
            dedf['context'] = 'DE in relevant Cells'

#             # extract scores for random genes from random tissues
#             signdf_alt = nondisease_depdf.loc[:, nondisease_depdf.columns.isin(sign_entrez_ids)]
#             signdf_alt = pd.melt(signdf_alt.reset_index(), id_vars=['DepMap_ID'], value_vars=signdf_alt.columns)
#             signdf_alt['context'] = 'Signatures in Others'

#             contextdf_alt = nondisease_depdf.loc[:, nondisease_depdf.columns.isin(context_entrez_ids)]
#             contextdf_alt = pd.melt(contextdf_alt.reset_index(), id_vars=['DepMap_ID'], value_vars=contextdf_alt.columns)
#             contextdf_alt['context'] = 'Context in Others'

#             random_entrez_ids = np.random.choice(entrez_ids, size=len(context_entrez_ids), replace=False)
#             randomdf = random_depdf.loc[:, depdf.columns.isin(random_entrez_ids)]
#             randomdf = pd.melt(randomdf.reset_index(), id_vars=['DepMap_ID'], value_vars=randomdf.columns)
#             randomdf['context'] = 'Random scores in random cells'

#             randomdf2 = disease_depdf.loc[:, depdf.columns.isin(random_entrez_ids)]
#             randomdf2 = pd.melt(randomdf2.reset_index(), id_vars=['DepMap_ID'], value_vars=randomdf2.columns)
#             randomdf2['context'] = 'Random scores in relevant Cells'

            scoredf = pd.concat([signdf, contextdf, dedf, discarddf, retaindf, hpadf])
            #print(scoredf.shape)
            sfile = f"{PROJECT_location}/depmap/{cancerType}3.xlsx"
            mode = 'a' if os.path.exists(sfile) else 'w'
            mode2 = 'replace' if os.path.exists(sfile) else None

            with pd.ExcelWriter(sfile, mode=mode, if_sheet_exists=mode2) as writer:  
                scoredf.to_excel(writer, sheet_name=signature, index=False)

    return None

def count_depmap_celllines(subType_option='.'):
    depdf =  pd.read_csv(f"{PROJECT_location}/depmap/CRISPR_gene_dependency.csv")
    depdf.set_index('DepMap_ID', inplace=True)
    entrez_ids  = [i.split(')')[0].split('(')[1] for i in depdf.columns]
    depdf.rename(columns=dict(zip(depdf.columns, entrez_ids)), inplace=True)
    
    metadata = pd.read_csv(f"{PROJECT_location}/depmap/metadata.csv", sep="\t")
    tcga_to_depmap_cancerType = TCGA_TO_DepMap_CANCERTYPE
    print(tcga_to_depmap_cancerType)
    
    cancerType_to_tissue = CANCERTYPE_TO_TISSUE
    if subType_option=='BRCA':
        cancerTypes = BREAST_CANCER_SUBTYPES
        subtype_df = pd.read_excel(f"{DB_location}/depmap/ccle_annotation.xlsx", sheet_name='subtypes')
        subtype_to_annot = {'BRCA_Basal': 'Basal', 'BRCA_LumA': 'Luminal A', 'BRCA_LumB': 'Luminal B', 'BRCA_Her2': 'Her2'}
    if subType_option=='Brain':
        cancerTypes = BRAIN_CANCER_SUBTYPES
        
        cancerTypes = ['Brain']
    else:
        cancerTypes = CANCER_TYPES
        subtype_info = None
    
    signatures = SIGNATURES #['SenMayo'] #
    
    celllines = set()
    for cancerType in cancerTypes:
        if cancerType not in tcga_to_depmap_cancerType.keys():
            continue
        disease = tcga_to_depmap_cancerType[cancerType]
        primary_disease = list(disease.keys())[0]
        subtype_disease = disease[primary_disease]
        
        if subType_option=='BRCA':
            disease_celllines = subtype_df[subtype_df['Annotation_My']==subtype_to_annot[cancerType]]['DepMap_ID'].values
        else:
            disease_celllines = metadata[(metadata['primary_disease']==primary_disease) & 
                                  (metadata['subtype_disease'].isin(subtype_disease))]['cell_line_name']
            print(disease_celllines.values)
            
        celllines = celllines.union(set(disease_celllines))
        
    return len(disease_celllines)



def count_cancer_risk_loci():
    cancerType_to_gwas = CANCERTYPE_to_GWAS
    
    gwas_df = pd.read_excel(f"{PROJECT_location}/cancer_risk_loci/41568_2017_BFnrc201782_MOESM22_ESM.xlsx",
                               sheet_name='Cancer Risk loci (March 2017)')
    
    cancerNames = cancerType_to_gwas.values()
    gwas_df = gwas_df[(gwas_df['Cancer type'].isin(cancerNames))]
    
    return gwas_df.shape[0]

def determine_cancer_risk_loci(cancerType, annotationdf):
    signatures = SIGNATURES
    
    cancerType_to_gwas = CANCERTYPE_to_GWAS
    
    gwas_df = pd.read_excel(f"{PROJECT_location}/cancer_risk_loci/41568_2017_BFnrc201782_MOESM22_ESM.xlsx",
                               sheet_name='Cancer Risk loci (March 2017)')
    
    cancerName = cancerType_to_gwas[cancerType]
    gwas_df = gwas_df[(gwas_df['Cancer type']==cancerName)]
    print(cancerType, cancerName)
    #print(gwas_df.head())
    gwas_df = gwas_df[(~gwas_df['Reported Genes'].isna()) & (gwas_df['Reported Genes']!='intergenic')]
    gwas_df = gwas_df[(gwas_df['Reported Genes'].apply(lambda x: isinstance(x, str)))]
    
    def check_gene(gene):
        genes = gene.split(',')
        genes = [i.strip() for i in genes]
        return genes
    
    genes = gwas_df['Reported Genes'].apply(lambda x: check_gene(x))
    genes = list(set(list(chain(*genes))))
    gene_df = annotationdf[annotationdf['Symbol'].isin(genes)]
    
    dfs = []
    for signature in signatures:
        genetab  = read_genes("original", signature, cancerType, annotationdf); originaldf = genetab
        original_ids = genetab['ENSG'].values
        original_count = gene_df[gene_df['ENSG'].isin(original_ids)].shape[0]

        genetab  = read_genes("context", signature, cancerType, annotationdf); contextdf = genetab
        context_ids = genetab['ENSG'].values
        context_count = gene_df[gene_df['ENSG'].isin(context_ids)].shape[0]

        genetab  = read_genes("de", signature, cancerType, annotationdf); dedf = genetab
        de_ids = genetab['ENSG'].values
        de_count = gene_df[gene_df['ENSG'].isin(de_ids)].shape[0]

        genetab  = originaldf.loc[~originaldf['ENSG'].isin(contextdf['ENSG']),:]
        discard_ids = genetab['ENSG'].values
        discard_count = gene_df[gene_df['ENSG'].isin(discard_ids)].shape[0]

        genetab  = originaldf.loc[originaldf['ENSG'].isin(contextdf['ENSG']),:]
        retain_ids = genetab['ENSG'].values
        retain_count = gene_df[gene_df['ENSG'].isin(retain_ids)].shape[0]
        retain_percentage = 0 if len(retain_ids)==0 else retain_count/len(retain_ids)
            
        if cancerType not in CANCERTYPE_TO_HPA_TISSUE.keys():
            hpa_percentage = 0
        else:
            hpadf  = read_genes("frost", signature, cancerType, annotationdf); 
            hpa_ids = hpadf['ENSG'].values
            hpa_count = gene_df[gene_df['ENSG'].isin(hpa_ids)].shape[0]
            hpa_percentage = 0 if len(hpa_ids)==0 else hpa_count/len(hpa_ids)
        
        genetab  = read_genes("ppi", signature, cancerType, annotationdf)
        ppi_ids = genetab['ENSG'].values
        ppi_count = gene_df[gene_df['ENSG'].isin(ppi_ids)].shape[0]
        
        # genetab  = pd.read_csv(f'/Users/msharmin/Dropbox/nih/context/data/TCGA/PPI_expanded/{signature}_Sig_PPI_Expanded3', sep="\t")
        # genetab = annotationdf[annotationdf['Symbol'].isin(genetab['x'].values)][['ENSG', 'Symbol']].drop_duplicates()
        # ppi_ids = genetab['ENSG'].values
        # ppi_count = gene_df[gene_df['ENSG'].isin(ppi_ids)].shape[0]
        
        dfs.append(pd.DataFrame({'signature': signature,
                                 'original_size': [len(original_ids)], 
                                 'ppi_size': [len(ppi_ids)],
                                 'context_size': [len(context_ids)],
                                 'discard_size': [len(discard_ids)],
                                 
                                 'original_rest': [len(original_ids)-original_count], 
                                 'ppi_rest': [len(ppi_ids)-ppi_count],
                                 'context_rest': [len(context_ids)-context_count],
                                 'discard_rest': [len(discard_ids)-discard_count],
                                 
                                 'original_count': [original_count], 
                                 'ppi_count': [ppi_count],
                                 'context_count': [context_count],
                                 'discard_count': [discard_count],
                                 
                                 'original_percentage': [original_count/len(original_ids)], 
                                 'ppi_percentage': [ppi_count/len(ppi_ids)],
                                 'context_percentage': [context_count/len(context_ids)],
                                 'de_percentage': [de_count/len(de_ids)],
                                 'discard_percentage': [discard_count/len(discard_ids)],
                                 'retain_percentage': [retain_percentage],
                                 'hpa_percentage': [hpa_percentage],
                                })
                  )

    df = pd.concat(dfs)
    return df


def run_gwas_overlaps(annotationdf):
    cancerTypes = list(CANCERTYPE_to_GWAS.keys())
    
    sfile = f"{PROJECT_location}/cancer_risk_loci/cancer_risk_summary.xlsx"
    for cancerType in cancerTypes:
        df = determine_cancer_risk_loci(cancerType, annotationdf)
        
        mode = 'a' if os.path.exists(sfile) else 'w'
        mode2 = 'replace' if os.path.exists(sfile) else None

        with pd.ExcelWriter(sfile, mode=mode, if_sheet_exists=mode2) as writer:  
            df.to_excel(writer, sheet_name=cancerType, index=False)
    
    return None

def load_expression_data(cancerType, annotationdf):
    expmat = pd.read_csv(f"{DB_location}/expression/TCGA-{cancerType}/matrix.csv.gz", 
                         header=0, index_col=0)
    genes = np.array([i.split(".")[0] for i in expmat.index])
    expmat.index = genes
    expmat = expmat[expmat.index.isin(annotationdf['ENSG'])]
    return expmat

def extend_clinical_info_with_subtypes():
    clinical = pd.read_csv(f"{DB_location}/tcga-clinical/TCGAclinicalDataFile-PFS.txt", sep='\t')
    subtypes = pd.read_csv(f"{DB_location}/tcga-clinical/BRCA.547.PAM50.SigClust.Subtypes.txt", sep='\t')
    subtypes['shortname'] = subtypes['Sample'].apply(lambda x: '-'.join(x.split('-')[:3]))
    extendeddf = clinical.merge(subtypes, left_on='samples', right_on='shortname')
    extendeddf = extendeddf.drop(columns=['Sample', 'Type', 'Siglust', 'shortname'])
    extendeddf = extendeddf.rename(columns={'PAM50': 'subtype'}) #didnot work
    extendeddf['subtype'] = 'BRCA_' + extendeddf['subtype']
    extendeddf
    extendeddf.to_csv(f"{PROJECT_location}/tcga-clinical/TCGAclinicalDataFileSubtypeExtended-PFS.txt", sep='\t')
    return None

def conduct_expression_survival_analysis(cancerType, genetab, expmat, gene_column_name):
    
    subTypes = BREAST_CANCER_SUBTYPES
    if cancerType in subTypes:
        clinical = pd.read_csv(f"{PROJECT_location}/tcga-clinical/TCGAclinicalDataFileSubtypeExtended-PFS.txt", sep='\t')
        clinical  = clinical[clinical['subtype']==cancerType]
    else:
        clinical = pd.read_csv(f"{DB_location}/tcga-clinical/TCGAclinicalDataFile-PFS.txt", sep='\t')
        if cancerType=='Brain':
            clinical  = clinical[clinical['type'].isin(['LGG', 'GBM'])]
        else:
            clinical  = clinical[clinical['type']==cancerType]
    
    clinical.loc[:,"race"] = clinical["race"].fillna(clinical["race"].mode()[0])
    #print(clinical.shape)
    survdf = []
    for k,idx  in enumerate(tqdm(genetab.index)):
        gene = genetab.loc[idx,][gene_column_name]
        #print(gene)
        if gene not in expmat.index:
            print(gene)
            continue
        
        gene_expmat = pd.DataFrame(expmat.loc[gene,])
        
        samples = ['-'.join(i.split('-')[0:3]) for i in gene_expmat.index]
        gene_expmat.index = samples
        clinical = clinical[clinical['samples'].isin(samples)]
        #print(clinical.shape)
        data = gene_expmat.reset_index().merge(clinical[clinical['samples'].isin(samples)], left_on='index', right_on='samples')
        #print(data.shape)
        data = data.set_index('index')
        data = data.replace({
            'sex': {key:val for val, key in enumerate(np.sort(np.unique(data['sex'])))},
            'race': {key:val for val, key in enumerate(np.sort(np.unique(data['race'])))}
        })
        
        cols_to_drop = ['age', 'race', 'sex', 'samples', 'type', 'stage', 'subtype', 'OS.time', 'OS.status']
        for covariate in cols_to_drop:
            if covariate not in data.columns:
                continue
            if covariate in ['age']:
                continue
            else:
                data = data.drop(columns=[covariate])
                
        data.rename(columns = {gene: 'score'}, inplace = True)
        data = data.apply(lambda  x: x.fillna(np.nanmedian(x)),  axis=0)
        
        data.loc[:,'age'] = scale(data['age'])
        #feature_columns = list(data.columns.difference(['time', 'status']))
        #sdata = scale(data[feature_columns])
        #sdata = pd.DataFrame(sdata, index=data.index, columns=data.columns[:len(feature_columns)])
        #data = sdata.join(data[['time', 'status']])
        
        if (data.var()<1e-4).any():
            continue
        
        cph = CoxPHFitter(penalizer=0.001)
        cph.fit(data, duration_col = 'PFI.time', event_col = 'PFI.status') #
        hrs = cph.hazard_ratios_.to_dict()
        hrs['pvalue'] = cph.summary.loc['score', 'p']
        hrs['ENSG'] = gene
        survdf.append(hrs)
        
    df = pd.DataFrame(survdf)
    df2 = df.merge(genetab, left_on='ENSG', right_on=gene_column_name)
    
    return df2

def record_background_survivals(annotationdf):
    sfile = f"{PROJECT_location}/survival/background_survivals.PFI.xlsx"
    
    cancerTypes = BREAST_CANCER_SUBTYPES #CANCER_TYPES #['Brain'] #
    # cancerTypes  = ['THCA']
    for cancerType  in cancerTypes:
        print(cancerType)
        if cancerType=='Brain':
            expmat1 = load_expression_data('GBM', annotationdf)
            expmat2 = load_expression_data('LGG', annotationdf)
            expmat =  pd.concat([expmat1, expmat2], ignore_index=False, axis=1)
        else:
            expmat = load_expression_data(cancerType, annotationdf)
        genetab = pd.DataFrame(expmat.index)
        df = conduct_expression_survival_analysis(cancerType, genetab, expmat, 
                                                  gene_column_name=0)
        
        mode = 'a' if os.path.exists(sfile) else 'w'
        mode2 = 'replace' if os.path.exists(sfile) else None
        
        with pd.ExcelWriter(sfile, mode=mode, if_sheet_exists=mode2) as writer:
            df.to_excel(writer, sheet_name=cancerType, index=False)

    return None

def run_ssgsea_commands(tissue_option='cancerType'):
    cancerTypes = CANCER_TYPES
    
    cancerType_to_tissue = CANCERTYPE_TO_TISSUE
    
    # cancerTypes = ['BRCA_Her2']
    
    for cancerType in cancerTypes:
        if cancerType!='ACC':
            continue
        tissueType = cancerType_to_tissue[cancerType]
        
        cmd  = f'Rscript {PROJECT_location}/R/run_gsea.R {cancerType} original'
        print(cmd)
        # os.system(cmd)
        
        cmd  = f'Rscript {PROJECT_location}/R/run_gsea.R {cancerType} ppi'
        print(cmd)
        os.system(cmd)
        
        cmd  = f'Rscript {PROJECT_location}/R/run_gsea.R {cancerType} context'
        print(cmd)
        # os.system(cmd)
        
        cmd  = f'Rscript {PROJECT_location}/R/run_gsea.R {cancerType} de {tissueType}'
        print(cmd)
        # os.system(cmd)
        #break
        
    return None


def conduct_joint_ssgsea_survival_analysis(cancerType, cancerType_2=None, tissueType=None):
    
    if cancerType_2 is None:
        if tissueType is None:
            sfile = f"{PROJECT_location}/data/results/common/ssgsea_scores/original.xlsx"
        else:
            sfile = f"{PROJECT_location}/data/results/common/ssgsea_scores/de.xlsx"
        scores_signature = pd.read_excel(sfile, sheet_name=cancerType, index_col=0)
    else:
        sfile = f"{PROJECT_location}/data/results/common/ssgsea_scores/original_{cancerType}.xlsx"
        scores_signature = pd.read_excel(sfile, sheet_name=cancerType_2, index_col=0)
    samples = ['-'.join(i.split('.')[0:3]) for i in scores_signature.index]
    scores_signature.index = samples

    if cancerType_2 is None:
        if tissueType is None:
            sfile = f"{PROJECT_location}/data/results/{CUR_DST}/ssgsea_scores/context.xlsx"
        else:
            sfile = f"{PROJECT_location}/data/results/{CUR_DST}/ssgsea_scores/context.xlsx"
        scores_context = pd.read_excel(sfile, sheet_name=cancerType, index_col=0)
    else:
        sfile = f"{PROJECT_location}/data/results/{CUR_DST}/ssgsea_scores/context_{cancerType}.xlsx"
        scores_context = pd.read_excel(sfile, sheet_name=cancerType_2, index_col=0)
    samples = ['-'.join(i.split('.')[0:3]) for i in scores_context.index]
    scores_context.index = samples
    
    # sfile = f"/Users/msharmin/Dropbox/nih/context/data/ssgsea_scores/ppi.xlsx"
    # scores_ppi = pd.read_excel(sfile, sheet_name=cancerType, index_col=0)
    # samples = ['-'.join(i.split('.')[0:3]) for i in scores_ppi.index]
    # scores_ppi.index = samples

    subTypes = BREAST_CANCER_SUBTYPES
    if cancerType in subTypes:
        clinical = pd.read_csv(f"{PROJECT_location}/tcga-clinical/TCGAclinicalDataFileSubtypeExtended-PFS.txt")
        if cancerType_2 is None:
            clinical  = clinical[clinical['subtype']==cancerType]
        else:
            clinical = clinical[clinical['subtype']==cancerType_2]
    else:
        clinical = pd.read_csv(f"{PROJECT_location}/tcga-clinical/TCGAclinicalDataFile-PFS.txt")
        clinical  = clinical[clinical['type']==cancerType]
    
    # print(clinical.shape)
    clinical.loc[:,"race"] = clinical["race"].fillna(clinical["race"].mode()[0])
    clinical = clinical[clinical['samples'].isin(samples)]
    # print(clinical.shape)
    survdf = []
    signatures = list(set(scores_signature.columns).intersection(set(SIGNATURES)))
    for signature in signatures:
        # print(signature)
        score_mat = pd.DataFrame({'scores_original': stats.zscore(scores_signature.loc[:,signature]),
                                  #'scores_ppi': scores_ppi.loc[:,signature],
                                  'scores_context': stats.zscore(scores_context.loc[:,signature]),
                                 })
        data = score_mat.reset_index().merge(clinical, left_on='index', right_on='samples')
        data = data.set_index('index')
        data = data.replace({
            'sex': {key:val for val, key in enumerate(np.sort(np.unique(data['sex'])))},
            'race': {key:val for val, key in enumerate(np.sort(np.unique(data['race'])))}
        })
        
        cols_to_drop = ['age', 'race', 'sex', 'samples', 'type', 'stage', 'subtype']
        for covariate in cols_to_drop:
            if covariate not in data.columns:
                continue
            if covariate in ['age']:
                continue
            else:
                data = data.drop(columns=[covariate])
            
        data = data.apply(lambda  x: x.fillna(np.nanmedian(x)),  axis=0)
        data.loc[:,'age'] = scale(data['age'])
        
        #feature_columns = list(data.columns.difference(['time', 'status']))
        #sdata = scale(data[feature_columns])
        #sdata = pd.DataFrame(sdata, index=data.index, columns=data.columns[:len(feature_columns)])
        #data = sdata.join(data[['time', 'status']])
        
        if (data.var()<1e-4).any():
            return None

        cph = CoxPHFitter(penalizer=0.001)
        cph.fit(data, duration_col = 'time', event_col = 'status')
        hrs = cph.hazard_ratios_.to_dict()
        hrs['pvalue_original'] = cph.summary.loc['scores_original', 'p']
        hrs['pvalue_context'] = cph.summary.loc['scores_context', 'p']
        hrs['signature'] = signature
        survdf.append(hrs)
    
    df = pd.DataFrame(survdf)
    
    return df

def record_joint_ssgsea_survivals(tissue_option='cancerType'):
    #
    cancerTypes = CANCER_TYPES
    
    cancerType_to_tissue = CANCERTYPE_TO_TISSUE
    tissueType = None
    for cancerType  in cancerTypes:
        # print(cancerType)
        if tissue_option=='tissue' or tissue_option=='tissue2':
            tissueType = cancerType_to_tissue[cancerType]
            
        df = conduct_joint_ssgsea_survival_analysis(cancerType, tissueType=tissueType)
        df.set_index('signature', inplace=True)
        if tissue_option=='cancerType':
            sfile = f"{PROJECT_location}/survival/joint_zssgsea_summary_context_vs_original.xlsx"
        elif tissue_option=='tissue':
            print(cancerType, tissueType)
            sfile = f"{PROJECT_location}/survival/joint_zssgsea_summary_de_vs_original.xlsx"
        elif tissue_option=='tissue2':
            print(cancerType, tissueType)
            sfile = f"{PROJECT_location}/survival/joint_zssgsea_summary_context_vs_de.xlsx"
        mode = 'a' if os.path.exists(sfile) else 'w'
        mode2 = 'replace' if os.path.exists(sfile) else None

        with pd.ExcelWriter(sfile, mode=mode, if_sheet_exists=mode2) as writer:  
            df.to_excel(writer, sheet_name=cancerType, index=True)
    
    return None

    