import numpy as np
import pandas as pd
import biomart
from annot import *

def get_annotation():
    genedf = pd.read_csv(f"{PROJECT_location}/gencode/gencode.v26.GRCh38.genes.gtf", 
                         comment="#", header=None, sep="\t")
    pcgenedf = genedf[(genedf[2]=='gene') & (genedf[8].str.contains("protein_coding"))]
    gene_id = pcgenedf.apply(lambda x: x[8].split(";")[0].split()[1].strip('\"').split(".")[0], axis=1)
    gene_name = pcgenedf.apply(lambda x: x[8].split(";")[3].split()[1].strip('\"'), axis=1)
    annotationdf = pd.DataFrame({'ENSG': gene_id, 'Symbol': gene_name})
    pcgenedf.columns = ['chrom', 'convention', 'gene', 'start', 'end', 'opt2', 'strand', 'opt2', 'details']
    annotationdf = pd.concat([pcgenedf[['chrom', 'start', 'end', 'strand']], annotationdf], axis=1)
    annotationdf = annotationdf.reset_index().drop(columns=['index'])
    annotationdf = annotationdf[annotationdf['chrom']!="chrM"]

    return annotationdf

def convert_ensembl_to_entrez(ensembl_gene_ids):
    """
    Converts a list of Ensembl Gene IDs to HGNC symbols using Ensembl BioMart.

    Args:
        ensembl_gene_ids (list): A list of Ensembl Gene IDs (e.g., ['ENSG00000162367', 'ENSG00000187048']).

    Returns:
        dict: A dictionary mapping Ensembl Gene IDs to their corresponding HGNC symbols.
              Returns None for IDs that could not be mapped.
    """
    # Connect to the Ensembl BioMart server
    server = biomart.BiomartServer('http://www.ensembl.org/biomart')

    # Select the human genes dataset
    # 'hsapiens_gene_ensembl' is for Homo sapiens (human) genes
    ensembl_dataset = server.datasets['hsapiens_gene_ensembl']

    # Define the attributes to retrieve and the filters to apply
    # 'ensembl_gene_id' is the filter for input IDs
    # 'hgnc_symbol' is the attribute to retrieve
    ensembl_to_entrez = {}
    
    n = 100
    for i in range(0, len(ensembl_gene_ids), n):
        index = i
        chunk = ensembl_gene_ids[i:i + n]
        
        response = ensembl_dataset.search({'filters': {'ensembl_gene_id': chunk},
                                            'attributes': ['ensembl_gene_id', 'entrezgene_id'],
                                         })

        # Parse the response and create the mapping dictionary
        for line in response.iter_lines():
            line = line.decode('utf-8')
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) == 2:
                ensembl_id, hgnc_symbol = parts[0], parts[1]
                ensembl_to_entrez[ensembl_id] = hgnc_symbol if hgnc_symbol else None
        # else:
        #     # Handle cases where HGNC symbol might be missing for an Ensembl ID
        #     ensembl_id = parts[0]
        #     ensembl_to_entrez[ensembl_id] = None

    # Add any input IDs that were not found in the BioMart query
    # for ensembl_id in ensembl_gene_ids:
    #     if ensembl_id not in ensembl_to_entrez:
    #         ensembl_to_entrez[ensembl_id] = None

    return ensembl_to_entrez
    
def read_genes(data_option, signature, cancerType, annotationdf):
        
    if data_option=='original':
        #sfile = f"{PROJECT_location}/data/{survival_option}/{cancerType}_original.xlsx"
        genetab  = pd.read_csv(f'{PROJECT_location}/cancerSEA/{signature}.txt', sep="\t")
        
    elif data_option=='ppi':
        # sfile = f"{PROJECT_location}/data/{survival_option}/{cancerType}_ppi.xlsx"
        gfile = f'{PROJECT_location}/Step1-PPI_expanded/{signature}_Sig_PPI_Expanded3'
        # print(gfile)
        genetab  = pd.read_csv(gfile, sep="\t")
        genetab = annotationdf[annotationdf['Symbol'].isin(genetab['x'].values)][['ENSG', 'Symbol']].drop_duplicates()
        
    elif data_option=='coexp':
        # sfile = f"{PROJECT_location}/data/{survival_option}/{cancerType}_ppi.xlsx"
        gfile = f'{PROJECT_location}/Step2-PPI_CoExp_expanded/{cancerType}/{cancerType}_{signature}_Sig_PPI_CoExpr_Expanded3'
        # print(gfile)
        genetab  = pd.read_csv(gfile, sep="\t")
        genetab = annotationdf[annotationdf['Symbol'].isin(genetab['x'].values)][['ENSG', 'Symbol']].drop_duplicates()
        
    elif data_option=='context':
        #sfile = f"{PROJECT_location}/data/{survival_option}/{cancerType}_context.xlsx"
        gfile = f'{PROJECT_location}/Step3-Prune_DE/{cancerType}/{cancerType}_{signature}_DE2'
        genetab  = pd.read_csv(gfile, sep="\t")
        genetab = annotationdf[annotationdf['Symbol'].isin(genetab['x'].values)][['ENSG', 'Symbol']].drop_duplicates()
        
    elif data_option=='de':
        cancerType_to_tissue = CANCERTYPE_TO_TISSUE
        tissueType = cancerType_to_tissue[cancerType]
        genetab  = pd.read_csv(f'{PROJECT_location}/DE genes/{tissueType}_de2.txt', sep="\t")
        #print(cancerType, tissueType)
        #print(genetab.head())
        genetab = annotationdf[annotationdf['Symbol'].isin(genetab['x'].values)][['ENSG', 'Symbol']].drop_duplicates()
        pass

    elif data_option=='frost':
        #genetab  = pd.read_csv(f'{PROJECT_location}/data/cancersea/hallmarks/{signature}.txt', sep="\t")
        #nsize = genetab.shape[0]
        
        cancerType_to_tissue = CANCERTYPE_TO_HPA_TISSUE
        tissueType = cancerType_to_tissue[cancerType]
        genetab  = pd.read_table(f'{PROJECT_location}/FROST/TissueSpecificWeightGeneration/computed_weights_{signature}.tsv', sep="\t")
        #genes = np.array(genetab[tissueType].sort_values(ascending=False)[:nsize].index.unique().tolist())
        # print(f'{PROJECT_location}/data/HPA/TissueSpecificWeightGeneration/computed_weights_{signature}.tsv')
        genes = np.array(genetab.loc[genetab[tissueType] > 0.6927579, :].index.unique().tolist())
        
        genetab = annotationdf[annotationdf['ENSG'].isin(genes)][['ENSG', 'Symbol']].drop_duplicates()
        pass
    
    return genetab



