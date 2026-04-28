import numpy as np
import pandas as pd

#import upsetplot
#import pyupset as pyu

import seaborn as sns
import matplotlib.cm as cm
from matplotlib import pyplot as plt
# from matplotlib import rc, rcParams
# import patchworklib as pw
import matplotlib.backends.backend_pdf as pdfbackend

from scipy.stats import spearmanr, pearsonr
from scipy.stats import mannwhitneyu, wilcoxon

from utils import read_genes
from annot import *


def aggregate_FC_data(annotationdf, subType_option='.'):
    signatures = SIGNATURES
    
    cancerTypes = CANCER_TYPES
    if subType_option=='BRCA':
        cancerTypes = BREAST_CANCER_SUBTYPES
    dfs = []
    for signature in signatures:
        df = []
        for cancerType in cancerTypes:
            original  = read_genes("original", signature, cancerType, annotationdf)
            ppi  = read_genes("ppi", signature, cancerType, annotationdf)
            coexp  = read_genes("coexp", signature, cancerType, annotationdf)
            pruned  = read_genes("context", signature, cancerType, annotationdf)
            n = original.shape[0]
            df.append({'Reference': original.shape[0]/n,#/n, 
                       'PPI Expanded': ppi.shape[0]/n,#/n,
                       'CoExp Expanded': coexp.shape[0]/n,#/n,
                       'DEG Pruned': pruned.shape[0]/n,#/n,
                       'signature': signature,
                       'cancerType': cancerType
                      })
        
        df = pd.DataFrame(df)
        df = pd.melt(pd.DataFrame(df), 
                     id_vars=['cancerType', 'signature'], 
                     value_vars=['Reference', 'PPI Expanded', 'CoExp Expanded', 'DEG Pruned']) # 
        dfs.append(df)
    df = pd.concat(dfs) 
    df['signature_variable'] = df['signature'] + '_' + df['variable']
    df['cancerType_variable'] = df['cancerType'] + '_' + df['variable']
    
    return df

def plot_FC_data_main(df, figno):
    fig, ax = plt.subplots(1, 1, figsize=(8,8), tight_layout=True, dpi=700)
    # for i, signature in enumerate(SIGNATURES[1:]):
        # ax = axes[int(i/5), i%5]
    sns.lineplot(ax=ax, data=df, lw=3,
                 x='variable', y='value', hue='cancerType')
    # ax.set_title(signature, fontsize=25)
    ax.set_ylabel('FC of geneset size', fontsize=25)
    ax.set_xlabel('', fontsize=25)
    ax.tick_params(axis='x', rotation=50, labelsize=15)
    ax.tick_params(axis='y', rotation=0, labelsize=15)
    # plt.setp(ax.get_legend().get_texts(), fontsize='13') # for legend text
    # plt.setp(ax.get_legend().get_title(), fontsize='15')

    print(f'{PROJECT_location}/Figures/Fig_{figno}.png')
    fig.savefig(f'{PROJECT_location}/Figures/Fig_{figno}.png')

    return ax

def plot_FC_data_supplementary(df, figno):
    fig, axes = plt.subplots(3, 5, figsize=(30,25), tight_layout=True, dpi=700)
    for i, signature in enumerate(SIGNATURES[1:]):
        ax = axes[int(i/5), i%5]
        sns.lineplot(ax=ax, data=df[df['signature']==signature], 
                     x='variable', y='value', hue='cancerType')
        ax.set_title(signature, fontsize=25)
        ax.set_ylabel('Fraction', fontsize=25)
        ax.set_xlabel('Step', fontsize=25)
        
        plt.setp(ax.get_legend().get_texts(), fontsize='13') # for legend text
        plt.setp(ax.get_legend().get_title(), fontsize='15')

    print(f'{PROJECT_location}/Figures/Fig_{figno}.png')
    fig.savefig(f'{PROJECT_location}/Figures/Fig_{figno}.png')

    return ax


def aggregate_depmap_figdata(survival_option='survival-age', subType_option='.', panel_option='top'):
    tcga_to_depmap_cancerType = TCGA_TO_DepMap_CANCERTYPE
    
    signatures = SIGNATURES

    cancerTypes = CANCER_TYPES
    
    if subType_option=='BRCA':
        cancerTypes = ['BRCA'] + BREAST_CANCER_SUBTYPES
    if subType_option=='Brain':
        cancerTypes = ['Brain'] + BRAIN_CANCER_SUBTYPES
    pruned_dfs = []
    original_dfs = []
    availableTypes = []
    for cancerType in cancerTypes:
        if cancerType not in tcga_to_depmap_cancerType:
            continue
        # print(cancerType)
        availableTypes.append(cancerType)
        pruned_records = {}
        original_records = {}
        for signature in signatures:
            sfile = f"{PROJECT_location}/data/results/{CUR_DST}/depmaps/{cancerType}3.xlsx"
            scoredf = pd.read_excel(sfile, sheet_name=signature)
            scoredf['signature'] = signature
            
            if panel_option=='top': #context vs reference
                original_df = scoredf[(scoredf['context']=='Signatures in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = original_df[(original_df['value']>0.5)]
                original_fraction = sub_df.shape[0]/original_df.shape[0]

                pruned_df = scoredf[(scoredf['context']=='Context in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = pruned_df[(pruned_df['value']>0.5)]
                pruned_fraction = sub_df.shape[0]/pruned_df.shape[0]
                
            elif panel_option=='de': # context vs DE
                original_df = scoredf[(scoredf['context']=='DE in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = original_df[(original_df['value']>0.5)]
                original_fraction = sub_df.shape[0]/original_df.shape[0]

                pruned_df = scoredf[(scoredf['context']=='Context in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = pruned_df[(pruned_df['value']>0.5)]
                pruned_fraction = sub_df.shape[0]/pruned_df.shape[0]
                
            elif panel_option=='hpa': # context vs hpa
                original_df = scoredf[(scoredf['context']=='HPA in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = original_df[(original_df['value']>0.5)]
                if original_df.shape[0]==0:
                    original_fraction = 0
                else:
                    original_fraction = sub_df.shape[0]/pruned_df.shape[0]
            
                pruned_df = scoredf[(scoredf['context']=='Context in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = pruned_df[(pruned_df['value']>0.5)]
                pruned_fraction = sub_df.shape[0]/pruned_df.shape[0]
                
            else: # retain vs discard
                original_df = scoredf[(scoredf['context']=='Discard in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = original_df[(original_df['value']>0.5)]
                original_fraction = sub_df.shape[0]/original_df.shape[0]
            
                pruned_df = scoredf[(scoredf['context']=='Retain in relevant Cells') & (scoredf['signature']==signature)]
                sub_df = pruned_df[(pruned_df['value']>0.5)]
                if pruned_df.shape[0]==0:
                    pruned_fraction = 0
                else:
                    pruned_fraction = sub_df.shape[0]/pruned_df.shape[0]
            
            pruned_records[signature] = pruned_fraction
            original_records[signature] = original_fraction
            
        pruned_dfs.append(pruned_records)
        original_dfs.append(original_records)
    
    pruned_df = pd.DataFrame(pruned_dfs, index=availableTypes).melt()
    pruned_df.rename(columns = {'value':'context-specific', 'variable':'hallmark'}, inplace = True)
    original_df = pd.DataFrame(original_dfs, index=availableTypes).melt()
    original_df.rename(columns = {'value':'context-agnostic'}, inplace = True)
    df = pd.concat([pruned_df, original_df], axis=1)
    
    # from scipy.stats import mannwhitneyu, wilcoxon
    print(wilcoxon(pruned_df['context-specific'].astype('float'), original_df['context-agnostic'].astype('float'), alternative='greater'))
    return df

def plot_risk_loci_overlaps(yColumn, yPercentage, xColumn, xPercentage, figno):
    # yColumn = 'Context' #
    # yPercentage = 'context' #
    # xColumn = 'DE' #Context, Reference, Frost-refined, Reference-retained, Reference-discarded
    # xPercentage = 'de' #context, original, hpa, retain, discard
    rc('font', weight='bold')
    # rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

    sfile = f"{PROJECT_location}/cancer_risk_loci/gwas/cancer_risk_summary.xlsx"
    gwas_dfs = pd.read_excel(sfile, sheet_name=None)
    cancerTypes = gwas_dfs.keys()
    
    df = {}
    scatterdfs = []
    for cancerType in cancerTypes:
        gwas_df = gwas_dfs[cancerType]
        df[cancerType] = gwas_df[f'{yPercentage}_percentage'] - gwas_df[f'{xPercentage}_percentage']
        scatterdf = gwas_df[['signature', f'{yPercentage}_percentage', f'{xPercentage}_percentage']]
        scatterdf['cancerType'] = cancerType
        scatterdfs.append(scatterdf)
        
    scatterdf = pd.concat(scatterdfs)
    scatterdf.columns = ['Hallmark', yColumn, xColumn, 'cancerType']
        
    df = pd.DataFrame(df)
    df.index = gwas_df['signature']
    # df = df.drop(['Neuroendocrine'])
    print(np.sum((df>0).sum()), np.sum((df==0).sum()), np.sum((df<0).sum()), df.shape[0]*df.shape[1])
    print(np.mean((df[df>0]).max()), np.mean((df[df<0]).min()))
    
    cmap = sns.diverging_palette(300, 145, s=60, as_cmap=True)
    if xColumn=='Reference':
        fig, axes = plt.subplots(1, 2, figsize=(16,7), tight_layout=True) #
        ax = axes[1]
    else:
        fig, axes = plt.subplots(1, 1, figsize=(10,7), tight_layout=True) #
        ax = axes

    sns.heatmap(ax=ax, data=df, cmap=cmap, vmin=-0.1, vmax=0.1)  
    ax.tick_params(axis='y', rotation=0, size = 5, labelsize=15)
    ax.tick_params(axis='x', rotation=90, size = 5, labelsize=15)
    #ax.set_ylabel('cancerType', fontsize=25, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    # ax.set_title('Fraction of genes with cancer specific gwas overlap\ncontext - reference')

    if xColumn=='Reference':
        ax = axes[0]
        sns.scatterplot(ax=ax, data=scatterdf, x=xColumn, y=yColumn, hue='Hallmark')
        ax.set_title('Fraction of genes with cancer risk loci')
        ax.set_ylim((-0.01, max(scatterdf[yColumn])+0.01))
        ax.set_xlim((-0.01, max(scatterdf[yColumn])+0.01))
        ax.axline((0,0), slope=1, color='grey', linestyle='--')
        # ax.set_xlabel('')
        # ax.set_ylabel('')
    
    print(wilcoxon(scatterdf[yColumn].astype('float'), scatterdf[xColumn].astype('float'), alternative='greater'))
    
    fig.savefig(f'{PROJECT_location}/Figures/Fig_{figno}.png', dpi=700)
    return scatterdf


def aggregate_survival_summary_data(annotationdf, subType_option='.', panel_option='top'):
    signatures = SIGNATURES

    cancerTypes = CANCER_TYPES
    
    if subType_option=='BRCA':
        cancerTypes = ['BRCA'] + BREAST_CANCER_SUBTYPES
    elif subType_option=='Brain':
        cancerTypes = ['Brain'] + BRAIN_CANCER_SUBTYPES
    if panel_option=='hpa':
        cancerTypes = CANCERTYPE_TO_HPA_TISSUE.keys()
    
    dfs = []
    reference_dfs = []
    context_dfs = []
    for cancerType in cancerTypes:
        survdf = pd.read_excel(f"{PROJECT_location}/survival/background_survivals.PFI.xlsx", 
                               sheet_name=cancerType)
        records = {}
        reference_records = {}
        context_records = {}
        for signature in signatures:

            #conrext
            if panel_option=='top': #context vs original (green expected)
                pruned  = read_genes("context", signature, cancerType, annotationdf)
                first_df = survdf[survdf['ENSG'].isin(pruned['ENSG'])]

                original  = read_genes("original", signature, cancerType, annotationdf)
                second_df = survdf[survdf['ENSG'].isin(original['ENSG'])]
                
            elif panel_option=='discard': #discard vs original (purple expected)
                pruned  = read_genes("context", signature, cancerType, annotationdf)
                original  = read_genes("original", signature, cancerType, annotationdf)

                retained = original.loc[original['ENSG'].isin(pruned['ENSG']),:]
                discarded = original.loc[~original['ENSG'].isin(pruned['ENSG']),:]; 
                # retained vs discarded
                first_df = survdf[survdf['ENSG'].isin(retained['ENSG'])]
                second_df = survdf[survdf['ENSG'].isin(discarded['ENSG'])]
                # 
                
            elif panel_option=='de':
                pruned  = read_genes("context", signature, cancerType, annotationdf)
                first_df = survdf[survdf['ENSG'].isin(pruned['ENSG'])]
                
                pruned = read_genes("de", signature, cancerType, annotationdf)
                second_df = survdf[survdf['ENSG'].isin(pruned['ENSG'])]

            else:
                pruned  = read_genes("context", signature, cancerType, annotationdf)
                first_df = survdf[survdf['ENSG'].isin(pruned['ENSG'])]

                pruned = read_genes("frost", signature, cancerType, annotationdf)
                second_df = survdf[survdf['ENSG'].isin(pruned['ENSG'])]
            
            sub_df = second_df[(second_df['pvalue']<0.05) & (np.log(second_df['score'])>0)]
            if second_df.shape[0]==0:
                second_fraction = 0
            else:
                second_fraction = sub_df.shape[0]/second_df.shape[0]
            
            sub_df = first_df[(first_df['pvalue']<0.05) & (np.log(first_df['score'])>0)]
            if first_df.shape[0]==0:
                first_fraction = 0
            else:
                first_fraction = sub_df.shape[0]/first_df.shape[0]
            
            records[signature] = first_fraction - second_fraction
            reference_records[signature] = second_fraction
            context_records[signature] = first_fraction
                
        dfs.append(records)
        reference_dfs.append(reference_records)
        context_dfs.append(context_records)
    
    df = pd.DataFrame(dfs, index=cancerTypes)
    reference_df = None
    context_df = None
    reference_df = pd.DataFrame(reference_dfs, index=cancerTypes) 
    context_df = pd.DataFrame(context_dfs, index=cancerTypes) 
    
    return df, reference_df, context_df


def plot_survival_supplementary(df, figno, hr_direction='.', subType_option='.'):
    # rc('font', weight='bold')
    # rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

    cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)
    if subType_option=='BRCA':
        figsize = (20, 7)
    elif subType_option=='Brain':
        figsize = (20, 3)
    else:
        figsize = (12,7)
        
    # ax = pw.Brick(figsize=figsize)
    fig, axes = plt.subplots(1, 1, figsize=figsize, tight_layout=True) #
    ax = axes
    
    if hr_direction=='+':
        cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)
        title = 'Positively Significant HR'
    elif hr_direction=='>':
        cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)
        title = 'Significant HR>1'
    elif hr_direction=='-':
        cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)
        title = 'Negatively Significant HR'
    elif hr_direction=='.':
        cmap = sns.diverging_palette(300, 145, s=60, as_cmap=True)
        title = 'Significant HR'
        
    
    sns.heatmap(ax=ax, data=df.T, cmap=cmap, vmin=-0.4, vmax=0.4)  
    ax.tick_params(axis='y', rotation=0, labelsize=20)
    ax.tick_params(axis='x', rotation=90, labelsize=22)
    #ax.set_title(title)
    #ax.set_ylabel('CancerType')
    #ax.set_xlabel('Signature')

    # fig.savefig(f'{PROJECT_location}/Figures/Fig_{figno}.png', dpi=700)
    
    return None
