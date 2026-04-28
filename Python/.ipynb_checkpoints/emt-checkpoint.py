import numpy as np
import pandas as pd
from collections import Counter
import pickle
import glob
import os
import random
from itertools import chain
import networkx as nx
from scipy import stats

import seaborn as sns
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf as pdfbackend
from matplotlib import rc,rcParams
from matplotlib.lines import Line2D
from matplotlib_venn import venn3
# from venn import venn
import patchworklib as pw

from annot import *
from utils import read_genes

# get_context(cell_line): for all three conditions get the DE with fdr<0.1 and union
# get_consensus_emt(): common for 4 cell lines: A
# get_unique_context(cell_line): context - A
# checkout the expression data

CELL_LINES = ['a549', 'du145', 'mcf7', 'ovca420']
TREATMENTS = ['tgfb1', 'egf', 'tnf']
CELL_LINE_to_CANCERTYPE = {'a549':'LUAD', 'du145':'PRAD', 'mcf7':'BRCA', 'ovca420':'OV'}
CELL_LINE_to_tissue = {'a549':'lung', 'du145':'prostate', 'mcf7':'breast', 'ovca420':'ovarian'}


def plot_emt_treatment_genes(annotationdf):
    dgedict = {}
    cell_line_to_color={'a549': '#440154', 'du145': '#31688e', 'mcf7': '#35b779', 'ovca420': '#fde725'}
    axes1 = {}
    for cell_line in CELL_LINES:
        dges = {}
        for treatment in TREATMENTS:
            dge = pd.read_excel(f"{PROJECT_location}/cooketal/41467_2020_16066_MOESM4_ESM.xlsx", 
                          sheet_name=f"{cell_line}_{treatment}_dge")
            dge = dge[dge['fdr']<0.01]['Gene']
            dges[treatment] = set(dge)
        
        fig, ax = plt.subplots(1, 1, figsize=(5,4), tight_layout=True, dpi=1200) 
        # ax1 = pw.Brick(figsize=(5,4))
        v = venn3(list(dges.values()), ('TGFB1', 'EGF', 'TNF'), ax=ax)
        
        v.get_patch_by_id('100').set_color('grey')
        v.get_patch_by_id('010').set_color('grey')
        v.get_patch_by_id('001').set_color('grey')
        v.get_patch_by_id('110').set_color('grey')
        v.get_patch_by_id('011').set_color('grey')
        v.get_patch_by_id('101').set_color('grey')
        
        v.get_patch_by_id('111').set_color(cell_line_to_color[cell_line])
        for text in v.set_labels:
            text.set_fontsize(18)
            text.set_fontweight('bold')
        for text in v.subset_labels:
            text.set_fontsize(16)
        # fig.show()
        # ax.set_title(cell_line)
        
        axes1[cell_line] = ax
        
        fig.savefig(f'{PROJECT_location}/Figures/Fig_1A_venn_{cell_line}.png', dpi=700)
        dges = Counter(list(chain(*dges.values())))
        
    return axes1

def read_emt_treatment_genes(annotationdf):
    dgedict = {}
    cell_line_to_color={'a549': '#440154', 'du145': '#31688e', 'mcf7': '#35b779', 'ovca420': '#fde725'}
    for cell_line in CELL_LINES:
        dges = {}
        for treatment in TREATMENTS:
            dge = pd.read_excel(f"{PROJECT_location}/cooketal/41467_2020_16066_MOESM4_ESM.xlsx", 
                          sheet_name=f"{cell_line}_{treatment}_dge")
            dge = dge[dge['fdr']<0.01]['Gene']
            dges[treatment] = set(dge)
        
        dges = Counter(list(chain(*dges.values())))
        
        dgedict[cell_line] = []
        k=3
        for key, c in dges.items():
            if c>=k: dgedict[cell_line].append(key)
        print(k)
        context = read_genes("context", 'EMT', CELL_LINE_to_CANCERTYPE[cell_line], annotationdf)
        print(cell_line, len(dges), context.shape, len(dgedict[cell_line]))
        
    return dgedict

def read_emt_consensus_genes(dgedict):
    global_dges = Counter(list(chain(*dgedict.values())))
    k=3
    emt_genes = []
    for gene, c in global_dges.items():
        # count = 0
        # for cell_line, dges in dgedict.items():
        #     if gene in dges:
        #         count += 1
        if c>=k:
            emt_genes.append(gene)
    print(len(emt_genes))
    
    return emt_genes

def read_emt_context_genes(g0_treatment, g0a, annotationdf):
    dgedict = {}
    for cell_line, cell_emts in g0_treatment.items():
        dgedict[cell_line] = set(cell_emts).difference(g0a)
        
        context = read_genes("context", 'EMT', CELL_LINE_to_CANCERTYPE[cell_line], annotationdf)
        print(cell_line, len(cell_emts), len(dgedict[cell_line]), context.shape)
        
    return dgedict



def get_cross_correlations(gset1, gset2, countdf, nsample=750):
    gset1 = np.array(list(gset1))
    if len(gset1) > nsample:
        np.random.shuffle(gset1)
        # gset1 = gset1[:nsample]
        
    gset2 = np.array(list(gset2))
    if len(gset2) > nsample:
        np.random.shuffle(gset2)
        # gset2 = gset2[:nsample]
        
    # keep = []
    # for i in gset1:
    #     if i in countdf.index:
    #         keep.append(i)
            
    # for i in gset2:
    #     if i in countdf.index:
    #         keep.append(i)
            
    # df = countdf.loc[keep,:]
    df = countdf.loc[(countdf.index.isin(gset1))|(countdf.index.isin(gset2)), :]
    corrM = df.T.corr(method='spearman')
    
    corrM = corrM.where(np.triu(np.ones(corrM.shape), k=1).astype(bool))
    corrM = corrM.stack().reset_index()

    # its for another cross-correlation - for single gene in gset1 and multiple genes in gset2, the following logic wont work
    # corrM = corrM[~((corrM['level_0'].isin(gset1)) & (corrM['level_1'].isin(gset1)))]
    # corrM = corrM[~((corrM['level_0'].isin(gset2)) & (corrM['level_1'].isin(gset2)))]
    
    corrM = corrM[((corrM['level_0'].isin(gset1)) & (corrM['level_1'].isin(gset2)))|(((corrM['level_1'].isin(gset1)) & (corrM['level_0'].isin(gset2))))]
    
    return corrM[0]

# p/q for the co-expression networks

def compute_pq_in_coexpression(cell_emt_genes, foreset):
    corrdict = dict()
    corCutOff = 0.1
    nsample = 1000
    
    for cell_line in CELL_LINES:
        emt_genes = cell_emt_genes[cell_line]
        for treatment in TREATMENTS:
            coldf = pd.read_csv(
                f"{PROJECT_location}/cooketal/GSE147405_{cell_line}_{treatment}_KinaseScreen_metadata.csv", 
                                index_col=0)
            countdf = pd.read_csv(
                f"{PROJECT_location}/cooketal/GSE147405_{cell_line}_{treatment}_KinaseScreen_UMI_matrix.csv", 
                index_col=0)
            treateddf = countdf.loc[:,coldf['ConditionBroad']==f'{treatment.upper()}_Inhibited']

            non_emt_genes = treateddf.loc[~treateddf.index.isin(list(emt_genes)+foreset), :].index
            non_emt_genes = np.array(list(non_emt_genes))
            if len(non_emt_genes) > nsample:
                np.random.shuffle(non_emt_genes)
                non_emt_genes = non_emt_genes[:nsample]
            
            corrdict[f"{cell_line}_{treatment}-pq"] = []
            print(cell_line, treatment)
            for j,gene_i in enumerate(emt_genes):
                # print(gene_i)
                if j%100==0:
                    print(j)
                aa = get_cross_correlations(foreset, [gene_i], treateddf, nsample=nsample)
                bb = get_cross_correlations(non_emt_genes, [gene_i], treateddf, nsample=nsample)
                # print(aa)
                # print(bb)
                leng_aa = np.sum(aa > corCutOff)
                leng_bb = np.sum(bb > corCutOff)
                size_aa = len(foreset)
                size_bb = len(non_emt_genes)
                p = leng_aa / size_aa
                q = leng_bb / size_bb
                
                # if q!=0:
                #     p_by_q = 100
                # else:
                #     p_by_q = p/q
                p_by_q = p/q
                corrdict[f"{cell_line}_{treatment}-pq"].append(p_by_q)
                # print(p, size_aa, leng_aa, q, size_bb, leng_bb, p_by_q)
                # return non_emt_genes
    
    return corrdict

def plot_pq_fractions_in_coexpression(fractions):
    # cmap = sns.color_palette("Set3")
    # cmap = np.array([c for i,c in enumerate(cmap)])
    # cmap = cmap[np.array([2,4, 0,1])]
    # rc('font', weight='bold')
    # #rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    # #fig, ax = plt.subplots(1, 1, figsize=(6,5), tight_layout=True, dpi=700) 
    # ax = pw.Brick(figsize=(6,5))
    
    # labels = fractions.keys()
    # bplot = ax.boxplot(fractions.values(), patch_artist=True)
    # ax.set_xticklabels(labels, rotation = 90)
    # for patch, label, c in zip(bplot['boxes'], labels, cmap):
    #     print(c)
    #     patch.set_facecolor(c)
    
    # # for i, key in enumerate(fractions.keys()):
    # #     sns.kdeplot(ax=ax, data=fractions[key], color=cmap[i])
    # # plt.legend(labels=list(fractions.keys()))
    # ax.set_ylabel('$\it{p/q}$', fontsize=18, fontweight='bold') #30
    # ax.set_ylim((-5, 30))
    # ax.axhline(y=1, linestyle='--', color='grey')
    # ax.tick_params(axis='both', which='major', labelsize=18)#30
    # #fig.savefig(f'{PROJECT_location}/plots/{CUR_DST}/Fig_1c_ppi.png')

    cmap = sns.color_palette("Set3")
    cmap = np.array([c for i,c in enumerate(cmap)])
    cmap = cmap[np.array([2,4, 0,1])]
    # cmap = sns.color_palette("Paired")
    labels = fractions.keys()
    fig, ax = plt.subplots(1, 1, figsize=(8,7), tight_layout=True)
    # sns.boxplot(ax=ax, data=fractions, palette=cmap)
    bplot = ax.boxplot(fractions.values(), patch_artist=True)
    for i, (patch, label) in enumerate(zip(bplot['boxes'], labels)):
        patch.set_facecolor(cmap[i])
    ax.set_ylim((-1, 20))
    ax.set_xticklabels(fractions.keys(), rotation = 90)
    ax.axhline(y=1, linestyle='--', color='grey')
    ax.set_ylabel('$\it{p/q}$', fontsize=18, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=18)
    fig.savefig(f'{PROJECT_location}/Figures/Fig_1c_coppi.png')
    
    
    return ax

def compute_cross_coorelations(g0a, cell_emts):
    corrdict = dict()
    for cell_line in CELL_LINES:
        for treatment in TREATMENTS:
            #df = {}
            coldf = pd.read_csv(f"{PROJECT_location}/cooketal/GSE147405_{cell_line}_{treatment}_KinaseScreen_metadata.csv",
                                index_col=0)
            countdf = pd.read_csv(f"{PROJECT_location}/cooketal/GSE147405_{cell_line}_{treatment}_KinaseScreen_UMI_matrix.csv",
                                  index_col=0)
            treateddf = countdf.loc[:,coldf['ConditionBroad']==f'{treatment.upper()}_Inhibited']
            print(f"{cell_line}_{treatment}", countdf.shape, treateddf.shape[1])
            
            grands = np.array(range(treateddf.shape[0]))
            np.random.shuffle(grands)
            grands = treateddf.index[grands[:1000]]
            
            corrdict[f"{cell_line}_{treatment}-within_consensus"] = get_correlations(g0a, treateddf)
            corrdict[f"{cell_line}_{treatment}-within_context"] = get_correlations(cell_emts[cell_line], treateddf)
            corrdict[f"{cell_line}_{treatment}-consensus_cross_context"] = get_cross_correlations(cell_emts[cell_line], g0a, treateddf)
            corrdict[f"{cell_line}_{treatment}-consensus_cross_random"] = get_cross_correlations(cell_emts[cell_line], grands, treateddf)
            
            # consensus_cross_treatments = []
            # for i, (cell, gset1) in enumerate(cell_emts.items()):
            #     consensus_cross_treatments.append(get_cross_correlations(gset1, g0a, treateddf))
            # corrdict[f"{cell_line}_{treatment}-consensus_cross_treatments"] = pd.concat(
            #     consensus_cross_treatments, ignore_index=True)

            # reference_cross_treatments = []
            # for i, (cell, gset1) in enumerate(cell_emts.items()):
            #     reference_cross_treatments.append(get_cross_correlations(gset1, g0b, treateddf))
            # corrdict[f"{cell_line}_{treatment}-reference_cross_treatments"] = pd.concat(
            #     reference_cross_treatments, ignore_index=True)

            # grands = np.array(range(treateddf.shape[0]))
            # np.random.shuffle(grands)
            # grands = treateddf.index[grands[:1000]]
            # random_cross_treatments = []
            # for i, (cell, gset1) in enumerate(cell_emts.items()):
            #     random_cross_treatments.append(get_cross_correlations(gset1, grands, treateddf))
            # corrdict[f"{cell_line}_{treatment}-random_cross_treatments"] = pd.concat(
            #     random_cross_treatments, ignore_index=True)
            
            # across_cells = []
            # for i, (cell, gset1) in enumerate(cell_emts.items()):
            #     for j, (cell2, gset2) in enumerate(cell_emts.items()):
            #         if j<=i:
            #             continue
            #         across_cells.append(get_cross_correlations(gset1, gset2, treateddf))
            # corrdict[f"{cell_line}_{treatment}-across_cells"] = pd.concat(across_cells, ignore_index=True)

            
                
    return corrdict


def plot_cross_coexpressions_sub(corrdict):
    cmap = sns.color_palette("Set3")
    cmap = np.array([c for i,c in enumerate(cmap)])
    cmap = cmap[np.array([2,4, 0,1])]
    rc('font', weight='bold')
    fig, ax = plt.subplots(1, 1, figsize=(6,5), tight_layout=True, dpi=700) 
    axes = {}
    treatment = 'tgfb1'
    for cell_type, c in zip(CELL_LINES, cmap):
        # ax = pw.Brick(figsize=(8,6))
        
        fractions = dict((k, corrdict[k]) for k in corrdict.keys() if cell_type in k and treatment in k)
        # print(fractions.keys())
        labels = ['Within Pan-context', 'Within Context-specific', 'Context-spec X Pan-context', #'Reference X Context-specific', 
                  'Context-spec X Random']
        hatches = {'Within Pan-context': "|", 'Within Context-specific': "/", 'Context-spec X Pan-context': "+", #'Reference X Context-specific': "-", 
                   'Context-spec X Random': "."} 
        # print(len(fractions), len(labels))
        bplot = ax.boxplot(fractions.values(), patch_artist=True)
        ax.set_xticklabels(labels, rotation = 55)
        for patch, label in zip(bplot['boxes'], labels):
            patch.set_hatch(hatches[label])
            patch.set_facecolor(c)

        ax.set_ylabel('Co-expression', fontsize=25, fontweight='bold')
        ax.set_ylim((-0.5, 1))
        # ax.axhline(y=1, linestyle='--', color='grey')
        ax.tick_params(axis='both', which='major', labelsize=22)
        ax.set_title(cell_type, fontsize=30, fontweight='bold')
        
        axes[cell_type] = ax
    fig.savefig(f'{PROJECT_location}/Figures/Fig_1B_crosscoexp.png')
    
    # ax_all = (axes['a549']|axes['du145'])/(axes['mcf7']|axes['ovca420'])
    # ax_all.savefig(f'{PROJECT_location}/plots/{CUR_DST}/Fig_1c_ppi.png', dpi=700)
    
    return axes

def get_hppin_network():
    hppin = pd.read_csv(f'{PROJECT_location}/hPPIN/human_PPIN.txt', sep="\t", header=None)
    g = nx.from_pandas_edgelist(hppin, source=0, target=1)
    
    return g

def compute_fractions(g_set, fore_gset, hppin):

    # what fraction of the consensus should be connected to my set: p
    # what fraction of the cell
    debug_genes = {}
    neighbor_fractions = []
    absent = 0
    for node in g_set:
        
        if node not in hppin:
            absent += 1
            continue
            
        giter = hppin.neighbors(node)
        
        # compute p
        neigh_p = 0
        neigh_q = 0
        for i in giter:
            if i in fore_gset:
                neigh_p += 1 #aa
            if i in hppin:
                neigh_q += 1 #bb
        
        if neigh_p==0:
            continue
        
        p = neigh_p/(len(fore_gset))
        q = neigh_q/(len(hppin)-absent)
        
        neighbor_fractions.append(p/q)
        if q>p:
            debug_genes[node] = (p, q)
        
    print(debug_genes)
    return np.array(neighbor_fractions)


def get_fraction_data(cell_emts, g0a, g0b):
    hppin = get_hppin_network()
    
    # fractions = {}
    # for i in cell_emts:
    #     #print(i, list(cell_emts[i])[:5])
    #     fractions[f"{i}-consensus"] = compute_fractions(cell_emts[i], g0a, hppin)
    
    # fractions['emt'] = compute_fractions(g0b['Symbol'], hppin)
    
    fractions = {}
    for i in cell_emts:
        #print(i, list(cell_emts[i])[:5])
        fractions[f"{i}"] = compute_fractions(cell_emts[i], g0a, hppin)
        # fractions[f"{i}-original"] = compute_fractions(cell_emts[i], g0b, hppin)
    
    # fractions['emt-consensus'] = compute_fractions(g0a, hppin)
    # fractions['emt-original'] = compute_fractions(g0b['Symbol'], hppin)
    
    return fractions

def plot_nw_fractions(fractions):
    cmap = sns.color_palette("Set3")
    cmap = np.array([c for i,c in enumerate(cmap)])
    cmap = cmap[np.array([2,4, 0,1])]
    # rc('font', weight='bold')
    # rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    fig, ax = plt.subplots(1, 1, figsize=(6,5), tight_layout=True, dpi=700) 
    # ax = pw.Brick(figsize=(6,5))
    
    labels = fractions.keys()
    bplot = ax.boxplot(fractions.values(), patch_artist=True)
    ax.set_xticklabels(labels, rotation = 90)
    for patch, label, c in zip(bplot['boxes'], labels, cmap):
        print(c)
        patch.set_facecolor(c)
    
    # for i, key in enumerate(fractions.keys()):
    #     sns.kdeplot(ax=ax, data=fractions[key], color=cmap[i])
    # plt.legend(labels=list(fractions.keys()))
    ax.set_ylabel('$\it{p/q}$', fontsize=18, fontweight='bold') #30
    ax.set_ylim((-5, 30))
    ax.axhline(y=1, linestyle='--', color='grey')
    ax.tick_params(axis='both', which='major', labelsize=18)#30
    fig.savefig(f'{PROJECT_location}/Figures/Fig_1c_ppi.png')
    
    return ax

