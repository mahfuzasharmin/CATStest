import numpy as np
import pandas as pd
from collections import Counter

import pickle
import glob
import os

from annot import *

from utils import read_genes
import networkx as nx
import networkx.algorithms.community as nx_comm

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, CirclePolygon
from matplotlib.lines import Line2D
from matplotlib import rc, rcParams


def build_network(cancerType, signature, annotationdf, outfile=None):
    hppin = pd.read_csv(f'{PROJECT_location}/network/human_PPIN.txt', sep="\t", header=None)

    original = read_genes('original', signature, cancerType, annotationdf)
    context = read_genes('context', signature, cancerType, annotationdf)
    all_genes = set(original['Symbol']).union(context['Symbol'])
    g = nx.from_pandas_edgelist(hppin, source=0, target=1)
    big_g = g.subgraph(all_genes)
    both_g = g.subgraph(set(original['Symbol']).union(context['Symbol']))
    
    return both_g

def show_network(g_, original, node_positions=None, legend_loc='upper right', figno='4b'):
    
    group_to_color = {0: 'turquoise', 1: 'deepskyblue'}
    label_to_group = {'Original': 0, 'Context': 1}
    custom_circles = [CirclePolygon((1,1), 3, color=group_to_color[0]),
                  CirclePolygon((1,1), 3, color=group_to_color[1])]
    # custom_lines = [Line2D([0], [0], lw=4, color=group_to_color[0]),
    #               Line2D([0], [0], lw=4, color=group_to_color[1])]
    
    color_map = []
    for node in g_:
        if node in original['Symbol'].values:
            color_map.append(group_to_color[label_to_group['Original']])
        else: 
            color_map.append(group_to_color[label_to_group['Context']])      

    if node_positions is None:
        node_positions = nx.kamada_kawai_layout(g_)
    fig, axes = plt.subplots(1, 1, figsize=(10,8), tight_layout=True) 
    # rc('font', weight='bold')
    # rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    
    ax = axes
    nx.draw(g_, ax=ax, pos=node_positions, node_color=color_map, node_size=150, with_labels=True, font_size=8, 
            edge_color='grey', width=0.5)
    
    for label in label_to_group:
        ax.plot([0],[0],color=group_to_color[label_to_group[label]],label=label)
    ax.legend(custom_circles, label_to_group.keys(), loc=legend_loc, fontsize = 12)
    fig.tight_layout()
    
    fig.savefig(f'{PROJECT_location}/Figures/Fig_{figno}.png', dpi=700)
    return node_positions


def show_network_without_subplots(g_, node_positions, original): #not used anymore
    color_map = []
    for node in g_:
        if node in original['Symbol'].values:
            color_map.append('turquoise')
        else: 
            color_map.append('deepskyblue')      
    
    nx.draw(g_, pos=node_positions, node_color=color_map, node_size=100, with_labels=True, font_size=8, 
            edge_color='grey', width=0.2)
    return None

def show_louvian_communities_together(both_g, node_positions, original):
    
    communities = nx_comm.louvain_communities(both_g, seed=123)
    fig, axes = plt.subplots(1, 1, figsize=(8,5), tight_layout=True) 
    for i in communities:
        if len(i)>1:
            sub_g = both_g.subgraph(i)
            show_network_without_subplots(sub_g, node_positions, original)

    return None

def show_network_with_subplots(g_, node_positions, original):
    color_map = []
    for node in g_:
        if node in original['Symbol'].values:
            color_map.append('turquoise')
        else: 
            color_map.append('deepskyblue')      

    fig, axes = plt.subplots(1, 1, figsize=(6,5), tight_layout=True) 
    nx.draw(g_, ax=axes, pos=node_positions, node_color=color_map, node_size=100, with_labels=True, font_size=8, 
            edge_color='grey', width=0.2)
    return None


def show_louvian_communities_separately(both_g, node_positions, original):
    
    communities = nx_comm.louvain_communities(both_g, seed=123)
    for i in communities:
        if len(i)>10:
            sub_g = both_g.subgraph(i)
            show_network_with_subplots(sub_g, node_positions, original)

    return None


def get_sizeable_communities(g_):
    communities = nx_comm.louvain_communities(g_, seed=123)
    retained_communities = []
    for i in communities:
        if len(i)>10:
            retained_communities.append(i)
            
    return retained_communities


def read_enrichments(cancerType, signature='EMT'):
    context_e = pd.read_excel(f"{PROJECT_location}/enrichments_context/{signature}_GO.xlsx", 
                              sheet_name=cancerType)

    original_e = pd.read_excel(f"{PROJECT_location}/enrichments_reference/cancersea_GO.xlsx", sheet_name=signature)
    context_e = context_e[~context_e['ID'].isin(original_e['ID'])]
    
    return context_e


def show_terms_for_nth_community(main_g, nth, node_positions, cancerType, signature, original):
    context_e = read_enrichments(cancerType, signature=signature)
    
    retained_communities = get_sizeable_communities(main_g)
    
    current_community = retained_communities[nth]
    sub_g = main_g.subgraph(current_community)
    
    #show_network_with_subplots(sub_g, node_positions, original)
    assigned_df = context_e[context_e['geneID'].apply(
        lambda x: len(current_community.intersection(x.split('/')))/len(x.split('/'))>0.50)]
    
    termid_to_genelist = assigned_df['geneID'].apply(lambda x: x.split('/')).to_dict()
    #print(termid_to_genelist)
    return assigned_df, termid_to_genelist

def get_nth_community(main_g, nth):
    retained_communities = get_sizeable_communities(main_g)
    current_community = retained_communities[nth]
    sub_g = main_g.subgraph(current_community)
    
    return sub_g

def show_network_focusing_enrichments(
    g_, node_positions, original, 
    termid_to_genelist, termid_to_color, term_to_termid, legend_loc, figno):
    
    color_map = []
    for node in g_:
        if node in original['Symbol'].values:
            color_map.append('turquoise')
        else: 
            community_id = None
            for i, community in termid_to_genelist.items():
                if node in community:
                    community_id = i
                    #break
            if community_id is None:
                color_map.append('deepskyblue')   
            else:
                color_map.append(termid_to_color[community_id])

    
    # rc('font', weight='bold')
    # rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    
    fig, axes = plt.subplots(1, 1, figsize=(10,8), tight_layout=True) 
    ax = axes
    nx.draw(g_, ax=ax, pos=node_positions, node_color=color_map, node_size=150, with_labels=True, font_size=15, 
            edge_color='grey', width=0.2)
    
    for label in term_to_termid:
        ax.plot([0],[0],color=termid_to_color[term_to_termid[label]],label=label)
    plt.axis('off')
    fig.set_facecolor('w')
    plt.legend(loc=legend_loc, fontsize = 12)
    fig.tight_layout()
    
    fig.savefig(f'{PROJECT_location}/Figures/Fig_{figno}.png', dpi=700)
    return None

def write_network_edges(main_g, outfile, nthlist, cancerType, signature, original, annotationdf):
    
    # context attribute
    context = read_genes('context', signature, cancerType, annotationdf)
    def is_context(n):
        if n in context['Symbol'].values:
            label = 'Context-specific' 
        else:
            label = 'Context-agnostic'
        return label
    
    #comminuty attribute
    gene_to_community = {}
    retained_communities = get_sizeable_communities(main_g)
    for n, nth in enumerate(retained_communities):
        for gene in nth:
            gene_to_community[gene] = n
    for gene in main_g:
        if gene not in gene_to_community:
            gene_to_community[gene] = -1
    
    #enrichment attribute under each community
    gene_to_enrichment = {}
    for n in nthlist:
        sub_g = get_nth_community(main_g, n)
        term_df, termid_to_genelist = show_terms_for_nth_community(main_g, n, None, cancerType, signature, original)
        for gene in sub_g:
            for termid, genelist in termid_to_genelist.items():
                term = term_df.loc[termid, 'Description']
                if gene in genelist:
                    gene_to_enrichment[gene] = term
                    break
    
    with open(f'{outfile}.edges.csv', 'w') as fp:
        fp.write(f'Source,Target,Attributes,Attributes,community,community,term,term\n')
        for e in main_g.edges:
            src_att = is_context(e[0])
            tgt_att = is_context(e[1])
            src_comm = gene_to_community[e[0]]
            tgt_comm = gene_to_community[e[1]]
            src_term = gene_to_enrichment[e[0]] if e[0] in gene_to_enrichment else ''
            tgt_term = gene_to_enrichment[e[1]] if e[1] in gene_to_enrichment else ''
            fp.write(f'{e[0]},{e[1]},{src_att},{tgt_att},{src_comm},{tgt_comm},{src_term},{tgt_term}\n')
            
    
    return None



