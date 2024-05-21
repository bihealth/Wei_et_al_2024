import scanpy as sc
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scvelo as scv
from scipy.sparse import issparse
# from GOEA import GOEA
import os


signatures_path_ = '../scrnaseq_signature_collection/'


def specify_genes(genes, species='human'):
    genes = genes if isinstance(genes, list) else list(genes) if isinstance(genes, np.ndarray) else [genes]
    if species is 'human':
        return [x.upper() for x in genes]
    elif species is 'mouse':
        return [x.capitalize() for x in genes]
    else:
        raise ValueError('Species '+species+' not known.')

def score_genes(adata, gene_list, score_name, species='human', **kwargs):
    gene_list_ = specify_genes(gene_list, species=species)
    sc.tl.score_genes(adata, gene_list_, score_name=score_name)

def shut_up_scanpy(func):
    def silent_wrapper(*args, **kwargs):
        verbosity = sc.settings.verbosity.value
        sc.settings.verbosity = 0
        func(*args, **kwargs)
        sc.settings.verbosity = verbosity
    silent_wrapper.__name__=func.__name__
    silent_wrapper.__annotations__=func.__annotations__
    silent_wrapper.__doc__=func.__doc__
    silent_wrapper.__defaults__=func.__defaults__
    return silent_wrapper

@shut_up_scanpy
def score_smillie_str_epi_imm(adata, signatures_path=signatures_path_, species='human'):
    tab=pd.read_excel(signatures_path+'/cell_type_markers/colonoid_cancer_uhlitz_markers_revised.xlsx', skiprows=1, index_col=0)
    score_genes(adata, np.array(tab.index[tab['Epithelial']==1].values, dtype='str'), score_name='epi_score', species=species)
    score_genes(adata, np.array(tab.index[tab['Stromal']==1].values, dtype='str'), score_name='str_score', species=species)
    score_genes(adata, np.array(tab.index[tab['Immune']==1].values, dtype='str'), score_name='imm_score', species=species)

def score_uhlitz_epithelial_cells(adata, signatures_path=signatures_path_, species='human'):
    tab=pd.read_excel(signatures_path+'/cell_type_markers/colonoid_cancer_uhlitz_markers_revised.xlsx', skiprows=1, index_col=0, sheet_name=2)
    annot = dict()
    for ct in pd.unique(tab.index):
        annot[ct.replace('/', '_')] = np.array(tab[tab.index==ct].gene.values, dtype='str')
    for ct in annot.keys():
        score_genes(adata, annot[ct], score_name=ct, species=species)
        
@shut_up_scanpy
def score_cell_cycle(adata, signatures_path=signatures_path_, species='human'):
    adatas = adata if isinstance(adata, list) else [adata]
    for i in range(len(adatas)):
        adata = adatas[i]
        # score cell cycle
        # cc score with genes from Kowalczyk, Monika S., et al. “Single-Cell RNA-Seq Reveals Changes in Cell Cycle and Differentiation Programs upon Aging of Hematopoietic Stem Cells.” Genome Research, vol. 25, no. 12, 2015, pp. 1860–72, doi:10.1101/gr.192237.115.
        cell_cycle_genes = [x.strip() for x in open(signatures_path+'cell_cycle_genes/regev_lab_cell_cycle_genes.txt')]
        cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
        # Split into 2 lists
        s_genes = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]

        # score
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        adatas[i] = adata
    return adatas[0] if len(adatas)==1 else adatas


def get_mito_percentage(adata, species='human'):
    key = 'MT-' if species == 'human' else 'mt-'
    get_genefamily_percentage(adata, key=key, start=True, name='mito')

def get_ribo_percentage(adata, species='human'):
    key = specify_genes(['RPS', 'RPL'], species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='ribo')

def get_hemo_percentage(adata, species='human'):
    key = specify_genes(['HBA', 'HBB'], species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='hemo')
