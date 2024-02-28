import scanpy as sc
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scvelo as scv
from scipy.sparse import issparse
from GOEA import GOEA
import os

signatures_path_=os.path.dirname(__file__)+'/'

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

def classify(adata, key, groupby='louvain', plot=True, average='mean', positive_selection=True, test=None, threshold=0.5, correct_for=None, **kwargs):
    """ Classifies groups into positive and negative acording to average of key.
    Key can be gene or obs.

    test can be 'l2fc' (default for genes), 'means' (default for signatures) or 't-test'
    """
    is_gene = True if key in adata.var_names else False if key in adata.obs else None
    if is_gene is None: raise ValueError('key not found in adata.var nor adata.obs')
    test = test if test is not None else 'l2fc' if is_gene else 'means'
    clusters = np.array(np.arange(len(pd.unique(adata.obs[groupby]))), dtype='str')
    n_clusters = len(clusters)
    mu = eval('np.'+average)

    if is_gene:
        Y = adata[:, key].X.A[:,0] if issparse(adata.X) else adata[:, key].X[:,0]
        Y = np.multiply(Y, adata.obs[correct_for] / np.mean(adata.obs[correct_for])) if correct_for is not None else Y
        X = [Y[adata.obs[groupby]==x] for x in clusters]
        X1 = [mu(Y[adata.obs[groupby]==x]) for x in clusters]
        M1 = mu(Y)
    else:
        Y = adata.obs[key]
        X = [Y[adata.obs.louvain==x] for x in clusters]
        X1 = [mu(Y[adata.obs.louvain==x]) for x in clusters]
        M1 = mu(Y)

    # # Smaller than mean
    if test=='means':
        selected_clusters = clusters[np.where(X1>M1)[0]] if positive_selection else clusters[np.where(X1<M1)[0]]
    elif test=='l2fc':
        # l2fc abs larger threshold
        l2fc = np.array([np.log2(x/M1) for x in X1])  # TODO: this throws NAN for negative Values, which can especially happen for scaled data or signatures
        l2fc[pd.isna(l2fc)]=0
        selected_clusters = clusters[np.where(l2fc>threshold)[0]] if positive_selection else clusters[np.where(l2fc<-threshold)[0]]
        # # log2fc compared to mean
        # pl.bar(np.arange(n_clusters), [np.log2(x/M1) for x in X1])
        # pl.xticks(np.arange(n_clusters), clusters)
        # pl.gca().axhline(0.5, color='r')
        # pl.gca().axhline(-0.5, color='r')
        # pl.show()
    else:
        raise NotImplementedError(test)
    # # 2-sided t-test
    # from scipy.stats import ttest_ind
    # pl.bar(np.arange(n_clusters), [ttest_ind(adata[adata.obs[groupby]==x, key].X.A[:,0], adata[adata.obs[groupby]!=x, key].X.A[:,0], equal_var=False)[1]  for x in clusters])
    # pl.xticks(np.arange(n_clusters), clusters)
    # pl.yscale('log')
    # pl.show()

    if plot:
        fig, axs = pl.subplots(1,3, figsize=[6*3,4*1]);
        if test=='means':
            axs[0].axhline(M1, color='red')
            axs[0].legend(['Mean '+key])
        elif test=='l2fc':
            axs[0].axhline(2**threshold*M1, color='red')
            axs[0].axhline(2**(-threshold)*M1, color='green')
            axs[0].legend(['L2FC '+key+' = '+str(threshold), 'L2FC '+key+' = -'+str(threshold)])
        axs[0].boxplot(X, showfliers=False, showmeans=True)
        axs[0].set_xticks(np.arange(n_clusters)+1)
        axs[0].set_xticklabels(clusters)
        axs[0].set_ylabel(key)
        axs[0].set_xlabel(groupby+' groups')
        sign = '+' if positive_selection else '-'
        scv.pl.scatter(adata, color=groupby, groups=selected_clusters, ax=axs[1], show=False, title=key+sign+' clusters', size=20, **kwargs)
        if 'perc' not in kwargs.keys(): kwargs['perc'] = [3,97]
        scv.pl.scatter(adata, color=key, ax=axs[2], show=False, title=key+ ' Umap', size=20, vmin=0, **kwargs)
        pl.show()
    return selected_clusters


def intersect_many(arrays):
    X = arrays[0]
    for Y in arrays:
        X = np.intersect1d(X,Y)
    return X


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

@shut_up_scanpy
def score_tumor_immune_cells(adata, signatures_path=signatures_path_, species='human'):
    # ImSigGenes immune tumor signatures
    tab=pd.read_excel(signatures_path+'/cell_type_markers/ImSigGenes_immunetumor.xlsx', skiprows=2, index_col=1)
    annot = dict()
    for ct in pd.unique(tab.Signature):
        annot[ct] = tab[tab.Signature==ct].index.values
    for ct in annot.keys():
        score_genes(adata, annot[ct], score_name=ct, species=species)

@shut_up_scanpy
def score_uhlitz_epithelial_cells(adata, signatures_path=signatures_path_, species='human'):
    tab=pd.read_excel(signatures_path+'cell_type_markers/colonoid_cancer_uhlitz_markers_revised.xlsx', skiprows=1, index_col=0, sheet_name=2)
    annot = dict()
    for ct in pd.unique(tab.index):
        annot[ct.replace('/', '_')] = np.array(tab[tab.index==ct].gene.values, dtype='str')
    for ct in annot.keys():
        score_genes(adata, annot[ct], score_name=ct, species=species)

@shut_up_scanpy
def score_ISC_stem(adata, signatures_path=signatures_path_, species='human'):
    # adds 'Stem_Lgr5_ISC-Munoz', 'Stem_Lgr5_ISC-Merlos' scores, mouse (!!!) derived and uppercased ISC signatures
    tab=pd.read_excel(signatures_path+'cell_type_markers/CRC-related_stem_cell_signatures.xlsx', header=0)
    tab = tab.drop(0)
    sigs = {'Stem_'+x: list(tab[x][~pd.isna(tab[x])].values) for x in tab.columns}
    for ct in ['Stem_Lgr5_ISC-Munoz', 'Stem_Lgr5_ISC-Merlos']:
        score_genes(adata, sigs[ct], score_name=ct, species=species)

@shut_up_scanpy
def score_hallmarks(adata, subset='organoid', signatures_path=signatures_path_, species='human'):
    # subset can be a list of hallmarks, 'organoid' (), 'CRC' (~18) or 'all' (50 scores)
    tab = pd.read_csv(signatures_path + 'msigdb/hallmark/h.all.v6.2.symbols.gmt', sep='\t', index_col=0, header=None).drop(1, axis=1).T
    hallsigs={hallmark : tab[hallmark][~pd.isna(tab[hallmark])].values for hallmark in tab.columns}
    if isinstance(subset, list):
        selection = subset
    elif subset == 'organoid':
        selection = ['HALLMARK_DNA_REPAIR', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_APOPTOSIS']
    elif subset == 'CRC':  # TODO this list is bugged, some entries do not exist
        selection = ['HALLMARK_DNA_REPAIR', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_APOPTOSIS',
        'HALLMARK_NOTCH_SIGNALING', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_HYPOXIA', 'HALLMARK_TGF_BETA_SIGNALING',
        'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_PROTEIN_SECRETION'
        'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
        'HALLMARK_P53_PATHWAY', 'HALLMARK_ANGIOGENESIS', 'HALLMARK_KRAS_SIGNALING_UP', 'HALLMARK_KRAS_SIGNALING_DN',
        'HALLMARK_GLYCOLYSIS']
    elif subset == 'all':
        selection = hallsigs.keys()
    else:
        raise ValueError('Please select a valid subset of hallmark to use. You can also choose "all".')
    for hm in selection:
        score_genes(adata, hallsigs[hm], score_name=hm, species=species)


@shut_up_scanpy
def score_Szabo_tcells(adata, subset='all', signatures_path=signatures_path_, species='human', key_add='Szabo_'):
    tab = pd.read_excel(signatures_path_+'cell_type_markers/Szabo_et_al_natcom_2019_S5_Tcellsig.xlsx', sheet_name=0, skiprows=1)
    signatures = {c: list(tab[c]) for c in tab.columns}
    for key in signatures.keys():
        score_genes(adata, signatures[key], score_name=f'{key_add}{key}', species=species)

@shut_up_scanpy
def score_CancerSEA(adata, subset='all', signatures_path=signatures_path_, species='human', key_add='CancerSEA_'):
    keys = ['Angiogenesis', 'Apoptosis', 'Cell_Cycle', 'Differentiation',
            'DNA_damage', 'DNA_repair', 'EMT', 'Hypoxia', 'Inflammation',
            'Invasion', 'Metastasis', 'Proliferation', 'Quiescence', 'Stemness']
    signatures = {s: list(pd.read_csv(f'{signatures_path}CancerSEA/{s}.txt', sep='\t')['GeneName']) for s in keys}
    if isinstance(subset, list):
        if all(np.isin(subset, keys)):
            selection = subset
        else:
            raise ValueError('Please select a valid subset of signature keys to use. You can also choose "all".')
    elif subset == 'all':
        selection = keys
    for key in keys:
        score_genes(adata, signatures[key], score_name=f'{key_add}{key}', species=species)

@shut_up_scanpy
def score_Guo_tcells(adata, signatures_path=signatures_path_, species='human', key_add='Guo_'):
    signatures = {}
    tab = pd.read_excel(signatures_path_+'cell_type_markers/Guo_DEGs_tcell_exhaustion.xlsx', sheet_name=0, skiprows=1)
    signatures['CD8_tcell_exhaustion'] = np.array(tab['Gene Symbol'][:-1], dtype=str)
    tab = pd.read_excel(signatures_path_+'cell_type_markers/Guo_DEGs_tcell_exhaustion.xlsx', sheet_name=1, skiprows=1)
    signatures['CD4_tcell_exhaustion'] = np.array(tab['Gene Symbol'][:-1], dtype=str)
    tab = pd.read_excel(signatures_path_+'cell_type_markers/Guo_Tcell_cluster_markers.xlsx', sheet_name=0, skiprows=1)
    signatures['CD8_tcells'] = np.array(tab['Gene Symbol'][:-1], dtype=str)
    tab = pd.read_excel(signatures_path_+'cell_type_markers/Guo_Tcell_cluster_markers.xlsx', sheet_name=1, skiprows=1)
    signatures['CD4_conventional_tcells'] = np.array(tab['Gene Symbol'][:-1], dtype=str)
    tab = pd.read_excel(signatures_path_+'cell_type_markers/Guo_Tcell_cluster_markers.xlsx', sheet_name=2, skiprows=1)
    signatures['CD4_regulatory_tcells'] = np.array(tab['Gene Symbol'][:-1], dtype=str)
    for key in signatures.keys():
        score_genes(adata, signatures[key], score_name=f'{key_add}{key}', species=species)

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


def get_genefamily_percentage(adata, key='MT-', start=True, name='mito'):
    keys = key if isinstance(key, list) else [key, '____ignore____']
    if start:
        family_genes = np.logical_or(*[adata.var_names.str.startswith(k) for k in keys])
    else:
        family_genes = np.logical_or(*[adata.var_names.str.endswith(k) for k in keys])
    if issparse(adata.X):
        adata.obs['percent_'+name] = np.sum(
            adata[:, family_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    else:
        adata.obs['percent_'+name] = np.sum(
            adata[:, family_genes].X, axis=1) / np.sum(adata.X, axis=1)

def get_mito_percentage(adata, species='human'):
    key = 'MT-' if species == 'human' else 'mt-'
    get_genefamily_percentage(adata, key=key, start=True, name='mito')

def get_ribo_percentage(adata, species='human'):
    key = specify_genes(['RPS', 'RPL'], species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='ribo')

def get_hemo_percentage(adata, species='human'):
    key = specify_genes(['HBA', 'HBB'], species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='hemo')
