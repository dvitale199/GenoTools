import pandas as pd
import subprocess
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from matplotlib import cm 
import numpy as np
import os
import shutil
from umap import UMAP
import umap.plot
from sklearn import preprocessing, metrics
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
import plotly.express as px
import plotly
# from sklearn.externals import joblib
# import sklearn.externals.joblib as extjoblib
import joblib

#local imports
# from gwas.qc import QC
from gwas.utils import shell_do



def ancestry_prune(geno_path, out_path=None):
    geno_ancestry_prune_tmp = f'{geno_path}_ancestry_prune_tmp'
    
    if out_path:
        geno_ancestry_prune_out = out_path
    else:
        geno_ancestry_prune_out = f'{geno_path}_ancestry_prune'

    # prune geno_path for geno, maf, hwe, and palindromes
    geno_bim = pd.read_csv(f'{geno_path}.bim', sep='\t', header=None)

    # find and drop palindromes in geno_path .bim file
    geno_bim.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    palindromes = geno_bim.loc[((geno_bim.a1 == 'A') & (geno_bim.a2 == 'T')) | ((geno_bim.a1 == 'T') & (geno_bim.a2 == 'A')) | ((geno_bim.a1 == 'C') & (geno_bim.a2 == 'G')) | ((geno_bim.a1 == 'G') & (geno_bim.a2 == 'C'))]
    palindromes['rsid'].to_csv(f'{geno_path}_palindromes.snplist', header=False, index=False, sep='\t')

    plink_cmd1 = f'plink --bfile {geno_path}\
     --maf 0.05\
     --geno 0.01\
     --hwe 0.0001\
     --autosome\
     --allow-no-sex\
     --exclude {geno_path}_palindromes.snplist\
     --make-bed\
     --out {geno_ancestry_prune_tmp}' 

    # exclude high-LD regions
    plink_cmd2 = f'plink --bfile {geno_ancestry_prune_tmp}\
     --exclude range ref_data/hg19_exclusion_regions.txt\
     --autosome\
     --allow-no-sex\
     --make-bed\
     --out {geno_ancestry_prune_out}'

    cmds = [plink_cmd1, plink_cmd2]

    for cmd in cmds:
        shell_do(cmd)
    
    for suffix in ['bed','bim','fam','log']:
        try:
            os.remove(f'{geno_ancestry_prune_tmp}.{suffix}')
        except OSError:
            pass

        
def flash_pca(geno_path, out_name, dim=20):
    flashpca_cmd = f'\
    flashpca --bfile {geno_path}\
     -d {dim}\
     --outpc {out_name}.pcs\
     --outvec {out_name}.vec\
     --outval {out_name}.val\
     --outpve {out_name}.pve\
     --outload {out_name}.loadings\
     --outmeansd {out_name}.meansd'
    
    shell_do(flashpca_cmd)

    
def pca_projection(geno_path, inmeansd, inload, outproj):
    
    '''
    project new samples onto other PCs in flashpca format
    Must use the EXACT same snps (use get_common_snps() below)'''
    
    project_geno_pcs_cmd = f'\
    flashpca --bfile {geno_path}\
     --project\
     --inmeansd {inmeansd}\
     --inload {inload}\
     --outproj {outproj}\
     -v'

    shell_do(project_geno_pcs_cmd)



# def plot_pcs(labeled_pcs_df, plot_out=None, dim1='PC1', dim2='PC2'):
#     # now plot PCs
#     fig = plt.figure(figsize = (8,8))
#     ax = fig.add_subplot(1,1,1) 
#     ax.set_xlabel('Principal Component 1', fontsize = 15)
#     ax.set_ylabel('Principal Component 2', fontsize = 15)
#     ax.set_title('2 component PCA', fontsize = 20)
#     cmap = cm.get_cmap('tab10', 10)
#     targets = list(labeled_pcs_df.label.unique())
#     # targets = ['new', 'EUR']
#     colors = cmap(np.linspace(0, 1, len(targets)))
    
#     for i, target in enumerate(targets):
#         indicesToKeep = labeled_pcs_df['label'] == target
#         ax.scatter(labeled_pcs_df.loc[indicesToKeep, dim1]
#                    , labeled_pcs_df.loc[indicesToKeep, dim2]
#                    , c = [colors[i]]
#                    , s = 50)
#     ax.legend(targets)
#     ax.grid()
    
#     if plot_out:
#         fig.savefig(plot_out)
def plot_3d(labeled_df, color, symbol=None, plot_out=None, x='PC1', y='PC2', z='PC3', title=None):
    '''
    Input: 
        labeled_df: Pandas dataframe. labeled ancestry dataframe
        color: String. color of ancestry label. column name containing labels for ancestry in labeled_pcs_df
        symbol: String. symbol of secondary label (for example, predicted vs reference ancestry). default: None
        plot_out: String. filename to output filename for .png and .html plotly images
        x: String. column name of x-dimension
        y: String. column name of y-dimension
        z: String. column name of z-dimension

    Output:
        3-D scatterplot (plotly.express.scatter_3d). If plot_out included, will write .png static image and .html interactive to plot_out filename
        
    '''    
    fig = px.scatter_3d(labeled_df, x=x, y=y, z=z, color=color, symbol=symbol, title=title, color_discrete_sequence=px.colors.qualitative.Bold)
    # fig.show()

    if plot_out:
        fig.write_image(f'{plot_out}.png', width=1980, height=1080)
        fig.write_html(f'{plot_out}.html')



