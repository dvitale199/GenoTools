import pandas as pd
import subprocess
import os
import sys
import gzip
import glob
import re
import numpy as np
import statsmodels.api as sm
import datetime
import matplotlib.pyplot as plt
!pip install qmplot
import qmplot
from sklearn import preprocessing
import gzip
# !pip install dask[dataframe]
import dask.dataframe as dd
import seaborn as sns
import matplotlib.pyplot as plt

# import from GenoTools.QC.qc.shell_do()
def shell_do(command, log=False, return_log=False, shell=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)
    if shell==True:
      res=subprocess.run(command.split(), stdout=subprocess.PIPE, shell=True)  
    else:
      res=subprocess.run(command.split(), stdout=subprocess.PIPE, shell=False)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))






# geno_path = <set path here>
geno_path = f'{geno_path}_pheno'
covar_path = 'DIAGNOMICS/ctcid_covariates-2021-07-28.csv'
covar_map_path = 'DIAGNOMICS/sample_ctcid_gender_map.csv'


# generate pcs to add to covariates
pca_cmd = f'{plink} --bfile {geno_path} --pca --out {geno_path}'
shell_do(pca_cmd)


# now create covariates
covar_map = pd.read_csv(covar_map_path)

covar_map.loc[:, 'IID'] = covar_map.sampleid.apply(lambda x: re.sub('0','', x, 1))
covar_map.loc[:,'FID'] = '0'
covar_map.loc[:,'sex_tmp'] = np.where(covar_map.gender.isna(),'0', covar_map.gender)
covar_map.loc[:,'sex'] = covar_map.loc[:,'sex_tmp'].map({'f': '2', 'm': '1','0':'0'})

covar = pd.read_csv(covar_path)
covar = covar.rename(columns={'ctcid':'IID'})
# read in pcs
pcs = pd.read_csv(f'{geno_path}.eigenvec',
                  sep='\s+',
                  header=None,
                  usecols=[0,1,2,3,4,5,6],
                  names=['FID','IID','PC1','PC2','PC3','PC4','PC5'],
                  dtype={'FID':str,'IID':str})

covar_merged = pcs.merge(covar, how='inner', on='IID')

now = pd.Timestamp('now')
covar_merged['dob'] = pd.to_datetime(covar_merged['date_of_birth'], format='%Y-%m-%d')
covar_merged['dob'] = covar_merged['dob'].where(covar_merged['dob'] < now, covar_merged['dob'] -  np.timedelta64(100, 'Y'))
covar_merged['age'] = (now - covar_merged['dob']).astype('<m8[Y]')
covar_merged.drop(columns=['smoking','date_of_birth','gender_x','sampleid','gender_y','sex_tmp','dob','ctcid'], inplace=True)



covar_merged2 = pcs.merge(covar_merged, on=['FID','ctcid'])

covar_merged2.to_csv(f'{geno_path}.covar', sep='\t', header=False, index=False)

# now let's run this thang!!!!!
# run logistic association on binary trait
assoc_cmd = f"\
{plink2} \
--bfile {geno_path} \
--covar {geno_path}.eigenvec \
--logistic \
--allow-no-sex \
--adjust cols=+qq \
--covar-variance-standardize \
--out {geno_path}"

shell_do(assoc_cmd)


logistic_adj = pd.read_csv(f'{geno_path}.PHENO1.glm.logistic.adjusted', sep='\s+')
print(logistic_adj[logistic_adj.BONF <= 0.05])
hits = logistic_adj.loc[logistic_adj.BONF <= 0.05, 'ID']


# clump
clump_cmd = f'\
{plink} --bfile {geno_path} \
--clump-p1 1e-3 \
--clump-r2 0.50 \
--clump-kb 250 \
--clump  {geno_path}.PHENO1.glm.logistic \
--clump-snp-field ID \
--clump-field P \
--allow-no-sex \
--out {geno_path}_clump'

shell_do(clump_cmd)

# Now, make scores
logistic = pd.read_csv(f'{geno_path}.PHENO1.glm.logistic', sep='\s+')
logistic_hits = logistic.loc[(logistic.TEST=='ADD')].copy()
logistic_hits.loc[:,'BETA'] = np.log10(logistic_hits.loc[:,'OR'])
logistic_hits = logistic_hits.dropna()
logistic_hits.loc[:,['ID','A1','BETA']].to_csv(f'{geno_path}_weights.tab', sep='\t', header=True, index=False)
clump = pd.read_csv(f'{geno_path}_clump.clumped', sep='\s+')
clump[['SNP']].to_csv(f'{geno_path}_clump.snps', sep='\t', header=False, index=False)
logistic_hits.loc[:,['ID','P']].to_csv(f'{geno_path}_pvals.snps', sep='\t', header=False, index=False)
!echo "s1 0 0.001" >> {geno_path}.range
!echo "s2 0 0.05" >> {geno_path}.range
!echo "s3 0 0.1" >> {geno_path}.range

prs_cmd = f'\
{plink} \
--bfile {geno_path} \
--score {geno_path}_weights.tab 1 2 3 \
--q-score-range {geno_path}.range {geno_path}_pvals.snps \
--extract {geno_path}_clump.snps \
--allow-no-sex \
--out {geno_path}.PRS \
--memory 50000'

shell_do(prs_cmd)

# read logistic for plotting
logistic = pd.read_csv(f'{geno_path}.PHENO1.glm.logistic', sep='\s+')
n_tests = logistic[logistic.TEST=='ADD'].shape[0]
corrected_significance = 5e-8
logistic_plot = logistic.loc[logistic.TEST=='ADD'].copy()

logistic_plot.dropna(how="any", axis=0, inplace=True)


# plot manhattan and QQ
# let's see if we can find a more customizable method to do this from-scratch
# here is something we could shoot for: https://dash.plotly.com/dash-bio/manhattanplot
ax = qmplot.manhattanplot(data=logistic_plot,
                          figname=f"{geno_path}.manhattan_plot.png",
                          genomewideline=corrected_significance,
                          is_annotate_topsnp=True,
                          suggestiveline=None,
                          hline_kws={"linestyle": "--", "lw": 1.3},
                          title="AMD WGS")

ax = qmplot.qqplot(data=logistic_plot["P"], figname=f"{geno_path}_QQ_plot.png")


sscore = pd.read_csv(f'{geno_path}.PRS.s1.profile', sep='\s+')
sscore_mean = sscore['SCORE'].mean()
sscore_std = sscore['SCORE'].std()
sscore['Z'] = (sscore['SCORE']-sscore_mean)/sscore_std
sscore.loc[:,'IID'] = sscore.IID.astype(str)


sns.set(style="darkgrid")

a = sscore[sscore['PHENO']==1]
b = sscore[sscore['PHENO']==2]

# plotting both distibutions on the same figure
fig = sns.kdeplot(a['Z'], shade=True, color="r")
fig = sns.kdeplot(b['Z'], shade=True, color="b")
plt.show()
fig.figure.savefig(f'{geno_path}_PRS_density.png')


##########################################
# linear regression for continuous trait #
##########################################
assoc_cmd = f"\
{plink2} \
--bfile {pheno_bed} \
--covar {pheno_bed}.covar \
--linear \
--adjust cols=+qq \
--covar-variance-standardize \
--out {pheno_bed}"

shell_do(assoc_cmd)

linear_adj = pd.read_csv(f'{pheno_bed}.PHENO1.glm.linear.adjusted', sep='\s+')
print(linear_adj[linear_adj.BONF <= 0.05])
hits = linear_adj.loc[linear_adj.BONF <= 0.05, 'ID']

# clump
clump_cmd = f'\
{plink} --bfile {pheno_bed} \
--clump-p1 1e-3 \
--clump-r2 0.50 \
--clump-kb 250 \
--clump  {pheno_bed}.PHENO1.glm.linear \
--clump-snp-field ID \
--clump-field P \
--out {pheno_bed}_clump'

shell_do(clump_cmd)


sscore = pd.read_csv(f'{pheno_bed}.PRS.s1.profile', sep='\s+')
sscore_mean = sscore['SCORE'].mean()
sscore_std = sscore['SCORE'].std()
sscore['Z'] = (sscore['SCORE']-sscore_mean)/sscore_std
sscore.loc[:,'IID'] = sscore.IID.astype(str)
covar_merged2.loc[:,'IID'] = covar_merged2.IID.astype(str)

sscore_merge = sscore.merge(covar_merged2, on='IID', how='inner')


# don't worry too much about this! useful to track linear relationship with PRS scores
this_formula = "PHENO ~ Z + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex + smoking_int"
res = sm.formula.glm(formula=this_formula, data=sscore_merge).fit() 
res.summary()

# plotting linear PRS. scatter instead of density plot
plt.scatter(sscore.Z, sscore.PHENO)
plt.title('Dry to Wet PRS')
plt.xlabel('PRS Z Score')
plt.ylabel('Phenotype')
plt.savefig(f'{pheno_bed}_PRS_scatter.png')

# same method to plot manhattan as above... let's make this better
n_tests = linear[linear.TEST=='ADD'].shape[0]
# corrected_significance = 0.05/n_tests
corrected_significance = 5e-8
linear_plot = linear.loc[linear.TEST=='ADD'].copy()
# linear_plot['-logp'] = -np.log10(linear_plot['P'])
# linear_plot
ax = qmplot.manhattanplot(data=linear_plot,
                          figname=f"{pheno_bed}.manhattan_plot.png",
                          genomewideline=corrected_significance,
                          is_annotate_topsnp=True,
                          suggestiveline=None,
                          hline_kws={"linestyle": "--", "lw": 1.3},
                          title="Dry to Wet Progression")

ax = qmplot.qqplot(data=linear_plot["P"], figname=f"{pheno_bed}_QQ_plot.png")