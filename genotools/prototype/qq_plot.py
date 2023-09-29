import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def genomic_qqplot(pval_array, title=None):
    """Generates a QQ-Plot given a numpy array of p-values + title
    
    When given numpy array of P-values (e.g. a column of P-value in a dataframe),
    this function draws a QQ plot. QQ Plot is used to visualize the differences in
    distribution between observed P-value and expected P-value assuming no causal
    relationship between the phenotype and the variants studied.
    
    Ideal QQ plot should begin by following a perfect 1:1 (y=x) line but with observed
    P-value breaking upwards around -log10p == 6 or so. If breaking upward earlier,
    that implies high false positive rates. If breaking downward, that implies high
    false negative rate.

    QQ Plot should be paired with genomic inflation (lambda), appropriately normalized
    if necessary, to quantify distributional difference.

    This function produces a FacetGrid object with a single QQ plot and a y=x line.
    Future iterations should look into generating FacetGrid with multiple QQ plots
    and/or adding support for multiline QQ plot seen in LocusZoom plot (e.g. lines
    divided by MAF, data without specific genes like APOE for Alz dataset).

    Parameters
    ----------
    pval_array : numpy array
            P-values in the GWAS summary statistic
    title : str, default None
            Optional. Title to be recorded on the plot

    Example
    ---------
    genomic_qqplot(dat['P-value'], title="QQ-Plot")
    # to save
    plt.savefig('PD_gwas_qqplot.png')
    """
    logp = -np.log(pval_array)
    logp = np.sort(logp)
    intermid = 1/logp.shape[0]
    expected_logp = -np.log(np.linspace(intermid, 1, logp.shape[0]))
    expected_logp = np.sort(expected_logp)
    dat = pd.DataFrame({'p': logp, 'expected_p': expected_logp}, columns=['p', 'expected_p'])
    del logp
    del expected_logp
    max_val = np.amax(dat['expected_p'])
    qqplot_res = sns.FacetGrid(data=dat, height = 6)
    qqplot_res.map(sns.scatterplot, 'expected_p', 'p', edgecolor="none")
    qqplot_res.ax.axline((0, 0), slope=1, color="r", ls="-", zorder=0)
    qqplot_res.fig.suptitle(title)
    return(qqplot_res)