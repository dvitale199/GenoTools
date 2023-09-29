from scipy.stats import ncx2
import numpy as np

def calculate_inflation(pval_array, normalize=False, ncases=None, ncontrols=None):
    """Calculate lambda/genomic inflation for a GWAS summary statistic

    When given numpy array of P-values (e.g. a column of P-value in a dataframe),
    this function calculates the lambda or genomic inflation of the GWAS. If the
    inflation is below 0.9, the study is likely to have false negatives. If the
    inflation is above 1.05, the study is likely to have false positives. If there
    is a significant discrepancy between number of cases and controls in the study
    of interest, the inflation should be normalized (normalize=True).

    Parameters
    ----------
    pval_array : numpy array
            P-values in the GWAS summary statistic
    normalize : bool, default False
            Normalizes the data to 1000 cases and 1000 controls. Also called lambda1000.
            Recommended for studies with large discrepancy in number between cases and controls
    ncases : int, default None
            Number of cases. Required if normalize=True
    ncontrols : int, default None
            Number of controls. Required if normalize=True

    Example
    ---------
    PD_GWAS = pd.read_csv('PD.sumstat', sep = '\t')
    calculate_inflation(PD_GWAS['P-value'], normalize=True, ncases=49053, ncontrols=1411006)
    """
    if normalize is True and (ncases is None or ncontrols is None):
        raise NotImplementedError("If normalizing, please add ncases and ncontrols")
    num = ncx2.ppf(1-pval_array, 1, nc=0)
    denom = ncx2.ppf(0.5, 1, nc = 0)
    inflation = np.median(num)/denom
    if normalize:
        inflation1000 = 1 + (inflation -1) * (1/ncases+ 1/ncontrols)/(1/1000 + 1/1000)
        inflation = inflation1000
    return(inflation)