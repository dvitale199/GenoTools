import numpy as np
from scipy.stats import norm


def zscore_pval_conversion(zscores=None, pvals=None, stats=None):

    # neither zscore or pvals provided provided
    if zscores is None and pvals is None:
        print('Conversion Failed!')
        print('Either p-values or z-scores must be provided')
    
    # both zscores and pvals provided
    elif zscores is not None and pvals is not None:
        print('Conversion Failed!')
        print('Provide only p-values or z-scores, not both')
    
    # pvals provided but stats not provided to determine sign of zscore
    elif pvals is not None and stats is None:
        print('Conversion Failed!')
        print('Stats must be provided when going from p-values to z-scores')
    
    else:
        # convert pvals to zscores using stats to get proper sign
        if zscores is None:
            z = np.where(stats > 0, norm.isf(pvals/2), -norm.isf(pvals/2))
            return z
        # convert zscores to pvals
        if pvals is None:
            p = 2*norm.sf(abs(zscores))
            return p