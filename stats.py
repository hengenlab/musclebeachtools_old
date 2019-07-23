import pandas as pd
import numpy as np
from scipy import stats as scstats
import pdb
from itertools import combinations as comb

def rmANOVA(df, alpha = 0.05):
    '''Code is taken directly from:
    http://pythonpsychologist.tumblr.com/post/139246503057/repeated-measures-anova-using-python
    on 05/17/2018 and transcribed [and edited where necessary] by KBH.
    Your pandas dataframe should be organized as follows (Subid is subject code, Xn are measurements):
                 Subid  X1  X2  X3
             0      1   6   8  10
             1      2   4   5   6
             2      3   5   5   5
             3      4   1   2   3
             4      5   0   1   2
             5      6   2   3   4'''

    def calc_grandmean(data, columns):
       '''Takes a pandas dataframe and calculates the grand mean data = dataframe columns = list of column names with the response variables.'''

       gm = np.mean(data[columns].mean())
       return gm


    # KBH ADDED: get a list of all column names:
    colnames = df.columns.values.tolist()

    #Grand mean
    grand_mean      = calc_grandmean(df, colnames)
    # KBH edited:
    #df['Submean']         = df[colnames].mean(axis=1)
    submean         = df[colnames].mean(axis=1)
    column_means    = df[colnames].mean(axis=0)

    n = df.shape[0] # len(df['Subid']) - KBH edit. this just needs the number of "subjects"
    k = len(colnames)
    #Degree of Freedom
    ncells = df[colnames].size


    dftotal = ncells - 1
    dfbw    = k - 1
    dfsbj   = n - 1 # len(df['Subid']) - 1 KBH EDIT
    dfw     = dftotal - dfbw
    dferror = dfw - dfsbj

    #We start with SS between. SS between is the sum of squared deviations of the sample means from the grand mean multiplied by the number of observations:
    SSbetween = sum(n*[(m - grand_mean)**2 for m in column_means])

    # We continue with SS within (the sum of squared deviations within each sample):
    SSwithin = sum(sum([(df[col] - column_means[i])**2 for i, col in enumerate(df[colnames])]))

    # SS subjects : The sum of squared deviations of the subject means from the grand mean multiplied by the number of conditions (k)
    SSsubject = sum(k*[(m -grand_mean)**2 for m in submean ]) #df['Submean']])

    #SS error : The sum of squared deviations due to sampling error
    SSerror = SSwithin - SSsubject

    #We can also calculate the SS total (i.e., The sum of squared deviations of all observations from the grand mean):
    SStotal = SSbetween + SSwithin

    # After we have calculated the Mean square error and Mean square between we can obtain the F-statitistica:
    #MSbetween
    msbetween = SSbetween/dfbw

    #MSerror
    mserror = SSerror/dferror

    #F-statistic
    F = msbetween/mserror

    #By using SciPy we can obtain a p-value. We start by setting our alpha to .05 and then we get our p-value.
    p_value = scstats.f.sf(F, 2, dferror)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Below is all by KBH:
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

    if p_value < alpha:
        # Generate all possible combinations of the conditions for post hoc tests
        pairs   = np.asarray(list(comb(colnames, 2)))
        pv      = np.zeros( ( np.shape(pairs)[0], 1) )

        icount = 0
        for i in pairs:
            # just run t tests and correct for multiple comparisons. Post hocs don't apply very well to rmANOVA
            st, pval    = scstats.ttest_rel(df[i[0]], df[i[1]])

            # correct for multiple comparisons
            pv[icount]  = pval*np.shape(pairs)[0]
            icount +=1

        feed = np.append( pairs, pv, 1 )
    else:
        feed = NaN


    return p_value, feed
