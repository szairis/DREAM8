import os
from scipy.io import savemat
from scipy.io import loadmat
import pandas as pd
import numpy as np
from pymatbridge import Matlab


def hill(df, prior_graph=[], lambdas=[], max_indegree=3, reg_mode='full', stdise=1, silent=0, maxtime=120):
    '''
    run_hill(df)

    input: dataframe
    should be a T x N dataframe with T time points and N samples.

    output: dict containing key 'e' and key 'i' from Hill's code
    '''

    mlab = Matlab(maxtime=maxtime)
    mlab.start()

    # slice out the data we want
    # just one time series for the cantone model
    D = df.values.T

    # .mat shuttle files
    inPath = os.path.join('..', 'cache', 'dbn_wrapper_in.mat')
    outPath = os.path.join('..', 'cache', 'dbn_wrapper_out.mat')

    # save the matlab object that the DBN wrapper will load
    # contains all the required parameters for the DBN code
    savemat(inPath, {'D' : D,
                     'max_indegree' : max_indegree,
                     'prior_graph' : prior_graph,
                     'lambdas' : lambdas,
                     'reg_mode' : reg_mode,
                     'stdise' : stdise,
                     'silent' : silent})

    # DBN wrapper just needs an input and output path
    args = {'inPath' : inPath, 'outPath' : outPath}

    # call DBN code
    res = mlab.run_func('dbn_wrapper.m', args, maxtime=maxtime)

    mlab.stop()

    out = loadmat(outPath)
    edge_prob = pd.DataFrame(out['e'], index=df.columns, columns=df.columns)
    edge_sign = pd.DataFrame(out['i'], index=df.columns, columns=df.columns)

    return (edge_prob, edge_sign)