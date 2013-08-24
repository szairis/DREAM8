import os
from scipy.io import savemat
from scipy.io import loadmat
import pandas as pd
import numpy as np
from pymatbridge import Matlab
from sklearn.ensemble import GradientBoostingRegressor

def hill(panel, prior_graph=[], lambdas=[], max_indegree=3, reg_mode='full', stdise=1, silent=0, maxtime=120):
    '''
    run_hill(panel)

    input: dataframe
    should be a T x N dataframe with T time points and N samples.

    output: dict containing key 'e' and key 'i' from Hill's code
    '''

    # start matlab
    mlab = Matlab(maxtime=maxtime)
    mlab.start()

    # slice out the data we want
    # just one time series for the cantone model
    D = np.transpose(panel.values, axes=(2,1,0))

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
    edge_prob = pd.DataFrame(out['e'], index=panel.minor_axis, columns=panel.minor_axis)
    edge_sign = pd.DataFrame(out['i'], index=panel.minor_axis, columns=panel.minor_axis)

    return (edge_prob, edge_sign)


def gbr(n_estimators=300, max_depth=3, learning_rate=0.1):
    '''
    run gradient boosting regression on design matrix.
    
    Input:
        design: TxN design matrix. should be a pandas dataframe with each network node
        as a column, and measured time points as rows.
        response: TxN response matrix. should be a pandas dataframe with each network node
        as a column, and each response variable as a row.

    Output
        Feature weights???
        Fit object???
    
    '''
    
    nodes = design.columns

    for target in nodes:
        regGBR[target] = {}
        for stim in stims:
            y = response[(response['Timepoint']>0) & (response['Stimulus']==stim)]
                .groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean()[target].values
            
            regGBR[target][stim] = GradientBoostingRegressor(n_estimators=n_estimators,
                                                             learning_rate=learning_rate,
                                                             max_depth=max_depth,
                                                             loss='ls')

            regGBR[target][stim].fit(X, y)



    return 0


def network_lasso(data, response_type='level', ground_truth=None, group_stimuli=False):
    '''
    lasso

    input:
        data
        response_type : (level, rate)
        ground_truth : adjacency matrix
        group_stimuli : binary
    '''
    from sklearn import preprocessing, linear_model, cross_validation, metrics

    training_dict = prepare_markov_data(data, response_type, group_stimuli)

    antibodies = data.columns[4:]
    stims = set(data['Stimulus'])

    if ground_truth:
        return A, scatter
    else:
        return A
