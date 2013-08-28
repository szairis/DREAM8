import pdb
import os
import pandas as pd
import numpy as np
from pymatbridge import Matlab

from utilities import prepare_markov_data, introduce_inhibs, score_network

def network_hill(panel, prior_graph=[], lambdas=[], max_indegree=3, reg_mode='full', stdise=1, silent=0, maxtime=120):
    '''
    run_hill(panel)

    input: dataframe
    should be a T x N dataframe with T time points and N samples.

    output: dict containing key 'e' and key 'i' from Hill's code
    '''
    from scipy.io import savemat
    from scipy.io import loadmat

    # start matlab
    mlab = Matlab(maxtime=maxtime)
    mlab.start()

    # .mat shuttle files
    # add path check
    inPath = os.path.join('..', 'cache', 'dbn_wrapper_in.mat')
    outPath = os.path.join('..', 'cache', 'dbn_wrapper_out.mat')
    D = np.transpose(panel.values)
    num_rows = np.shape(D)[0]
    num_cols = np.shape(D)[1]
    D = np.reshape(D, (num_rows, num_cols, 1))

    # save the matlab object that the DBN wrapper will load
    # contains all the required parameters for the DBN code
    savemat(inPath, {"D" : D,
                     "max_indegree" : max_indegree,
                     "prior_graph" : prior_graph,
                     "lambdas" : lambdas,
                     "reg_mode" : reg_mode,
                     "stdise" : stdise,
                     "silent" : silent})

    # DBN wrapper just needs an input and output path
    args = {"inPath" : inPath, "outPath" : outPath}

    # call DBN code
    res = mlab.run_func('dbn_wrapper.m', args, maxtime=maxtime)

    mlab.stop()

    out = loadmat(outPath)
    edge_prob = pd.DataFrame(out['e'], index=panel.columns, columns=panel.columns)
    edge_sign = pd.DataFrame(out['i'], index=panel.columns, columns=panel.columns)

    return (edge_prob, edge_sign)


def gbr(data, n_estimators=300, max_depth=3, learning_rate=0.1, inhib_targets=None, perfect=True):
    '''
    run gradient boosting regression on design matrix.
    
    Input:
        data
        design: TxN design matrix. should be a pandas dataframe with each network node
        as a column, and measured time points as rows.
        response: TxN response matrix. should be a pandas dataframe with each network node
        as a column, and each response variable as a row.

    Output
        Feature weights???
        Fit object???
    
    '''
    from sklearn.ensemble import GradientBoostingRegressor
    
    # model interventions if supplied an inhib_targets dict
    if inhib_targets:
        training_dict = prepare_markov_data(introduce_inhibs(data, inhib_targets=inhib_targets, perfect=perfect), response_type, group_stimuli)
    else:
        training_dict = prepare_markov_data(data, response_type, group_stimuli)

    nodes = design.columns

    for target in nodes:
        regGBR[target] = {}
        for stim in stims:
            y = response[(response['Timepoint']>0) & (response['Stimulus']==stim)] \
                .groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean()[target].values
            
            regGBR[target][stim] = GradientBoostingRegressor(n_estimators=n_estimators,
                                                             learning_rate=learning_rate,
                                                             max_depth=max_depth,
                                                             loss='ls')

            regGBR[target][stim].fit(X, y)



    return 0


def network_lasso(data, response_type='level', ground_truth=None, inhib_targets=None, perfect=True, group_stimuli=False):
    '''
    do lasso. automatically do CV to find best alpha.

    input:
        data
        response_type : (level, rate)
        ground_truth : adjacency matrix
        group_stimuli : binary
    '''
    from sklearn import preprocessing, linear_model, cross_validation, metrics

    # model interventions if supplied an inhib_targets dict
    if inhib_targets:
        training_dict = prepare_markov_data(introduce_inhibs(data, inhib_targets=inhib_targets, perfect=perfect), response_type, group_stimuli)
    else:
        training_dict = prepare_markov_data(data, response_type, group_stimuli)

    antibodies = [col for col in data.columns if col not in ['Cell Line', 'Inhibitor', 'Stimulus', 'Timepoint']]
    stims = set(data['Stimulus'])

    # fit lasso for each (X,Y) pair
    A = {}
    for key in training_dict:
        X = training_dict[key][0]
        Y = training_dict[key][1]
        preprocessing.StandardScaler().fit_transform(X)

        A[key] = pd.DataFrame(np.zeros((X.shape[1], X.shape[1])), columns=X.columns, index=X.columns)

        for col in Y.columns:
            #print col
            # check if col is not all the identical
            if len(set(Y[col])) > 1:
                rgn = linear_model.LassoCV(verbose=False).fit(X, Y[col])
                if np.max(rgn.coef_) != 0:
                    A[key].ix[:,col] = np.abs(rgn.coef_) / np.abs(rgn.coef_).max()
            else:
                A[key].ix[:,col] = np.zeros((X.shape[1],))

    if ground_truth:
        auc = {}
        for key in training_dict:
            auc[key] = score_network(A[key], ground_truth)
        return A, auc 
    else:
        return A
