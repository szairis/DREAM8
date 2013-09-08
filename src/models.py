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

def do_gbr(X, Y, n_estimators=100, learning_rate=0.1, max_depth=5, verbose=False):
    '''does gbr on design matrix.
    returns dict regGBR, one GBR for each column in the target (Y) matrix

    do this and then do network_gbr, which will give you the a-matrix from the feature importances
    '''
    regGBR = {}

    for target in Y.columns:
        
        if verbose:
            print target
        
        # get target values
        y = Y[target].values
        
        regGBR[target] = GradientBoostingRegressor(n_estimators=n_estimators,
                                                   learning_rate=learning_rate,
                                                   max_depth=max_depth,
                                                   loss='ls')
        regGBR[target].fit(X, y)
        
    return regGBR

def build_adj_matrix(regGBR, node_list, stims):   
    '''take as input the regGBR object, build an adjacency matrix out of each stimulus
    '''
    num_nodes = len(node_list)
    adj_dict = {}
    for stim in stims:
        adj = np.zeros((num_nodes, num_nodes), dtype='f')
        for nidx, node in enumerate(node_list):
            adj[:, nidx] = regGBR[stim][node].feature_importances_[:num_nodes]
        adj_dict[stim] = pd.DataFrame(adj, index=node_list, columns=node_list)
        
    return adj_dict

def timeseries_gbr(regGBR, scalar, data, node_list, stims, inhibs, times):
    '''takes a regGBR object (output of do_GBR)
    and predicts timeseries for inhibited nodes

    you supply the normalization scalar, original data, list of stimulii,
    list of inhibitors present in data, and prediction timepoints

    returns pred_dict, after which you would generally call write_midas, looping
    over the node_list
    '''


    num_inhibs = len(inhibs)

    pred_dict = {}


    for test_inhib in node_list:
        print test_inhib
        pred_dict[test_inhib] = {}
        for stim in stims:
            print stim
            # set up new df to use, and fill t=0 values
            pred_df = pd.DataFrame(np.zeros((len(times), len(node_list))), index=times, columns=node_list)
            pred_df.ix[0, :] = data.groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix['None', stim, 0]
            pred_df.ix[0, test_inhib] = 0
            
            # loop over times
            for tidx in range(1,len(times)):
                time = times[tidx]
                
                # get covariates for this time step and scale
                # covariates_df = scaler.transform(pred_df.ix[times[tidx-1], :])
                covariates_df = ((pred_df.ix[times[tidx-1], :]) - scalar.mean_[:-num_inhibs]) / scalar.std_[:-num_inhibs]

                # zero out covariate we are inhibiting
                try:
                    covariates_df.ix[test_inhib_targets[test_inhib]] = 0
                except:
                    pass
                
                num_cov = len(pred_df.columns)
                covariates = np.zeros((num_cov + num_inhibs,))
                covariates[:num_cov] = covariates_df.values

                # loop over proteins to get values for current time step
                for p in node_list:
                    pred_df.ix[time, p] = insilico_regGBR[p].predict(covariates)
                
                # zero out covariate we are inhibiting, again
                pred_df.ix[time, test_inhib] = 0
        
            pred_dict[test_inhib][stim] = pred_df

    return pred_dict

def network_gbr_cv(X, Y, verbose=False, n_estimators=100, learning_rate=0.1, max_depth=5):
    
    n_folds = 5
    kf = list(cross_validation.KFold(X.shape[0], n_folds=n_folds, shuffle=True))
    
    regGBR = {}
    test_score = {}

    for target in Y.columns:
        
        if verbose:
            print target
        
        # get target values
        y = Y[target].values
        
        regGBR[target] = []
        test_score[target] = []
        
        for fold in range(n_folds):
            if verbose:
                print 'cv fold ', fold
            X_train, y_train = X.ix[kf[fold][0],:], y[kf[fold][0]]
            X_test, y_test = X.ix[kf[fold][1],:], y[kf[fold][1]]
        
            regGBR[target].append(GradientBoostingRegressor(n_estimators=n_estimators,
                                                            learning_rate=learning_rate,
                                                            max_depth=max_depth,
                                                            loss='ls'))
            regGBR[target][fold].fit(X_train, y_train)
        
            test_score[target].append(np.zeros((n_estimators,), dtype=np.float64))
            
            for i, y_pred in enumerate(regGBR[target][fold].staged_decision_function(X_test)):
                test_score[target][fold][i] = regGBR[target][fold].loss_(y_test, y_pred)
        
        #mse = mean_squared_error(y_test, regGBR[target].predict(logB_test))
        
    return regGBR, test_score


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
