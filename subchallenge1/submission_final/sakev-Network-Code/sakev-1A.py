import os
import pandas as pd
import numpy as np
import utilities
from sklearn import preprocessing
from sklearn.ensemble import GradientBoostingRegressor


def build_adj_matrix(regGBR, node_list, stims):   
    num_nodes = len(node_list)
    adj_dict = {}
    for stim in stims:
        adj = np.zeros((num_nodes, num_nodes), dtype='f')
        for nidx, node in enumerate(node_list):
            adj[:, nidx] = regGBR[stim][node].feature_importances_[:num_nodes]
        adj_dict[stim] = pd.DataFrame(adj, index=node_list, columns=node_list)
        
    return adj_dict


def do_gbr(X, Y, n_estimators=100, learning_rate=0.1, max_depth=5, verbose=False):

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

if not os.path.exists('./results'):
    os.makedirs('./results')

bt20 = pd.read_csv('data/BT20_main.csv')
bt549 = pd.read_csv('data/BT549_main.csv')
mcf7 = pd.read_csv('data/MCF7_main.csv')
uacc812 = pd.read_csv('data/UACC812_main.csv')

prior = pd.read_csv('data/experimental_prior.csv', index_col=0, header=0)

## BT20
##############################################
print '----------- ' + 'BT20' + ' ------------'
inhibs = set(bt20['Inhibitor'])
stims = set(bt20['Stimulus'])
node_list = bt20.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}
regGBR= {}
scalar = {}
td_bt20 = utilities.prepare_markov_data(utilities.introduce_inhibs(bt20, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in td_bt20:
    X, Y = td_bt20[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3
    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)

adj = build_adj_matrix(regGBR, node_list, stims)
A_true_bt20 = prior
for stim in adj:
    path = 'results/sakev-BT20-{0}-Network'.format(stim)
    for row in A_true_bt20.index:
        for col in A_true_bt20.columns:
            if A_true_bt20.ix[row,col] == 1:
                adj[stim].ix[row,col] = 1
    utilities.write_SIF_EDA(adj[stim], path)

## BT549
##############################################
print '----------- ' + 'BT549' + ' ------------'
inhibs = set(bt549['Inhibitor'])
stims = set(bt549['Stimulus'])
node_list = bt549.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}
regGBR= {}
scalar = {}
td_bt549 = utilities.prepare_markov_data(utilities.introduce_inhibs(bt549, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in td_bt549:
    X, Y = td_bt549[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3
    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)

adj = build_adj_matrix(regGBR, node_list, stims)
A_true_bt549 = A_true_bt20.ix[td_bt549['EGF'][0].columns, td_bt549['EGF'][0].columns]
for stim in adj:
    path = 'results/sakev-BT549-{0}-Network'.format(stim)
    for row in A_true_bt549.index:
        for col in A_true_bt549.columns:
            if A_true_bt549.ix[row,col] == 1:
                adj[stim].ix[row,col] = 1
    utilities.write_SIF_EDA(adj[stim], path)

## MCF7
##############################################
print '----------- ' + 'MCF7' + ' ------------'
inhibs = set(mcf7['Inhibitor'])
stims = set(mcf7['Stimulus'])
node_list = mcf7.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}
regGBR= {}
scalar = {}
td_mcf7 = utilities.prepare_markov_data(utilities.introduce_inhibs(mcf7, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in td_mcf7:
    X, Y = td_mcf7[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3
    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)

adj = build_adj_matrix(regGBR, node_list, stims)
A_true_mcf7 = A_true_bt20.ix[td_mcf7['EGF'][0].columns, td_mcf7['EGF'][0].columns]
for stim in adj:
    path = 'results/sakev-MCF7-{0}-Network'.format(stim)
    for row in A_true_mcf7.index:
        for col in A_true_mcf7.columns:
            if A_true_mcf7.ix[row,col] == 1:
                adj[stim].ix[row,col] = 1
    utilities.write_SIF_EDA(adj[stim], path)

## UACC812
##############################################
print '----------- ' + 'UACC812' + ' ------------'
inhibs = set(uacc812['Inhibitor'])
stims = set(uacc812['Stimulus'])
node_list = uacc812.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}
regGBR= {}
scalar = {}
td_uacc812 = utilities.prepare_markov_data(utilities.introduce_inhibs(uacc812, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in td_uacc812:
    X, Y = td_uacc812[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3
    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)

adj = build_adj_matrix(regGBR, node_list, stims)
A_true_uacc812 = A_true_bt20.ix[td_uacc812['EGF'][0].columns, td_uacc812['EGF'][0].columns]
for stim in adj:
    path = 'results/sakev-UACC812-{0}-Network'.format(stim)
    for row in A_true_uacc812.index:
        for col in A_true_uacc812.columns:
            if A_true_uacc812.ix[row,col] == 1:
                adj[stim].ix[row,col] = 1
    utilities.write_SIF_EDA(adj[stim], path)
