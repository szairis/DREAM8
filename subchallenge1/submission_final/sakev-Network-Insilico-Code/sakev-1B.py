import os
import pandas as pd
import numpy as np
import utilities
from sklearn import preprocessing
from sklearn.ensemble import GradientBoostingRegressor

# function which builds the adjacency matrix from the feature rankings
def build_adj_matrix(regGBR, node_list, stims):   
    num_nodes = len(node_list)
    adj_dict = {}
    for stim in stims:
        adj = np.zeros((num_nodes, num_nodes), dtype='f')
        for nidx, node in enumerate(node_list):
            adj[:, nidx] = regGBR[stim][node].feature_importances_[:num_nodes]
        adj_dict[stim] = pd.DataFrame(adj, index=node_list, columns=node_list)
        
    return adj_dict

# function which fits ensembles of gradient boosted trees to training data (X, Y)
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

# create the results directory if it does no already exist
if not os.path.exists('./results'):
    os.makedirs('./results')

# load in the insilico data
insilico_data = pd.read_csv('data/insilico.csv', header=0)
inhibs = set(insilico_data['Inhibitor'])
stims = set(insilico_data['Stimulus'])

node_list = ['AB{0}'.format(i) for i in range(1, 21)]
inhib_targets = {'INH1' : 'AB12', 'INH2' : 'AB5', 'INH3' : 'AB8'}

regGBR= {}
scalar = {}
td = utilities.prepare_markov_data(utilities.introduce_inhibs(insilico_data, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=True)
X, Y = td['all_stimuli']
scalar['all_stimuli'] = preprocessing.StandardScaler()
scalar['all_stimuli'].fit_transform(X)
X.ix[X.ix[:,'Inhib_INH1']>0,'Inhib_INH1'] = 1
X.ix[X.ix[:,'Inhib_INH2']>0,'Inhib_INH2'] = 1
X.ix[X.ix[:,'Inhib_INH3']>0,'Inhib_INH3'] = 1
X.ix[X.ix[:,'Inhib_INH1']<0,'Inhib_INH1'] = 0
X.ix[X.ix[:,'Inhib_INH2']<0,'Inhib_INH2'] = 0
X.ix[X.ix[:,'Inhib_INH3']<0,'Inhib_INH3'] = 0

# Step 2 : Fit
n_estimators = 100
max_depth = 3

regGBR['all_stimuli'] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)

# Step 3 : build and write adjacency matrix
adj = build_adj_matrix(regGBR, node_list, ['all_stimuli'])
utilities.write_SIF_EDA(adj['all_stimuli'], 'results/sakev-Network-Insilico')
