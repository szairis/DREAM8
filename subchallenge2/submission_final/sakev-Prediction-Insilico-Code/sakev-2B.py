import os
import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor
from sklearn import preprocessing
import utilities

### FUNCTIONS

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

# load data
insilico_data = pd.read_csv('data/insilico.csv', header=0)
inhibs = set(insilico_data['Inhibitor'])
stims = set(insilico_data['Stimulus'])

node_list = ['AB{0}'.format(i) for i in range(1, 21)]
inhib_targets = {'INH1' : 'AB12', 'INH2' : 'AB5', 'INH3' : 'AB8'}
INH_targets = ['_'.join(['Inhib', node]) for node in node_list]

# wrangle data
regGBR= {}
scalar = {}
td = utilities.prepare_markov_data(utilities.introduce_inhibs(insilico_data, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in stims:
    X, Y = td[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_INH1']>0,'Inhib_INH1'] = 1
    X.ix[X.ix[:,'Inhib_INH2']>0,'Inhib_INH2'] = 1
    X.ix[X.ix[:,'Inhib_INH3']>0,'Inhib_INH3'] = 1
    X.ix[X.ix[:,'Inhib_INH1']<0,'Inhib_INH1'] = 0
    X.ix[X.ix[:,'Inhib_INH2']<0,'Inhib_INH2'] = 0
    X.ix[X.ix[:,'Inhib_INH3']<0,'Inhib_INH3'] = 0

    # Step 2 : Fit
    n_estimators = 100
    max_depth = 3

    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)
    #plot_regGBR(X, regGBR, test_score, n_estimators=n_estimators)


# predict
times = [0, 1, 2, 4, 6, 10, 15, 30, 60, 120]
pred_dict = {}
for test_inhib in node_list:
    print test_inhib
    pred_dict[test_inhib] = {}
    for stim in stims:
        scalarblah = scalar[stim]
        # set up new df to use, and fill t=0 values
        pred_df = pd.DataFrame(np.zeros((len(times), len(node_list))), index=times, columns=node_list)
        pred_df.ix[0, :] = insilico_data.groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix['None', stim, 0]
        pred_df.ix[0, test_inhib] = 0
        
        # loop over times
        for tidx in range(1,len(times)):
            time = times[tidx]
            
            # get covariates for this time step and scale
            # covariates_df = scaler.transform(pred_df.ix[times[tidx-1], :])
            covariates_df = ((pred_df.ix[times[tidx-1], :]) - scalarblah.mean_[:-3]) / scalarblah.std_[:-3]

            # zero out covariate we are inhibiting
            try:
                covariates_df.ix[test_inhib_targets[test_inhib]] = 0
            except:
                pass
            
            num_cov = len(pred_df.columns)
            covariates = np.zeros((num_cov + 3,))
            covariates[:num_cov] = covariates_df.values

            # loop over proteins to get values for current time step
            for p in node_list:
                pred_df.ix[time, p] = regGBR[stim][p].predict(covariates)
            
            # zero out covariate we are inhibiting, again
            pred_df.ix[time, test_inhib] = 0
    
        pred_dict[test_inhib][stim] = pred_df
        
### SAVE TO MIDAS
for node in node_list:
    out_path = 'results/sakev-{0}-Prediction-Insilico'.format(node)
    utilities.write_MIDAS(pred_dict[node], node, out_path, datatype='inSilico', cell_line='inSilico')
