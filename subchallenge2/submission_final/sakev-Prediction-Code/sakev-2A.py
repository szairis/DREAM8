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

# BT20
###########################################
cell_line = 'BT20'

test_inhib_targets = {'TestInhib1' : ['EGFR_pY1068', 'EGFR_pY1173', 'EGFR_pY992', 'HER2_pY1248'],
                        'TestInhib2' : ['Src_pY416','Src_pY527'],
                        'TestInhib3' : ['mTOR_pS2448'],
                        'TestInhib4' : ['EGFR_pY1068', 'EGFR_pY1173', 'EGFR_pY992'],
                        'TestInhib5' : []}

print '----------- ' + cell_line + ' ------------'
data = pd.read_csv('data/{0}_main.csv'.format(cell_line), header=0)
inhibs = set(data['Inhibitor'])
stims = set(data['Stimulus'])

node_list = data.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}

regGBR= {}
scalar = {}
td = utilities.prepare_markov_data(utilities.introduce_inhibs(data, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in stims:
    X, Y = td[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3

    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)


# predict
times = [0, 5, 15, 30, 60, 120, 240]
pred_dict = {}
for test_inhib in test_inhib_targets:
    print test_inhib
    pred_dict[test_inhib] = {}
    for stim in stims:
        scalarblah = scalar[stim]
        # set up new df to use, and fill t=0 values
        pred_df = pd.DataFrame(np.zeros((len(times), len(node_list))), index=times, columns=node_list)
        pred_df.ix[0, :] = data.groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix['DMSO', stim, 0]
        pred_df.ix[0, test_inhib_targets[test_inhib]] = 0
        
        # loop over times
        for tidx in range(1,len(times)):
            time = times[tidx]
            
            # get covariates for this time step and scale
            # covariates_df = scaler.transform(pred_df.ix[times[tidx-1], :])
            covariates_df = ((pred_df.ix[times[tidx-1], :]) - scalarblah.mean_[:-2]) / scalarblah.std_[:-2]

            # zero out covariate we are inhibiting
            try:
                covariates_df.ix[test_inhib_targets[test_inhib]] = 0
            except:
                pass
            
            num_cov = len(pred_df.columns)
            covariates = np.zeros((num_cov + 2,))
            covariates[:num_cov] = covariates_df.values

            # loop over proteins to get values for current time step
            for p in node_list:
                pred_df.ix[time, p] = regGBR[stim][p].predict(covariates)
            
            # zero out covariate we are inhibiting, again
            pred_df.ix[time, test_inhib_targets[test_inhib]] = 0
    
        pred_dict[test_inhib][stim] = pred_df
        
    out_path = 'results/sakev-{0}-{1}-Prediction'.format(cell_line, test_inhib)
    utilities.write_MIDAS(pred_dict[test_inhib], test_inhib_targets[test_inhib], out_path, datatype='Experimental', cell_line=cell_line)

# BT549
###########################################
cell_line = 'BT549'

test_inhib_targets = {'TestInhib1' : ['EGFR_pY1068', 'EGFR_pY1173', 'HER2_pY1248'],
                        'TestInhib2' : ['Src_pY416','Src_pY527'],
                        'TestInhib3' : ['mTOR_pS2448'],
                        'TestInhib4' : ['EGFR_pY1068', 'EGFR_pY1173'],
                        'TestInhib5' : []}

print '----------- ' + cell_line + ' ------------'
data = pd.read_csv('data/{0}_main.csv'.format(cell_line), header=0)
inhibs = set(data['Inhibitor'])
stims = set(data['Stimulus'])

node_list = data.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}

regGBR= {}
scalar = {}
td = utilities.prepare_markov_data(utilities.introduce_inhibs(data, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in stims:
    X, Y = td[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3

    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)


# predict
times = [0, 5, 15, 30, 60, 120, 240]
pred_dict = {}
for test_inhib in test_inhib_targets:
    print test_inhib
    pred_dict[test_inhib] = {}
    for stim in stims:
        scalarblah = scalar[stim]
        # set up new df to use, and fill t=0 values
        pred_df = pd.DataFrame(np.zeros((len(times), len(node_list))), index=times, columns=node_list)
        pred_df.ix[0, :] = data.groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix['DMSO', stim, 0]
        pred_df.ix[0, test_inhib_targets[test_inhib]] = 0
        
        # loop over times
        for tidx in range(1,len(times)):
            time = times[tidx]
            
            # get covariates for this time step and scale
            # covariates_df = scaler.transform(pred_df.ix[times[tidx-1], :])
            covariates_df = ((pred_df.ix[times[tidx-1], :]) - scalarblah.mean_[:-2]) / scalarblah.std_[:-2]

            # zero out covariate we are inhibiting
            try:
                covariates_df.ix[test_inhib_targets[test_inhib]] = 0
            except:
                pass
            
            num_cov = len(pred_df.columns)
            covariates = np.zeros((num_cov + 2,))
            covariates[:num_cov] = covariates_df.values

            # loop over proteins to get values for current time step
            for p in node_list:
                pred_df.ix[time, p] = regGBR[stim][p].predict(covariates)
            
            # zero out covariate we are inhibiting, again
            pred_df.ix[time, test_inhib_targets[test_inhib]] = 0
    
        pred_dict[test_inhib][stim] = pred_df
        
    out_path = 'results/sakev-{0}-{1}-Prediction'.format(cell_line, test_inhib)
    utilities.write_MIDAS(pred_dict[test_inhib], test_inhib_targets[test_inhib], out_path, datatype='Experimental', cell_line=cell_line)

# MCF7
###########################################
cell_line = 'MCF7'

test_inhib_targets = {'TestInhib1' : ['EGFR_pY1068', 'EGFR_pY1173', 'EGFR_pY992', 'HER2_pY1248'],
                        'TestInhib2' : ['Src_pY416','Src_pY527'],
                        'TestInhib3' : ['mTOR_pS2448'],
                        'TestInhib4' : ['EGFR_pY1068', 'EGFR_pY1173', 'EGFR_pY992'],
                        'TestInhib5' : []}

print '----------- ' + cell_line + ' ------------'
data = pd.read_csv('data/{0}_main.csv'.format(cell_line), header=0)
inhibs = set(data['Inhibitor'])
stims = set(data['Stimulus'])

node_list = data.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}

regGBR= {}
scalar = {}
td = utilities.prepare_markov_data(utilities.introduce_inhibs(data, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in stims:
    X, Y = td[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3

    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)


# predict
times = [0, 5, 15, 30, 60, 120, 240]
pred_dict = {}
for test_inhib in test_inhib_targets:
    print test_inhib
    pred_dict[test_inhib] = {}
    for stim in stims:
        scalarblah = scalar[stim]
        # set up new df to use, and fill t=0 values
        pred_df = pd.DataFrame(np.zeros((len(times), len(node_list))), index=times, columns=node_list)
        pred_df.ix[0, :] = data.groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix['DMSO', stim, 0]
        pred_df.ix[0, test_inhib_targets[test_inhib]] = 0
        
        # loop over times
        for tidx in range(1,len(times)):
            time = times[tidx]
            
            # get covariates for this time step and scale
            # covariates_df = scaler.transform(pred_df.ix[times[tidx-1], :])
            covariates_df = ((pred_df.ix[times[tidx-1], :]) - scalarblah.mean_[:-2]) / scalarblah.std_[:-2]

            # zero out covariate we are inhibiting
            try:
                covariates_df.ix[test_inhib_targets[test_inhib]] = 0
            except:
                pass
            
            num_cov = len(pred_df.columns)
            covariates = np.zeros((num_cov + 2,))
            covariates[:num_cov] = covariates_df.values

            # loop over proteins to get values for current time step
            for p in node_list:
                pred_df.ix[time, p] = regGBR[stim][p].predict(covariates)
            
            # zero out covariate we are inhibiting, again
            pred_df.ix[time, test_inhib_targets[test_inhib]] = 0
    
        pred_dict[test_inhib][stim] = pred_df
        
    out_path = 'results/sakev-{0}-{1}-Prediction'.format(cell_line, test_inhib)
    utilities.write_MIDAS(pred_dict[test_inhib], test_inhib_targets[test_inhib], out_path, datatype='Experimental', cell_line=cell_line)

# UACC812
###########################################
cell_line = 'UACC812'

test_inhib_targets = {'TestInhib1' : ['EGFR_pY1068', 'EGFR_pY1173', 'HER2_pY1248'],
                        'TestInhib2' : ['Src_pY416','Src_pY527'],
                        'TestInhib3' : ['mTOR_pS2448'],
                        'TestInhib4' : ['EGFR_pY1068', 'EGFR_pY1173'],
                        'TestInhib5' : []}

print '----------- ' + cell_line + ' ------------'
data = pd.read_csv('data/{0}_main.csv'.format(cell_line), header=0)
inhibs = set(data['Inhibitor'])
stims = set(data['Stimulus'])

node_list = data.columns[4:]
inhib_targets = {'GSK690693' : ['AKT_pT308','AKT_pS473'],
                 'GSK690693_GSK1120212' : ['AKT_pT308','AKT_pS473','MEK1_pS217_S221']}

regGBR= {}
scalar = {}
td = utilities.prepare_markov_data(utilities.introduce_inhibs(data, inhib_targets=inhib_targets, perfect=True), 'level', group_stimuli=False)
for stim in stims:
    X, Y = td[stim]
    scalar[stim] = preprocessing.StandardScaler()
    scalar[stim].fit_transform(X)
    X.ix[X.ix[:,'Inhib_GSK690693']>0,'Inhib_GSK690693'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']>0,'Inhib_GSK690693_GSK1120212'] = 1
    X.ix[X.ix[:,'Inhib_GSK690693']<0,'Inhib_GSK690693'] = 0
    X.ix[X.ix[:,'Inhib_GSK690693_GSK1120212']<0,'Inhib_GSK690693_GSK1120212'] = 0
    n_estimators = 100
    max_depth = 3

    regGBR[stim] = do_gbr(X, Y, n_estimators=n_estimators, max_depth=max_depth)


# predict
times = [0, 5, 15, 30, 60, 120, 240]
pred_dict = {}
for test_inhib in test_inhib_targets:
    print test_inhib
    pred_dict[test_inhib] = {}
    for stim in stims:
        scalarblah = scalar[stim]
        # set up new df to use, and fill t=0 values
        pred_df = pd.DataFrame(np.zeros((len(times), len(node_list))), index=times, columns=node_list)
        pred_df.ix[0, :] = data.groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix['DMSO', stim, 0]
        pred_df.ix[0, test_inhib_targets[test_inhib]] = 0
        
        # loop over times
        for tidx in range(1,len(times)):
            time = times[tidx]
            
            # get covariates for this time step and scale
            # covariates_df = scaler.transform(pred_df.ix[times[tidx-1], :])
            covariates_df = ((pred_df.ix[times[tidx-1], :]) - scalarblah.mean_[:-2]) / scalarblah.std_[:-2]

            # zero out covariate we are inhibiting
            try:
                covariates_df.ix[test_inhib_targets[test_inhib]] = 0
            except:
                pass
            
            num_cov = len(pred_df.columns)
            covariates = np.zeros((num_cov + 2,))
            covariates[:num_cov] = covariates_df.values

            # loop over proteins to get values for current time step
            for p in node_list:
                pred_df.ix[time, p] = regGBR[stim][p].predict(covariates)
            
            # zero out covariate we are inhibiting, again
            pred_df.ix[time, test_inhib_targets[test_inhib]] = 0
    
        pred_dict[test_inhib][stim] = pred_df
        
    out_path = 'results/sakev-{0}-{1}-Prediction'.format(cell_line, test_inhib)
    utilities.write_MIDAS(pred_dict[test_inhib], test_inhib_targets[test_inhib], out_path, datatype='Experimental', cell_line=cell_line)
