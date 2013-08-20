import pandas as pd
import numpy as np
import pdb

def drift(x, x_lag):
    
    # prevent blow-up
    x = np.maximum(x, np.zeros(x.shape))
    x_lag = np.maximum(x_lag, np.zeros(x_lag.shape))
    
    # parameters
    alpha = np.array([0, 0.000149, 0.003, 0.00074, 0.00061])
    v = np.array([0.04, 0.000882, 0.0201, 0.0147, 0.0182])
    k = np.array([1, 0.0356, 0.0372, 0.01, 1.814, 1.814])
    h = np.array([1, 1, 1, 4, 1, 1])
    d = np.array([0.0222, 0.0478, 0.4217, 0.098, 0.05])
    gamma = 0.6
    
    # drift function
    x_dot = np.zeros((5,))
    x_dot[0] = alpha[0] + v[0] * (x_lag[2]**h[0] / ((k[0]**h[0] + x_lag[2]**h[0])*(1 + (x[4]/k[1])**h[1]))) - d[0]*x[0]
    x_dot[1] = alpha[1] + v[1] * (x[0]**h[2] / (k[2]**h[2] + x[0]**h[2])) - d[1]*x[1]
    x_dot[2] = alpha[2] + v[2] * (x[1]**h[3] / (k[3]**h[3] + x[1]**h[3]*(1+(x[3]/gamma)**4))) - d[2]*x[2]
    x_dot[3] = alpha[3] + v[3] * (x[2]**h[4] / (k[4]**h[4] + x[2]**h[2])) - d[3]*x[3]
    x_dot[4] - alpha[4] + v[4] * (x[2]**h[5] / (k[5]**h[5] + x[2]**h[5])) - d[4]*x[4]
    
    return x_dot


def SDDE(drift, diffusion, lag, x_0, t):
    
    res = 1
    t = res * np.round(t/res)
    lag = res * np.round(lag/res)
    P = len(x_0)
    
    # initial condition
    x = np.zeros((P, lag+t))
    x[:, :lag] = np.tile(x_0, (lag, 1)).T
    
    # numerically solve
    for n in range(t):
        x[:, lag+n] = x[:, lag+n-1] + drift(x[:, lag+n-1], x[:, n]) * res
        x[:, lag+n] = np.random.multivariate_normal(x[:, lag+n].T, diffusion(x[:, lag+n-1])).T
    
    sol = x[:, lag+t-1]
    return sol


def SDDE_ave(t, snr_cell, mod_S):
    
    x_0 = np.array([0.0045, 0.0324, 0.0075, 0.0221, 0.012])

    # empirical distribution threshold
    mod_S_star = 30
    
    # multiplicative cell noise
    diffusion = lambda X: np.abs(np.diag(X)) / (1000000 * snr_cell)
    
    # destructively sample cells at time t
    P = x_0.shape[0]
    sample_size = min(mod_S, mod_S_star)
    sample = np.zeros((P, sample_size))
    lag = 100
    noise = np.mean(x_0) / snr_cell
    for cell in range(sample_size):
        initial = np.random.multivariate_normal(x_0, noise**2 * np.eye(P))
        sample[:, cell] = SDDE(drift, diffusion, lag, initial, t)
        
    # average mod_S cells
    y = np.zeros((P,))
    for cell in range(mod_S):
        y = y + sample[:, np.random.randint(sample_size)]
        
    y = y / mod_S
    
    return y
    

def gen_cantone(time_points, snr_cell=10, snr_meas=10, mod_S=30):
    '''
    Generate data according to the Cantone 5-gene network model.
    Originally written in MATLAB by Chris Oates.
    Adapted for Python by kje.
    
    Input:
    time_points : array of sample timepoints
    snr_cell : signal to noise ratio in the cells (default=10)
    snr_meas : signal to noise ratio in the measurement (default=10)
    mod_S : number of cells to average over (default=30)

    Output:
    A : true adjacency matrix of Cantone model (dataframe)
    dat : simulated time series (dataframe)
    '''
    
    genes = ['CBF1', 'GAL4', 'SWI5', 'GAL80', 'ASH1']

    # true adjacency matrix
    adj = np.array([[0, 1, 0, 0, 0],
                  [0, 0, 1, 0, 0],
                  [1, 0, 0, 1, 1],
                  [0, -1, 0, 0, 0],
                  [-1, 0, 0, 0, 0]])
    
    A = pd.DataFrame(adj, index=genes, columns=genes)
    
    
    # simulate the experiment
    P = A.shape[0]
    y = np.zeros((len(time_points), P))
    for tidx, time in enumerate(time_points):
        y[tidx, :] = SDDE_ave(time, snr_cell, mod_S)
    
    # add measurement error
    noise = np.mean(y) / snr_meas
    for tidx, time in enumerate(time_points):
        y[tidx, :] = np.random.multivariate_normal(y[tidx, :], noise**2 * np.eye(P))
    dat = pd.DataFrame(y, index=time_points, columns=genes)
    
    return A, dat
