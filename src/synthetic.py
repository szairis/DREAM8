import numpy as np
import pandas as pd
from scipy.integrate import odeint
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


def gen_xu(t_max, snr_meas=None, EGF=0, Gap=2400, Cilostamide=0, EPACA=0, PKAA=0):
    '''generate data from xu model

    Inputs:
    time_points: (list). list of times (in minutes) to sample data from
    snr_meas: (float). signal to noise ratio
    
    initial conditions (as keyword arguments):
    EGF
    Gap
    Cilostamide
    EPACA
    PKAA

    Output:
    A: true adjacency matrix
    dat: time series of network nodes
    '''
    
    # global parameters
    # should pass this as an additional argument to f(t, y)
    p1 = 0.00389106
    p2 = 8800.83
    p3 = 12.8748
    p4 = 6017.33
    p5 = 5.43844
    p6 = 7686.8
    p7 = 8378.11
    p8 = 2123.21
    p9 = 0.0889915
    p10 = 5422.75
    p11 = 25.0482
    p12 = 4950.17
    p13 = 17.658
    p14 = 5671.34
    p15 = 4511.32
    p16 = 13119.7
    p17 = 7.19206
    p18 = 5007.53
    p19 = 6965.92
    p20 = 632.512 
    p21 = 0.648428
    p22 = 450.436
    p23 = 5.1648
    p24 = 12463.3
    p25 = 12.0102
    p26 = 12976.4
    p27 = 1.66052
    p28 = 3102.31
    p29 = 9681.2
    p30 = 13519
    p31 = 4.07534
    p32 = 8741.85
    p33 = 9.34508
    p34 = 9699.87
    p35 = 8320.35
    p36 = 872.757
    p37 = 0.017411
    p38 = 2013.96
    p39 = 12.4142
    p40 = 3120.93
    p41 = 17.2111
    p42 = 9665.27
    p43 = 571.284
    p44 = 7.59316
    p45 = 0.20545
    p46 = 685.67
    p47 = 2.94579
    p48 = 4991.84
    p49 = 2109.86
    p50 = 2.08238
    p51 = 14421.3
    p52 = 5.85846
    p53 = 2432.18
    p54 = 462.007
    p55 = 193.304

    # solve the system dy/dt = f(y, t)
    def f(y, t):

        # fix negatives in current state
        y[y < 0] = 0

        unboundEGFR = y[0]
        boundEGFR = y[1]
        EGF = y[2]
        activeC3G = y[3]
        inactiveC3G = y[4]
        activeSOS = y[5]
        inactiveSOS = y[6]
        removedSOS = y[7]
        ERK = y[8]
        ERKPP = y[9]
        activeRas = y[10]
        inactiveRas = y[11]
        MEK = y[12]
        MEKPP = y[13]
        BRaf = y[14]
        BRafPP = y[15]
        Raf1 = y[16]
        Raf1PP = y[17]
        RemovedRaf1 = y[18] 
        PKA = y[19]
        inactivePKA = y[20]
        EPAC = y[21]
        inactiveEPAC = y[22]
        activeRap1 = y[23]
        inactiveRap1 = y[24]
        Gap = y[25]
        Cilostamide = y[26]
        EPACA = y[27]
        PKAA = y[28]

        
        # unboundEGFR
        f0 = -p5 * EGF * unboundEGFR + p6 * boundEGFR
        # boundEGFR
        f1 = p5 * EGF * unboundEGFR - p6 * boundEGFR
        # EGF
        f2 = -p5 * EGF * unboundEGFR + p6 * boundEGFR
        # activeC3G
        f3 = ((p47 * boundEGFR * inactiveC3G) / (p48 + inactiveC3G)) - ((p49 * activeC3G) / (p50 + activeC3G))
        # inactiveC3G
        f4 = -((p47 * boundEGFR * inactiveC3G)) / (p48 + inactiveC3G) + ((p49 * activeC3G) / (p50 + activeC3G))
        # activeSOS
        f5 = (p3*boundEGFR*inactiveSOS/(p4+inactiveSOS)) - (p8*activeSOS/(p7+activeSOS)) - (p1*ERKPP*activeSOS/(p2+activeSOS))
        # inactiveSOS
        f6 = -(p3*boundEGFR*inactiveSOS/(p4+inactiveSOS)) + (p8*activeSOS)/(p7+activeSOS) - (p1*ERKPP*inactiveSOS/(p2+inactiveSOS))
        # removed SOS
        f7 = (p1*ERKPP*inactiveSOS/(p2+inactiveSOS)) + (p1*ERKPP*activeSOS/(p2+activeSOS))
        # ERK
        f8 = -(p21*MEKPP*ERK/(p22+ERK)) + (p55*ERKPP/(p54+ERKPP))
        # ERKPP
        f9 = (p21*MEKPP*ERK/(p22+ERK)) - (p55*ERKPP/(p54+ERKPP))
        # activeRAS
        f10 = (p9*activeSOS*inactiveRas/(p10+inactiveRas)) - (p11*Gap*activeRas/(p12+activeRas))
        # inactiveRAS
        f11 = -(p9*activeSOS*inactiveRas/(p10+inactiveRas)) + (p11*Gap*activeRas/(p12+activeRas))
        # MEK
        f12 = -(p17*Raf1PP*MEK/(p18+MEK)) + (p20*MEKPP/(p19+MEKPP)) - (p45*BRafPP*MEK/(p46+MEK))
        # MEKPP
        f13 = (p17*Raf1PP*MEK/(p18+MEK)) - (p20*MEKPP/(p19+MEKPP)) + (p45*BRafPP*MEK/(p46+MEK))
        # BRaf
        f14 = -(p41*activeRap1*BRaf/(p42+BRaf)) + (p44*BRafPP/(p43+BRafPP)) - (p52*activeRas*BRaf/(p53+BRaf))
        # BRafPP
        f15 = (p41*activeRap1*BRaf/(p42+BRaf)) - (p44*BRafPP/(p43+BRafPP)) + (p52*activeRas*BRaf/(p53+BRaf))
        # Raf1
        f16 = -(p13*activeRas*Raf1/(p14+Raf1)) + (p16*Raf1PP/(p15+Raf1PP))
        # Raf1PP
        f17 = (p13*activeRas*Raf1/(p14+Raf1)) - (p16*Raf1PP/(p15+Raf1PP))
        # RemovedRaf1
        f18 = (p23*PKA*Raf1/(p24 + Raf1))
        # PKA
        f19 = (p25*PKAA*inactivePKA/(p26+inactivePKA)) + (p27*Cilostamide*inactivePKA/(p28+inactivePKA)) - (p30*PKA/(p29+PKA))
        # inactivePKA
        f20 = -(p25*PKAA*inactivePKA/(p26+inactivePKA)) - (p27*Cilostamide*inactivePKA/(p28+inactivePKA)) + (p30*PKA/(p29+PKA))
        # EPAC
        f21 = (p31*EPACA*inactiveEPAC/(p32+inactiveEPAC)) + (p33*Cilostamide*inactiveEPAC/(p34 + inactiveEPAC)) - (p36*EPAC/(p35 + EPAC))
        # inactiveEPAC
        f22 = -(p31*EPACA*inactiveEPAC/(p32+inactiveEPAC)) - (p33*Cilostamide*inactiveEPAC/(p34 + inactiveEPAC)) + (p36*EPAC/(p35 + EPAC))
        # activeRap1
        f23 = (p37*EPAC*inactiveRap1/(p38+inactiveRap1)) - (p39*Gap*activeRap1/(p40+activeRap1)) + (p50*activeC3G*inactiveRap1/(p51+inactiveRap1))
        # inactiveRap1
        f24 = -(p37*EPAC*inactiveRap1/(p38+inactiveRap1)) + (p39*Gap*activeRap1/(p40+activeRap1)) - (p50*activeC3G*inactiveRap1/(p51+inactiveRap1))
        # Gap
        f25 = 0
        # Cilostamide
        f26 = 0
        # EPACA
        f27 = 0
        # PKAA
        f28 = 0

        return [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,
                 f11, f12, f13, f14, f15, f16, f17, f18, f19,
                 f20, f21, f22, f23, f24, f25, f26, f27, f28]

    nodes = ['unboundEGFR','boundEGFR','EGF',
             'activeC3G','inactiveC3G',
             'activeSOS','inactiveSOS','removedSOS',
             'ERK','ERKPP',
             'activeRas','inactiveRas',
             'MEK','MEKPP',
             'BRaf','BRafPP',
             'Raf1','Raf1PP','RemovedRaf1',
             'PKA','inactivePKA',
             'EPAC','inactiveEPAC',
             'activeRap1','inactiveRap1',
             'Gap', 'Cilostamide', 'EPACA', 'PKAA']
    N = len(nodes)
    network_nodes = ['EGF', 'boundEGFR',
                     'activeSOS','activeC3G', 'activeRas',
                     'Raf1PP','MEKPP', 'ERKPP',
                     'PKA', 'EPAC', 'activeRap1', 'BRafPP',
                     'Gap', 'Cilostamide', 'EPACA', 'PKAA']
    NN = len(network_nodes)
    
    # true adjacency matrix
    # this is in a 'reduced form'
    A = pd.DataFrame(np.zeros((NN, NN)), index=network_nodes, columns=network_nodes)
    A.ix['EGF', 'boundEGFR'] = 1
    A.ix['boundEGFR', 'activeSOS'] = 1
    A.ix['boundEGFR', 'activeC3G'] = 1
    A.ix['activeSOS', 'activeRas'] = 1
    A.ix['activeRas', 'Raf1PP'] = 1
    A.ix['activeRas', 'BRafPP'] = 1
    A.ix['Raf1PP', 'MEKPP'] = 1
    A.ix['BRafPP', 'MEKPP'] = 1
    A.ix['MEKPP', 'ERKPP'] = 1
    A.ix['EPACA', 'EPAC'] = 1
    A.ix['PKAA', 'PKA'] = 1
    A.ix['Cilostamide', 'EPAC'] = 1
    A.ix['Cilostamide', 'PKA'] = 1
    A.ix['EPAC', 'activeRap1'] = 1
    A.ix['Gap', 'activeRap1'] = 1
    A.ix['Gap', 'BRafPP'] = 1
    A.ix['activeC3G', 'activeRap1'] = 1
    A.ix['activeRap1', 'BRafPP'] = 1
    A.ix['ERKPP', 'activeSOS'] = 1
    
    # initial conditions
    unboundEGFR0 = 500
    boundEGFR0 = 0
    EGF0 = EGF
    activeC3G0 = 0
    inactiveC3G0 = 1200
    activeSOS0 = 0
    inactiveSOS0 = 1200
    removedSOS0 = 0
    ERK0 = 10000
    ERKPP0 = 0
    activeRas0 = 0
    inactiveRas0 = 1200
    MEK0 = 3000
    MEKPP0 = 0
    Raf10 = 1500
    Raf1PP0 = 0
    removedRaf1 = 0
    BRaf0 = 1500
    BRafPP0 = 0
    PKA0 = 0
    inactivePKA0 = 1000
    EPAC0 = 0
    inactiveEPAC0 = 1000
    activeRap10 = 0
    inactiveRap10 = 1200

    y0 = [unboundEGFR0, boundEGFR0, EGF0,
          activeC3G0, inactiveC3G0,
          activeSOS0, inactiveSOS0, removedSOS0,
          ERK0, ERKPP0,
          activeRas0, inactiveRas0,
          MEK0, MEKPP0,
          BRaf0, BRafPP0,
          Raf10, Raf1PP0, removedRaf1,
          PKA0, inactivePKA0,
          EPAC0, inactiveEPAC0,
          activeRap10, inactiveRap10,
          Gap, Cilostamide, EPACA, PKAA]
    
    t_range = np.arange(0, t_max * 60 + 1, 1)
    soln = pd.DataFrame(odeint(f, y0, t_range), index=t_range, columns=nodes)
    
    # add measurement noise to solution
    if snr_meas:
        noise = soln.mean(axis=0) / snr_meas
        for tidx, time in enumerate(t_range):
            #pdb.set_trace()
            soln.ix[tidx, :] = np.random.multivariate_normal(soln.ix[tidx, :], noise.values**2 * np.eye(N))
        soln[soln < 0] = 0

    # we only take out the times in time_points
    #pdb.set_trace()
    dat = soln.ix[t_range, network_nodes]
    dat['Timepoint'] = t_range
    dat['Inhibitor'] = 'None'
    dat['Cell Line'] = 'Xu'
    dat['Stimulus'] = 'EGF_{0}_Gap_{1}_Cilostamide_{2}_EPACA_{3}_PKAA_{4}'.format(EGF,Gap,Cilostamide,EPACA,PKAA)
    
    return A, dat

