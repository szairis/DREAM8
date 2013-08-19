import numpy as np
import pandas as pd
from scipy.integrate import odeint
import pdb

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
def f(y, t, Gap, Cilostamide, EPACA, PKAA):

    # current state
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
    
    return [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,
             f11, f12, f13, f14, f15, f16, f17, f18, f19,
             f20, f21, f22, f23, f24]


def gen_xu(time_points, snr_meas=None, EGF=10, Gap=2400, Cilostamide=0, EPACA=0, PKAA=0):
    
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
             'activeRap1','inactiveRap1']
    N = len(nodes)
    
    network_nodes = ['EGF', 'EGFR',
                     'SOS', 'C3G', 'Ras', 'Raf1',
                     'MEK', 'ERK',
                     'PKA', 'EPAC', 'Rap1', 'BRaf',
                     'Gap', 'Cilostamide', 'EPACA', 'PKAA']
    NN = len(network_nodes)
    
    # true adjacency matrix
    # this is in a 'reduced form'
    A = pd.DataFrame(np.zeros((NN, NN)), index=network_nodes, columns=network_nodes)
    A.ix['EGF', 'EGFR'] = 1
    A.ix['EGFR', 'SOS'] = 1
    A.ix['EGFR', 'C3G'] = 1
    A.ix['SOS', 'Ras'] = 1
    A.ix['Ras', 'Raf1'] = 1
    A.ix['Ras', 'BRaf'] = 1
    A.ix['Raf1', 'MEK'] = 1
    A.ix['BRaf', 'MEK'] = 1
    A.ix['MEK', 'ERK'] = 1
    A.ix['EPACA', 'EPAC'] = 1
    A.ix['PKAA', 'PKA'] = 1
    A.ix['Cilostamide', 'EPAC'] = 1
    A.ix['Cilostamide', 'PKA'] = 1
    A.ix['EPAC', 'Rap1'] = 1
    A.ix['Gap', 'Rap1'] = 1
    A.ix['Gap', 'BRaf'] = 1
    A.ix['C3G', 'Rap1'] = 1
    A.ix['Rap1', 'BRaf'] = 1

    A.ix['ERK', 'SOS'] = 1
    
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
          activeRap10, inactiveRap10]
    
    # agonists
    stim = (Gap, Cilostamide, EPACA, PKAA)
    soln = odeint(f, y0, time_points, args=stim)
    
    # add measurement noise to solution
    if snr_meas:
        noise = soln.mean(axis=0) / snr_meas
        for tidx, time in enumerate(time_points):
            soln[tidx, :] = np.random.multivariate_normal(soln[tidx, :], noise**2 * np.eye(N))

    dat = pd.DataFrame(soln, index=time_points, columns=nodes)
    
    return A, dat
