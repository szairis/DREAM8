import numpy as np
import pandas as pd

cell_lines = ['BT20','BT549','MCF7','UACC812']
store = pd.HDFStore('pandas/experimental_panels', mode='w')

for cell_line in cell_lines:
    fname = 'experimental/CSV/' + cell_line + '_main.csv'
    f = open(fname)
    if cell_line == 'UACC812':
        f.readline(); f.readline(); f.readline(); f.readline()
    else:
        f.readline(); f.readline(); f.readline()
    nodes = f.readline().split(',')[4:]
    inhibs = ['GSK690693','GSK690693_GSK1120212','PD173074','DMSO']
    stims = ['FGF1','Insulin','EGF','IGF1','HGF','Serum','NRG1','PBS']
    timepts = [0,5,15,30,60,120,240]
    
    tmp = np.zeros((4,8,7,len(nodes)))
    p4d = pd.Panel4D(tmp)
    p4d.labels = inhibs
    p4d.items = stims
    p4d.major_axis = timepts
    p4d.minor_axis = nodes
    time_conversions = {'0min':0, '5min':5, '15min':15, '30min':30, '60min':60, '2hr':120, '4hr':240}
    
    ss_1 = np.zeros((len(nodes),))
    ss_1_counter = 0
    ss_2 = np.zeros((len(nodes),))
    ss_2_counter = 0
    ss_3 = np.zeros((len(nodes),))
    ss_3_counter = 0
    ss_4 = np.zeros((len(nodes),))
    ss_4_counter = 0
    
    for line in f:
        frags = line.split(',')
        if frags[3] != '0min':
            p4d.ix[frags[1], frags[2], time_conversions[frags[3]], nodes] = frags[4:]
        else:
            if frags[1] == 'GSK690693':
                ss_1 = ss_1 + np.array(map(float, frags[4:]))
                ss_1_counter += 1
            if frags[1] == 'GSK690693_GSK1120212':
                ss_2 = ss_2 + np.array(map(float, frags[4:]))
                ss_2_counter += 1
            if frags[1] == 'PD173074':
                ss_3 = ss_3 + np.array(map(float, frags[4:]))
                ss_3_counter += 1
            if frags[1] == 'DMSO':
                ss_4 = ss_4 + np.array(map(float, frags[4:]))
                ss_4_counter += 1
    
    ss_1 = ss_1 / ss_1_counter
    ss_2 = ss_2 / ss_2_counter
    ss_3 = ss_3 / ss_3_counter
    ss_4 = ss_4 / ss_4_counter
    
    for stim in stims:
        p4d.ix['GSK690693',stim,0,:] = ss_1
        p4d.ix['GSK690693_GSK1120212',stim,0,:] = ss_2
        p4d.ix['PD173074',stim,0,:] = ss_3
        p4d.ix['DMSO',stim,0,:] = ss_4
    
    store.append(cell_line, p4d)
    print cell_line.upper() + ':'
    for inhib in inhibs:
        for stim in stims:
            for time in timepts:
                if p4d.ix[inhib, stim, time, 'EIF4EBP1_pS65'] == 0:
                    print (inhib, stim, time)
    print '\n'
    f.close()

store.close()

