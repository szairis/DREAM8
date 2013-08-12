import numpy as np
import pandas as pd
#import pdb

store = pd.HDFStore('pandas/insilico_panels', mode='w')

fname = 'insilico/CSV/insilico.csv'
f = open(fname)

nodes = f.readline().strip().split(',')[4:]
print nodes
inhibs = ['No_INH','INH1','INH2','INH3']
lig1 = ['no_lig1','lo_lig1','hi_lig1']
lig2 = ['no_lig2','lo_lig2','hi_lig2']
timepts = [0,1,2,4,6,10,15,30,45,60,120]

p4d_none = pd.Panel4D(np.zeros((3, 3, 11, len(nodes))))
p4d_inh1 = pd.Panel4D(np.zeros((3, 3, 11, len(nodes))))
p4d_inh2 = pd.Panel4D(np.zeros((3, 3, 11, len(nodes))))
p4d_inh3 = pd.Panel4D(np.zeros((3, 3, 11, len(nodes))))
p4d_none.labels = lig1
p4d_none.items = lig2
p4d_none.major_axis = timepts
p4d_none.minor_axis = nodes
p4d_inh1.labels = lig1
p4d_inh1.items = lig2
p4d_inh1.major_axis = timepts
p4d_inh1.minor_axis = nodes
p4d_inh2.labels = lig1
p4d_inh2.items = lig2
p4d_inh2.major_axis = timepts
p4d_inh2.minor_axis = nodes
p4d_inh3.labels = lig1
p4d_inh3.items = lig2
p4d_inh3.major_axis = timepts
p4d_inh3.minor_axis = nodes


time_conversions = {'0min':0, '1min':1, '2min':2, '4min':4, '6min':6, '10min':10, 
                    '15min':15, '30min':30, '45min':45, '60min':60, '120min':120}

lig_conversions = {'loLIG1':('lo_lig1','no_lig2'), 'hiLIG1':('hi_lig1','no_lig2'), 'loLIG2':('no_lig1','lo_lig2'),
                    'hiLIG2':('no_lig1','hi_lig2'), 'loLIG1+loLIG2':('lo_lig1','lo_lig2'), 'loLIG1+hiLIG2':('lo_lig1','hi_lig2'),
                    'hiLIG1+loLIG2':('hi_lig1','lo_lig2'), 'hiLIG1+hiLIG2':('hi_lig1','hi_lig2')}

for i in range(220):
    frags1 = f.readline().strip().split(',')
    frags2 = f.readline().strip().split(',')
    frags3 = f.readline().strip().split(',')
    arr1 = np.array(map(float,frags1[4:]))
    arr2 = np.array(map(float,frags2[4:]))
    arr3 = np.array(map(float,frags3[4:]))
    avg_arr = (arr1 + arr2 + arr3) / 3
    
    print frags1[1], lig_conversions[frags1[2]], time_conversions[frags1[3]]

    if frags1[1] == '':
        p4d_none.ix[lig_conversions[frags1[2]][0], lig_conversions[frags1[2]][1], time_conversions[frags1[3]], :] = avg_arr
    elif frags1[1] == 'INH1':
        p4d_inh1.ix[lig_conversions[frags1[2]][0], lig_conversions[frags1[2]][1], time_conversions[frags1[3]], :] = avg_arr
    elif frags1[1] == 'INH2':
        p4d_inh2.ix[lig_conversions[frags1[2]][0], lig_conversions[frags1[2]][1], time_conversions[frags1[3]], :] = avg_arr
    elif frags1[1] == 'INH3':
        p4d_inh3.ix[lig_conversions[frags1[2]][0], lig_conversions[frags1[2]][1], time_conversions[frags1[3]], :] = avg_arr
    

store.append('No_INH', p4d_none)
store.append('INH1', p4d_inh1)
store.append('INH2', p4d_inh2)
store.append('INH3', p4d_inh3)
f.close()
store.close()

