import numpy as np
import pandas as pd

#according to cokelaer we should omit unconnected nodes always
def write_SIF_EDA(A, path):
    f_sif = open(path + '.sif', 'w')
    f_eda = open(path + '.eda', 'w')
    f_eda.write('EdgeScore\n')
    prots = list(A.index)

    for parent in prots:
        if A.ix[parent,:].sum() == 0 and A.ix[:,parent].sum() == 0:
            continue
        for child in prots:
            if A.ix[parent,child] > 0:
                f_sif.write(parent + ' 1 ' + child + '\n')
                f_eda.write(parent + ' (1) ' + child + ' = {0}\n'.format(A.ix[parent,child]))
    f_sif.close()
    f_eda.close()


def write_MIDAS(pred, path):



def score_network(A):



def score_predictions(pred):

