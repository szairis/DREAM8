import numpy as np
import pandas as pd
from copy import deepcopy


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


def write_MIDAS(pred_dict, inhib_target, path, datatype='inSilico', cell_line='inSilico'):
    f_midas = open(path + '.csv', 'w')
    
    if datatype == 'inSilico':
        stims = ['loLIG1', 'hiLIG1', 'loLIG2', 'hiLIG2', 'loLIG1_loLIG2', 'loLIG1_hiLIG2', 'hiLIG1_loLIG2', 'hiLIG1_hiLIG2']
        antibodies = list(pred_dict['loLIG1'].columns)
        time_points = [0,1,2,4,6,10,15,30,45,60,120]

    elif dataype == 'Experimental':
        stims = ['Serum','PBS','EGF','Insulin','FGF1','HGF','NRG1','IGF1']
        antibodies = list(pred_dict['Serum'].columns)
        time_points = [0,5,15,30,60,120,240]
    
    header = 'TR:{0}:CellLine'.format(cell_line)
    for stim in stims:
        header += ',TR:{0}:Stimuli'.format(stim)
    header += ',DA:ALL'
    for ab in antibodies:
        header += ',DV:{0}'.format(ab)
    header += '\n'
    f_midas.write(header)

    for tp in time_points:
        for stim in stims:
            time_courses = pred_dict[stim]
            line = '1,'
            indvec = tuple(map(int, [stim==stims[0], stim==stims[1], stim==stims[2],
                    stim==stims[3], stims==stims[4], stim==stims[5], stim==stims[6], stim==stims[7]]))
            line += '{0},{1},{2},{3},{4},{5},{6},{7}'.format(*indvec)
            line += ',{0}'.format(tp)
            for ab in antibodies:
                if ab == inhib_target:
                    line += ',NA'
                else:
                    line += ',{0}'.format(time_courses.ix[tp,ab])
            line += '\n'
            f_midas.write(line)
    
    f_midas.close()


def im2bw(Adj, threshold):
    A = deepcopy(Adj)
    A[A < threshold] = 0
    A[A >= threshold] = 1
    return A


def score_network(Adj, Adj_true):
    A = np.abs(np.array(deepcopy(Adj), dtype='f'))
    A_true = np.abs(np.array(deepcopy(Adj_true), dtype='f'))
    n = len(A)
    A = np.reshape(A, (n**2,))
    idx = np.argsort(A)
    for i in range(1, n**2 + 1):
        A[idx[i-1]] = float((2*i - 1)) / float((2 * n**2))
    A = np.reshape(A, (n,n))
    
    tp = np.zeros((n**2 + 1))
    fp = np.zeros((n**2 + 1))
    tn = np.zeros((n**2 + 1))
    fn = np.zeros((n**2 + 1))
    
    for k in range(n**2 + 1):
        thresh = float(k) / float(n**2)
        inferred = im2bw(A, thresh)
#        print 'THRESHOLD: ' + str(thresh)
#        print 'INFERRED:'
#        print inferred
#        print 'TRUTH:'
#        print A_true
        true_pos = np.logical_and(inferred, A_true).astype('float')
        false_pos = inferred - true_pos
        true_neg = np.logical_not(np.logical_or(inferred, A_true)).astype('float')
        false_neg = np.logical_not(inferred).astype('float') - true_neg
        tp[k] = np.sum(true_pos)
        fp[k] = np.sum(false_pos)
        tn[k] = np.sum(true_neg)
        fn[k] = np.sum(false_neg)
#        print 'TP: {0}, FP: {1}, FN: {2}, TN: {3}'.format(tp[k], fp[k], fn[k], tn[k])
#        print 'FPR: {0}, TPR: {1}'.format(fp[k] / (fp[k] + tn[k]), tp[k] / (tp[k] + fn[k]))
#        print '-------------------------------------------'
        
    fpr = fp / (tn + fp)
    tpr = tp / (tp + fn)
#    plot(fpr, tpr)
#    show()
    auc = np.trapz(tpr[::-1], fpr[::-1])
    return auc


def score_predictions(pred_course, true_course):
    diff_log_squared = np.log(true_course / pred_course) ** 2
    rms_error = np.sqrt(np.array(diff_log_squared).mean())
    return rms_error

