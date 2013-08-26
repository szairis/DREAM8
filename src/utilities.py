import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pydot as dot

from copy import deepcopy


def plot_adjacency_matrix(A, title=None):
    '''Plots an adjacency matrix using pcolor. Input should be a dataframe with
    index and column names corresponding to nodes'''
    
    fig = plt.figure()

    plt.pcolor(A)
    plt.yticks(np.arange(0.5, len(A.index), 1), A.index)
    plt.xticks(np.arange(0.5, len(A.columns), 1), A.columns, rotation=90)
    plt.colorbar()
    plt.xlabel('Child')
    plt.ylabel('Parent')
    plt.xlim((0, len(A.index)))
    plt.ylim((0, len(A.columns)))

    if title is not None:
        plt.title(title)

    return fig


def plot_network(edge_weights_df, int_signs_df, graph_label, threshold):
    '''Generate a graphviz object from input posterior edge probabilities as 
    provided from the DBN code of Hill. This is also general enough to work
    with the output of any other model, just pass the edge_weights in twice.
    
    Inputs:
    edge_weight_df : square dataframe of edge weights
    int_signs_df : square dataframe of interaction signs (+/- 1);
    graph_label : (str) title of graph
    threshold : (float) threshold to use for drawing edges

    Output
    graph : dot object, use Ipython.display.Image(graph) for inline viz, or
            save file using graph.write_pdf(filename) if working in console
    '''

    node_list = edge_weights_df.columns
    nodes = dict((node, dot.Node(node, style='filled', fontcolor='blue')) for node in node_list)

    graph = dot.Dot(graph_type='digraph',
                label=graph_label, 
                fontsize=24,
                fontname='Helvetica',
                labelloc='t')
    
    map(graph.add_node, nodes.values())
    
    for node_p in node_list:
        for node_c in node_list:
            if edge_weights_df.ix[node_p, node_c] > threshold:
                if int_signs_df.ix[node_p, node_c] > 0:
                    ic = 'green'
                    at = 'normal'
                else:
                    ic = 'red'
                    at = 'tee'
                graph.add_edge(dot.Edge(nodes[node_p], nodes[node_c], 
                                        penwidth=edge_weights_df.ix[node_p,node_c]*5, 
                                        arrowhead=at,
                                        color=ic,
                                        label=edge_weights_df.ix[node_p, node_c],
                                        fontcolor='#0000ff'))
    
    return graph


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


def write_MIDAS(pred_dict, inhib_targets, path, datatype='inSilico', cell_line='inSilico'):
    f_midas = open(path + '.csv', 'w')
   
    print datatype
    if datatype is 'inSilico':
        stims = ['loLIG1', 'hiLIG1', 'loLIG2', 'hiLIG2', 'loLIG1_loLIG2', 'loLIG1_hiLIG2', 'hiLIG1_loLIG2', 'hiLIG1_hiLIG2']
        antibodies = list(pred_dict['loLIG1'].columns)
        time_points = [0,1,2,4,6,10,15,30,60,120]

    elif datatype is 'Experimental':
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
                    stim==stims[3], stim==stims[4], stim==stims[5], stim==stims[6], stim==stims[7]]))
            line += '{0},{1},{2},{3},{4},{5},{6},{7}'.format(*indvec)
            line += ',{0}'.format(tp)
            for ab in antibodies:
                if ab in inhib_targets:
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


def introduce_inhibs(data, inhib_targets, perfect=True):

    inhibs = inhib_targets.keys()

    inhib_df = data.copy()

    for inhib in inhibs:
        inhib_col = 'Inhib_{0}'.format(inhib)
        inhib_df[inhib_col] = 0
        inhib_df.ix[inhib_df['Inhibitor']==inhib, inhib_col] = 1
        if perfect:
            inhib_df.ix[inhib_df['Inhibitor']==inhib, inhib_targets[inhib]] = 0

    return inhib_df


def prepare_markov_data(data, response_type, group_stimuli):
    '''prepare_markov data
    '''

    training_dict = {}

    time_points = np.array(np.sort(list(set(data['Timepoint']))))
    t_min = time_points[0]
    t_max = time_points[-1]

    covariates = data.columns[4:]
    antibodies = [col for col in covariates if 'Inhib_' not in col]
    inhibs = set(data['Inhibitor'])
    stims = set(data['Stimulus'])

    num_conditions = len(data.groupby(['Inhibitor', 'Stimulus']).mean())
    num_antibodies = len(antibodies)


    if group_stimuli:
        #we need a training_dict with 1 key: 'all_stimuli'
        if response_type is 'level':
            # (level, True)
            X = data[data['Timepoint'] < t_max].groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().values
            Y = data[data['Timepoint'] > 0].groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix[:, antibodies].values
        else:
            # (rate, True)
            X = data[data['Timepoint'] < t_max].groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().values
            yf = data[data['Timepoint'] > 0].groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix[:, antibodies].values
            yb = data[data['Timepoint'] < t_max].groupby(['Inhibitor', 'Stimulus', 'Timepoint']).mean().ix[:, antibodies].values
            tf = time_points[1:]
            tb = time_points[:-1]
            dt = np.tile((tf-tb)[:, None], (num_conditions, num_antibodies))
            Y = (yf - yb) / dt
        training_dict['all_stimuli'] = (pd.DataFrame(X, columns=covariates), pd.DataFrame(Y, columns=antibodies))
    else:
        #we need a training_dict with 8 keys: the 8 different stimuli
        for stim in stims:
            if response_type is 'level':
                # (level, False)
                X = data[(data['Stimulus']==stim) & (data['Timepoint'] < t_max)].groupby(['Inhibitor', 'Timepoint']).mean().values
                Y = data[(data['Stimulus']==stim) & (data['Timepoint'] > 0)].groupby(['Inhibitor', 'Timepoint']).mean().ix[:, antibodies].values
            else:
                # (rate, False)
                X =  data[(data['Timepoint'] < t_max) & (data['Stimulus']==stim)].groupby(['Inhibitor', 'Timepoint']).mean().values
                yf = data[(data['Timepoint'] > 0) & (data['Stimulus']==stim)].groupby(['Inhibitor', 'Timepoint']).mean().ix[:, antibodies].values
                yb = data[(data['Timepoint'] < t_max) & (data['Stimulus']==stim)].groupby(['Inhibitor', 'Timepoint']).mean().ix[:, antibodies].values
                tf = time_points[1:]
                tb = time_points[:-1]
                dt = np.tile((tf-tb)[:, None], (len(set(data[data['Stimulus']==stim]['Inhibitor'])), num_antibodies))
                Y = (yf - yb) / dt
            training_dict[stim] = (pd.DataFrame(X, columns=covariates), pd.DataFrame(Y, columns=antibodies))
    
    return training_dict
