import numpy as np
import pandas as pd
import pydot as dot

def plot_network(edge_weights_df, int_signs_df, graph_label, threshold):
    '''Generate a graphviz object from input posterior edge probabilities as 
    provided from the DBN code of Hill.
    
    Inputs:
    edge_weight_df : square dataframe of edge weights
    int_signs_df : square dataframe of interaction signs (+/- 1);
    graph_label : (str) title of graph
    threshold : (float) threshold to use for drawing edges

    Output
    graph : dot object, use Ipython.display.Image(graph) for inline viz
    '''

    node_list = edge_weight_df.columns
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
