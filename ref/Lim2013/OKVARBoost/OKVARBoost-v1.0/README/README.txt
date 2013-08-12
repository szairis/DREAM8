 This file explains how to use OKVAR-Boost code for the inference of gene regulatory 
networks from time-course gene expression data.

Authors : Nehemy LIM (AMIS group, IBISC, dept of computer sciences, University of Evry, Evry, FRANCE), Yasin Senbabaoglu (Dept statistics, University of Michigan, Ann Arbor, USA)
Date: 11/6/2012
This codes comes with the paper entitled "OKVARBoost: a novel algorithm"



---------------------------------------------------------------------------



Follow the instructions below :



1) Import your data in the format of a cell array of matrices of size N * p where N is the 
number of time points and p is the number of genes.


2) See example in script scriptOKVARboostTcell.m to get consensus network from multiple runs and/or multiple time series.


3) See example in script scriptBlockStabilityTcell.m to compute Block-Stability from bootstrap samples.