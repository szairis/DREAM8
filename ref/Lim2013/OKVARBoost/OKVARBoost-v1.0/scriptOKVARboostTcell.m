%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Script for running OKVAR-Boost on the T-cell data set 
%
% AUTHOR   : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version   : 2012/08/24
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% OKVAR-Boost inputs
load('DATA\TCELL\samples.mat');
dataTcell=sample;
nRun=10;
M=20;
randFrac=6/58;
lambda2=1;
lambda1=1;

%% Run OKVAR-Boost
[Adj_thresh]=runMtsOKVARboost(dataTcell,nRun,M,randFrac,lambda2,lambda1);