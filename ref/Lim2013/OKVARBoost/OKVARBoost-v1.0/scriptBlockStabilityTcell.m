%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Script for computing Block-Stability on the T-cell data set 
%
% AUTHOR   : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version   : 2012/08/24
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% Inputs
load('DATA\TCELL\samples.mat');
dataTcell=sample;
M=20;
randFrac=4/58;
lambda2Vec=10.^(-1:1);
lambda1Vec=10.^(-2:0);
Bsize=7;
nBS=10;

%% Compute Block-Stability
[distJ] = blockStabilityOKVARboost(dataTcell,M,randFrac,lambda2Vec,lambda1Vec,Bsize,nBS);