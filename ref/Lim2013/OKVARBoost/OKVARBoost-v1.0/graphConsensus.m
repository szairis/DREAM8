function [Gcons] = graphConsensus(listG,freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Get consensus network from a committee of graphs for different
% edge frequency thresholds
%
% USAGE :
% [Gcons] = graphConsensus(listG,freq)
%
% INPUTS :
% listG     : ([1,n] cell) list of adjacency matrices
% freq      : ([nFreq] vector) edge frequency threshold for consensus 
%
% OUTPUT :
% Gcons     : ([1,nFreq] cell) adjacency matrices for each edge frequency 
%             threshold
%
% AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version : 2012/02/22
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nGraph=length(listG);
    nFreq=length(freq);
    Gcons=cell(1,nFreq);
    tmp=zeros(length(listG{1}));
    
    for j=1:length(listG)
        tmp = tmp+abs(listG{j});
    end
    
    for i=1:nFreq
        Gcons{i} = (tmp >= freq(i)*nGraph);
    end
    
end
