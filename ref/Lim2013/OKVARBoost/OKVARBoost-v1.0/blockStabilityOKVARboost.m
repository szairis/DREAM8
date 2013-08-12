function [distJ] = blockStabilityOKVARboost(data,M,randFrac,lambda2Vec,lambda1Vec,Bsize,nBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Compute Block-instability criterion
%
% USAGE : [distJ] = blockStabilityOKVARboost(data,M,randFrac,lambda2Vec,lambda1Vec,Bsize,nBS)
%
% INPUTS
% data          : cell array, data{ts} contains time-series number 'ts' 
%                data{ts} is an [n,p] matrix (n time points/p variables,genes) 
% M             : maximum number of iterations in a single OKVAR-Boost run
% randFrac      : size of random subsets (in percentage of total number of genes)
% lambda2Vec    : vector (size nL2) of tested ridge parameters
% lambda1Vec    : vector (size nL1) of tested L1-regularization parameters
% Bsize         : size of one block bootstrap
% nBS           : number of pairs of block bootstraps
%
% OUTPUTS
% distJ         : [nL2,nL1] matrix, distJ(iL2,iL2) is the block-instability
%                 value for lambda2=lambda2Vec(iL2) and
%                 lambda1=lambda1Vec(iL1)
%
% AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version   : 2012/08/24
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nTS=length(data); % number of time series

    [N,p]=size(data{1});

    gaussParam = 0.2; % parameter of the Gaussian kernel

    eps = 10^-2; % stopping criterion threhold for residuals
    flagRes = 1; % use the above stopping criterion if 1 
    doPlot = 0; % shift to 1 to plot various figures (mean squared errors)

    % Other parameters
    level_alpha = 0.05; % level of statistical test
    NE_pick=0; 

    t0=0; tf=size(data{1},1)-1; % indices of initial and final times 

    numBS = 2*nBS;

    nL2 = length(lambda2Vec);
    nL1 = length(lambda1Vec);

    distJ = zeros(nL2,nL1);

    %% Generate block-bootstraps
    BSvec = cell(numBS,nTS);        
    timeInd = nan(numBS,Bsize+1);
    for j=1:numBS
        timeInd(j,:) = bb_ind(N,Bsize);
        for iTS=1:nTS
            sample = data{iTS};
            BSvec{j,iTS} = sample(timeInd(j,:),:);
        end
    end

    %% Run OKVAR-Boost on block-bootstraps
    for iL2=1:nL2
        for iL1=1:nL1
            for b=1:2:numBS
                J_b1 = zeros(p); J_b2 = zeros(p);
                for iTS=1:nTS
                    display(strcat(' lambda2=', num2str(lambda2Vec(iL2)),' lambda1=', num2str(lambda1Vec(iL1)), ' blockBS=', num2str(b), ' ts=', num2str(iTS)));

                    %% Block-bootstrap b1
                    name = strcat('sample',num2str(iTS),'_blockBS',num2str(b));

                    % Execute OKVAR-Boost
                    [listW,listC,rho,h,H_m,mse,idx_subset,lastIter]  = OKVARboost(BSvec{b,iTS},M,t0,tf,gaussParam,name,lambda2Vec(iL2),lambda1Vec(iL1),doPlot,randFrac,level_alpha,NE_pick,eps,flagRes);

                    % Compute the Jacobian for the OKVAR-Boost model
                    [J1] = computeJacobianOKVARboost(BSvec{b,iTS},gaussParam,listW,listC,rho,lastIter,idx_subset);

                    %% Block-bootstrap b2
                    name = strcat('sample',num2str(iTS),'_blockBS',num2str(b));

                    % Execute OKVAR-Boost
                    [listW,listC,rho,h,H_m,mse,idx_subset,lastIter]  = OKVARboost(BSvec{b+1,iTS},M,t0,tf,gaussParam,name,lambda2Vec(iL2),lambda1Vec(iL1),doPlot,randFrac,level_alpha,NE_pick,eps,flagRes);

                    % Compute the Jacobian for the OKVAR-Boost model
                    [J2] = computeJacobianOKVARboost(BSvec{b+1,iTS},gaussParam,listW,listC,rho,lastIter,idx_subset);

                    J_b1 = J_b1 + J1;
                    J_b2 = J_b2 + J2;
                end
                J_b1 = J_b1/nTS; J_b2 = J_b2/nTS;
                distJ(iL2,iL1) = distJ(iL2,iL1) + norm(J_b1-J_b2,'fro')^2;
            end

            distJ = distJ/nBS;
        end
    end
end