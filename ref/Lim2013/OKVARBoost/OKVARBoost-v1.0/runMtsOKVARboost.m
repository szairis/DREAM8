function [Adj_thresh]=runMtsOKVARboost(data,nRun,M,randFrac,lambda2,lambda1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Multiple Time Series OKVAR-Boost run on T-cell data set
% T-cell data: 
% USAGE : [Adj_thresh]=runMtsOKVARboost(data,nRun,M,randFrac,lambda2,lambda1)
%
% INPUTS
% data         : cell array, data{ts} contains time-series number 'ts' 
%                data{ts} is an [n,p] matrix (n time points/p variables,genes) 
% nRun         : number of OKVAR-Boost runs
% M            : maximum number of iterations in a single OKVAR-Boost run
% randFrac     : size of random subsets (in percentage of total number of genes)
% lambda2(ridge), lambda1(L1) : regularization parameters
%
% OUTPUTS
% Adj_thresh   : cell array, Adj_thresh{freq} contains the estimated
%                adjacency matrix for the consensus threshold corresponding
%                to frequency 'freq'
%
% AUTHOR   : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% VERSION   : 2012/03/22 modified 2012/10/01
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nTS=length(data); % number of time series

    gaussParam = 0.2; % parameter of the Gaussian kernel

    freq=0.01:0.01:1; % edge frequency consensus threshold

    eps = 10^-2; % stopping criterion threhold for residuals
    flagRes = 1; % use the above stopping criterion if 1 
    doPlot = 0; % shift to 1 to plot various figures (mean squared errors)

    % Other parameters
    level_alpha = 0.05; % level of statistical test
    NE_pick=0; 

    listA=cell(1,nTS*nRun); % list of estimated adjacency matrices

    slope = exp(1); % slope in the tanh thresholding function
    t0=0; tf=size(data{1},1)-1; % indices of initial and final times 

    for iTS=1:nTS
        display(iTS);

    %% EXECUTE MULTIPLE RUN OKVAR-Boost

       for iRun=1:nRun
            display(strcat(' lambda1=', num2str(lambda1), ' lambda2=', num2str(lambda2),' run=',num2str(iRun),' ts=',num2str(iTS)));

            name = strcat('timeSeries',num2str(iTS),'_run',num2str(iRun));

            % Execute OKVAR-Boost
            [listW,listC,rho,h,H_m,mse,idx_subset,lastIter]  = OKVARboost(data{iTS},M,t0,tf,gaussParam,name,lambda2,lambda1,doPlot,randFrac,level_alpha,NE_pick,eps,flagRes);

            % Compute the Jacobian for the OKVAR-Boost model
            [J] = computeJacobianOKVARboost(data{iTS},gaussParam,listW,listC,rho,lastIter,idx_subset);


            J(1:length(J)+1:end) = 0; % remove self-edges

            Jn= J / norm(J,'fro'); % normalize the Jacobian

            % Threshold the Jacobian to get the estimated adjecency matrix
             listA{nRun*(iTS-1)+iRun} = round(tanh_param(slope,0,0,Jn));
       end
    end

    %% GET CONSENSUS NETWORKS
    Adj_thresh = graphConsensus(listA,freq); % list of unsigned adjacency matrices for each particular threshold
end