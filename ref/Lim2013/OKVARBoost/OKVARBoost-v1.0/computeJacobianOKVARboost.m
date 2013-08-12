function [J,Jm] = computeJacobianOKVARboost(data,gaussParam,listW,listC,rho,mSTOP,subsets)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Compute the Jacobian for the Gaussian matrix-valued kernel
%
% USAGE :
% [J,Jm] = computeJacobianOKVARboost(data,gaussParam,listW,listC,rho,mSTOP,subsets)
%
% INPUTS :
% data          : ([N,p] matrix) N is the number of observations 
%                 and p is the number of genes in the network
% gaussParam    : parameter of the scalar Gaussian kernel 
% listW         : ([M,1] cell) list of interaction matrices
% listC         : ([M,1] cell) list of C_m parameters
% rho           : (vector of length M), coefficients of the "line search" 
%                 along the "steepest-descent" 
% mSTOP         : stopping iteration in OKVAR-Boost
% subsets       : indices of selected genes at each iteration
%
% OUTPUTS :
% J [p,p] : Jacobian matrix 
% Jm      : cell array, Jm{m} is the Jacobian matrix for model h_m
%
% AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% CO-AUTHOR   : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% Version : 2012/02/15
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nFrac = length(subsets{1});
    [N,p] = size(data);
    tf = N - 1;%number of time points when predictions are made 
    Jm = cell(1,mSTOP);
    J = zeros(p);
    data = data(1:tf,:);
    
    betaParam = .2; % diffusion parameter for the Laplacian
    display('Jacobian');
    
    for m=1:mSTOP
        W_m = listW{m};
        W_m_sub = W_m(subsets{m},subsets{m});
        
        deg_W = diag(sum(W_m_sub));
        L_m_sub = deg_W - W_m_sub; %standard Laplacian
        
        B_m_sub = expm(betaParam*L_m_sub);
        C_m = listC{m}; %reshape(listC{m},p,tf); % just to make sure it is correct size
        C_m_sub = C_m(subsets{m},:);
        tmpJ_sub = zeros(nFrac);
        tmpJ = zeros(p);
        % iterate over time points 1:tf
        for t=1:tf
            for l=1:tf
                K_lt = computeGaussKernel(gaussParam,data(l,subsets{m})',data(t,subsets{m})');
                tmpJ_sub = tmpJ_sub + (repmat(data(l,subsets{m})',1,nFrac)-repmat(data(t,subsets{m}),nFrac,1)).*K_lt.*repmat(C_m_sub(:,l)',nFrac,1);
            end
        end
        
        tmpJ(subsets{m},subsets{m}) = tmpJ_sub;
        B_m = zeros(p);
        B_m(subsets{m},subsets{m}) = B_m_sub;
        
        J = J + (2/tf)*gaussParam*rho(m)*B_m.*tmpJ;
        Jm{m} = J;
    end
end