function K = computeGaussKernel(gamma,X,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Compute the Gaussian matrix-valued kernel
%               K(i,j) = exp( -gamma * ( X(i) - Y(j) )^2 )
% 
% USAGE :
% K = computeGaussKernel(gamma,X,Y)
%
% INPUTS :
% gamma : positive parameter
% X     : vector of length p
% Y     : vector of length p
%
% OUTPUT :
% K [p*p] : Gaussian kernel matrix
%
% AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version   : 2011/06/21
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = length(X);
    K = exp(-gamma*((repmat(X,1,p)-repmat(Y',p,1)).^2));
end