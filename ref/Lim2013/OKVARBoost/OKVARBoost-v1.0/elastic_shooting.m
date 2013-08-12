function [c,cv,diff]  = elastic_shooting(K,y,lambda2,lambda1,c_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT :
% Solve the optimization problem 
%           argmin ||y-Kc||^2 + lambda2 * ||h||^2 + lambda1 * ||c||_1
% using subgradient and coordinate descent
%
% REFERENCE :
% W. Fu, Penalized regressions: The bridge versus the lasso, 1996
%
% USAGE :
% [c,cv,diff]  = elastic_shooting(K,y,lambda2,lambda1,c_init)
%
% INPUTS :
% K      : [n*p,n*p] Gram matrix of the kernels K(x_i,x_j) 
% y      : [n*p] vector of observations
% lambda2: positive parameter penalizing the complexity of the regressor
% lambda1: positive parameter penalizing non-sparse solutions
% c_init : [n*p] initialization for c
%
% OUTPUTS :
% c       : [n*p] estimation of the solution
% cv      : 1 if convergence, 0 if not
% diff    : (scalar) l1-norm of the difference between the current solution 
%           and the previous one 
%
% AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version : 2012/01/31
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    np = length(K);

    c = c_init;

    m=0;
    M=100;
    eps = 0.01*norm(y,1)/np;

    diff = eps+1;

    d = zeros(M,1);

    Z = K*(K+lambda2*eye(np));

    KTerm = 2*Z;
    yTerm = 2*K*y;

    while (m < M && diff > eps) % Stop after M iterations or when the solution is stable
        c_old = c;
        for j=1:np
            c(j) = 0;
            grad = KTerm*c - yTerm; %gradient of the Ridge-penalized Residual Sum of Squares
    %         Grad(:,j) = grad;
             if grad(j) > lambda1
                 c(j) = (lambda1-grad(j))/(2*Z(j,j));
             elseif grad(j) < -lambda1 
                 c(j) = (-lambda1-grad(j))/(2*Z(j,j));
             else
                 c(j) = 0;
             end
        end


        diff = norm(c - c_old,1);
        d(m+1) = diff;
        m=m+1;
    end

    cv = (diff <= eps);
end

    