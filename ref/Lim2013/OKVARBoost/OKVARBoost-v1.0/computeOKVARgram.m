function [K] = computeOKVARgram(data_train,data_test,B,gauss_param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Compute the OKVAR matrix-valued kernel from data
%
% USAGE :
% [K] = computeBlockKernel_v3(data_train,data_test,B,gauss_param)
%
% INPUTS :
% data_train    : ([n,p] matrix) training set
% data_test     : ([N,p] matrix) test set
% B             : ([p,p] matrix) symmetric positive semidefinite matrix 
% gauss_param   : (positive scalar) parameter of the Gaussian kernel  
%
% OUTPUT :
% K [N*p,n*p]   : Block Gram matrix
%
% AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version : 2011/07/14
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = size(data_train,1);
    T = size(data_test,1);
    K_cell = cell(T,n);

    for t=1:T
        for i=1:n
            K_cell{t,i} = B.*computeGaussKernel(gauss_param,data_train(i,:)',data_test(t,:)');
        end
    end

    K=cell2mat(K_cell);
end