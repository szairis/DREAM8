function [ind] = bb_ind(origLen,Bstar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : bb_ind gives indices for block bootstrap
%   Time goes from 0 to origLen-1
%   Sampling is done modulo(origLen)
%
% USAGE : 
% [ind] = bb_ind(origLen,Bstar)
%
% INPUTS :
%   origLen : length of the original time-series
%   Bstar   : optimal block length
%   N       : number of bootstap samples you want
%
% OUTPUTS :
%   ind     : indices for block bootstrap
%
% AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% CO-AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% Version : 2012/05/16
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = randi(origLen,1,1);
    ind = mod(r:r+Bstar,origLen);
    ind(ind==0) = origLen;

end

