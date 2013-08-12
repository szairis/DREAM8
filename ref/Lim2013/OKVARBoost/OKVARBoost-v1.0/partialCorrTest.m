function [Res,stop] = partialCorrTest(data,alpha,NE_pick)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : Test conditional independence between variables using partial
%            correlations
%
% USAGE :
% [Res,stop] = partialCorrTest(data,alpha,NE_pick,flag)
%
% INPUTS :
% data      : ([N,p] matrix) N is the number of p-dimensional observations
% alpha     : significance level of the test
% NE_pick   : number of edges to pick in each random subset
%             NE_pick=0 means that the all the significant edges are picked
%
% OUTPUTS :
% Res       : ([p,p] matrix) Res(i,j) = 1 if the partial correlation
%             between X_i and X_j given the other variables is significant
% stop      : true if the covariance matrix is singular
%
%
% AUTHOR    : Néhémy LIM, AMIS group, IBISC, Université d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% Version : 2012/05/16
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 3
        NE_pick = 0;
    end
    
    [N,p] = size(data);
    covMat = zeros(p); 
    Res = zeros(p);
    
    % compute the covariance matrix with time lag 1
    for i=1:p
        for j=1:i
            cov_ij = covariance(data(1:N-1,i),data(2:N,j)); 
            cov_ji = covariance(data(1:N-1,j),data(2:N,i));
            if abs(cov_ij(1)) >= abs(cov_ji(1))
                covMat(i,j) = cov_ij(1);
            else
                covMat(i,j) = cov_ji(1);
            end
            covMat(j,i) = covMat(i,j);
        end
    end
        
    stop = false;

    if (det(covMat) == 0 || isnan(det(covMat)))
        stop = true;
    else
        prMat = inv(covMat); % precision matrix
        partCorr = -repmat(1./sqrt(diag(prMat)),1,p).*prMat.*repmat(1./sqrt(diag(prMat)'),p,1); % partial correlations
        Stat = abs(atanh(partCorr)) * sqrt(N-(p-2)-3); % test statistic
        Test = Stat > norminv(1-0.5*alpha,0,1);
        
        if NE_pick == 0
            Res = Test;
        else
            SigStat = Stat.* Test; % significant statistic values


            ind = reshape(1:p*p,p,p);
            idxEdge = ind(tril(true(size(ind)),-1));

            sigStat = SigStat(find(tril(SigStat,-1)));

            [value,order] = sort(sigStat);
            idxSigEdge = idxEdge(order); % significant edges

            % If the network is sparse and can't choose m edges, then choose |sigEdge|
            last_m = min(NE_pick,length(idxSigEdge));

            if(NE_pick > last_m)
                warning(strcat('Chose_',num2str(last_m),'_edges instead of the requested_',num2str(NE_pick)));
            end

            idxPickedEdge = idxSigEdge(1:last_m);
            Res(idxPickedEdge) = 1;
            Res = symmetric(Res);
        end
    
    end
end