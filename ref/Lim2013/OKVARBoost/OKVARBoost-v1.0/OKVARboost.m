function [W,C,rho,h,H_m,mse,subsets,stop]  = OKVARboost(data,M,t0,tf,gaussParam,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code comes with the paper entitled
% OKVAR-Boost : a novel boosting 
% algorithm to infer nonlinear dynamics and interactions in gene regulatory 
% networks
% from: Senbabaoglu, Y., Lim, N., Michailidis, M. and d'Alché-Buc, F.
% The two first authors have equally contributed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSTRACT : OKVAR-Boost procedure
%
% USAGE :
% [W,C,rho,h,H_m,mse,subsets,stop]  = OKVARboost(data,M,t0,tf,gaussParam,varargin)
%
% INPUTS :
% data      : ([N,p] matrix) N is the number of observations 
%             and p is the number of genes in the network
% M         : number of boosting iterations
% t0,tf     : initial and final times
% gaussParam   : positive parameter of the scalar Gaussian kernel  
%   
% varargin  : 9 optional inputs
%   name        : (string) name of saved files
%   lambda2     : (positive scalar) ridge penalty parameter
%   lambda1     : (positive scalar) l1 penalty parameter
%   doPlot      : plot and save regression results
%   randFrac    : size of random subset as a percentage of the network size
%   alpha       : level of the partial correlation test
%   NE_pick     : number of edges to pick in each random subset
%                 NE_pick=0 means that the all the significant edges are picked
%   eps         : stopping criterion threshold for residuals
%   flagRes     : if not 0, variables whose residuals are too low are
%                 removed at each iteration 
%
%
% OUTPUTS :
% W         : ([M,1] cell) list of [p,p] interaction matrices W_m
%             W_m[i,j] = 1 -> genes i and j interact
%             W_m[i,j] = 0 -> gene i and j do not interact
% C         : ([M,1] cell) list of C_m parameters
% rho       : (vector of length M), coefficients of the "line search" 
%             along the "steepest-descent"
% h         : ([M,1] cell) list of base models learned from each subset
% H_m       : ([N-1,p] matrix) estimated boosting model
% mse       : ([M,p] matrix) mean squared error for each gene at each iteration
% subsets   : ([M,1] cell) indices of selected genes at each iteration
% stop      : stopping iteration
%
% AUTHOR   : Néhémy LIM, AMIS group, IBISC, UniversitŽ d'Evry, France
% CO-AUTHOR : Yasin Senbabaoglu, Dept of statistics, University of Michigan, Ann Arbor, USA.
% VERSION   : 2012/02/15 
% This code is part of the ODESSA project - ANR-09-SYSC-009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Variables and default values
    
    numvarargs = length(varargin);
    if numvarargs > 9
        error('functions:OKVARboost:TooManyInputs','requires at most 9 optional inputs');
    end
    optargs = {'unknown' 0.001 1 0 1 0.05 0 0.0001 1};
    [optargs{1:numvarargs}] = varargin{:};
    [name,lambda2,lambda1,doPlot,randFrac,alpha,NE_pick,eps,flagRes] = optargs{:};

    %% Execution

    [N,p] = size(data);
    W = cell(M,1);
    C = cell(M,1);
    rho = zeros(M,1);
    h = cell(M,1);
    H_m = repmat(mean(data),N-1,1); % average vector of gene expressions, at the beginning, data are centered

    tspan = transpose(linspace(t0,tf,N));

    mse = zeros(M,p);

    subsets = cell(M,1);

    nFrac = round(randFrac*p); % size of random subsets 

    genes=1:p; nGenes=p;
    
    tot = 100; % maximum number of trials for the partial correlation test
    
    betaParam = .2; % diffusion parameter for the Laplacian
    
    stop = M; % stopping iteration

  
    for m=1:M
        if mod(m,10) == 0
            display(m);
        end
        
        U_m = data(2:N,:)-H_m; % (N-1)*p matrix of residuals

        % If the residual for gene i is too low then remove it
        if flagRes
            genesOut = [];
            for j=1:nGenes
                if norm(U_m(:,genes(j)))^2 < eps
                    genesOut = [genesOut j];
                end
            end

            genes(genesOut) = [];

            if(isempty(genes))
                stop = m-1;
                display(strcat('Stop at iteration_',num2str(stop),': No more significant residuals'));
                break
            end
        end
        
        
        
        %% Learn the interaction matrix with partial correlations

        terminate = true;
        if nFrac <= length(genes)
            nTry=0;
            while(nTry<tot && terminate)
                % Select a random subset
                idx_rand_m = sort(randsample(genes,nFrac,false)); % indices of genes selected at random
                
                % Partial Correlation test
                [W_m_sub,terminate] = partialCorrTest(U_m(:,idx_rand_m),alpha,NE_pick);

                if terminate % if no significant edge was found in the subnetwork, choose another one
                    nTry=nTry+1;
                else
                    W_m=zeros(p);
                    W_m(idx_rand_m,idx_rand_m) = W_m_sub;
                end
            end
        end

        if(terminate) %  if no significant edge was found after 'tot' trials
            stop = m-1;
            display(strcat('Stop at iteration_',num2str(stop),': Could not find significant edges'));
            break
        end

        nGenes = length(genes); % number of remaining genes
        
        %% Compute the block Gram matrix
       
        L = diag(sum(W_m_sub)) - W_m_sub; %Laplacian
        B_m = expm(betaParam*L);
        K_m = computeOKVARgram(U_m(:,idx_rand_m),U_m(:,idx_rand_m),B_m,gaussParam); % Gram matrix


        %% Learn the C coefficients
        Z = (K_m+lambda2*eye(size(K_m)));
        if (det(Z) ~= 0 && ~isnan(det(Z)))
            y = (reshape(U_m(1:N-1,idx_rand_m)',size(K_m,1),1));
            c_init = Z \ y; % Kernel Ridge solution as initial guess
            
            % Solve the Elastic Net with Fu's strategy
            c_m_k = elastic_shooting(K_m,y,lambda2,lambda1,c_init);
            
            if all(c_m_k == 0)
                stop=m-1;
                display(strcat('Stop at iteration_',num2str(stop),': All regression coefficients are zero'));
                break
            end
        else
            stop=m-1;
            display(strcat('Stop at iteration_',num2str(stop),': Matrix K_m+lambda*Id is singular'));
            break
        end

        h_m_k = transpose(reshape(K_m*c_m_k,nFrac,N-1));
        h_m = zeros(N-1,p); % genes that do not belong to the subset are assigned a 0 prediction
        h_m(:,idx_rand_m) = h_m_k;
        h{m} = h_m;

        W{m} = W_m;

        C_m_k = reshape(c_m_k,nFrac,N-1);

        C_m = zeros(p,N-1); % rows of C_m = 0 for genes that do not belong to the subset
        C_m(idx_rand_m,:) = C_m_k;

        C{m} = C_m;
        subsets{m} = idx_rand_m; 

        %% Line search
        rho(m) = trace(h_m'*U_m)/norm(h_m,'fro')^2; %rho(m) = \arg\min_\rho ||D_m - \rho * h_m||_2^2;

        %% Update the boosting model
        H_m = H_m + rho(m)*h_m;
   
        %% Compute the Mean Squared Errors
        mse(m,:) = (1/(N-1))*sum((data(2:N,:)-H_m).^2);
    
        %% Resize the outputs
        W = W(1:stop);
        C = C(1:stop); 
        rho = rho(1:stop);
        mse = mse(1:stop,:);
        subsets = subsets(1:stop);
        h = h(1:stop);
    end
    
    %% Plots
    if (doPlot == 1)
        %plot the data and the regressor
        figure(1)
        plot(tspan,data)
        title('Trajectory fit');
        xlabel('Index of time points');
        hold on
        plot(tspan(2:N),H_m,'+');
        hold off
        saveas(gcf,strcat('OKVAR-Boost_fit_',name,'.png'));
        close(1);

        %plot the mean squared error for each iteration
        figure(2)
        plot(1:stop,mse)
        title('Plot of genes''mean squared errors');
        xlabel('Iteration'); ylabel('mse')
        saveas(gcf,strcat('OKVAR-Boost_mse_',name,'.png'));
        close(2);

        % plot of the avearge mean squared error for each iteration
        figure(3)
        plot(1:stop,mean(mse,2))
        title('Plot of average mean squared error');
        xlabel('Iteration'); ylabel('mse')
        saveas(gcf,strcat('OKVAR-Boost_meanMse_',name,'.png'));
        close(3);
    end    
end
    