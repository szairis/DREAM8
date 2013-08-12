% PERFORM BAYESIAN LINEAR REGRESSION
% Inputs:
% Y = P*n matrix of responses at t_1,...,T_T.
% X = P*n matrix of predictors at t_0,...,t_{T-1}.
% prods = number of products to look at in the design matrix. Takes
%         values 'linear' or 'interaction'.
% lag = binary variable indicating use of lagged predictors.
% Output:
% A = inferred adjacency matrix.
% © Chris Oates 2011.

function A = Bayesian(Y,X,prods,lag)

% enumerate parent sets
[P,~] = size(X);
models = enumerate_parent_sets(P);
    
% model likelihoods
Q = model_likelihoods(Y,X,models,prods,lag);
[P,M] = size(Q);
    
% sensibly weighted parent priors
edges = sum(models);
Inv_Pr = zeros(1,M);
for m = 1:M
    Inv_Pr(m) = nchoosek(P,edges(m));
end
Pr = bsxfun(@rdivide,ones(P,M),Inv_Pr);

% model posteriors
Po = Q .* Pr;
Po = bsxfun(@rdivide,Po,sum(Po,2));
    
% edge posterior probabilities
A = edge_post(Po,models);
end

% COMPUTE MODEL LIKELIHOODS
% Inputs:
% models = P*M list of possible parent sets.
% Output:
% Q = P*M matrix of model marginal likelihoods.
% © Chris Oates 2011.

function Q = model_likelihoods(Y,X,models,prods,lag)

[~,T] = size(X);
[P,M] = size(models);
Q = zeros(P,M);

Y = standardise(Y,'rows');

for m = 1:M
    
    % compute the design matrix
    B = design_matrix(X,m,models,prods,lag);
    [~,N] = size(B);
    
    % precompute some terms
    C = eye(T) - T/(T+1) * B*pinv(B);
    k = 1/(T+1)^(N/2);
    
    % calculate marginal likelihood
    for g = 1:P
        y_g = Y(g,:)';
        temp = y_g'*(C*y_g);
        Q(g,m) = k * temp^(-T/2);    
    end
end
end

% PERFORM BAYESIAN MODEL AVERAGING
% Inputs:
% Po = P*M matrix of model posterior probabilities.
% Output:
% A = inferred adjacency matrix.
% © Chris Oates 2011.

function A = edge_post(Po,models)

[P,M] = size(models);
A = zeros(P,P);

for g = 1:P
    for m = 1:M
        A(models(:,m)==1,g) = A(models(:,m)==1,g) + Po(g,m);
    end
end
end



    
    
    
    