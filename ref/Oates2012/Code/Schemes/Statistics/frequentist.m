% VARIABLE SELECTION USING AICc
% Inputs:
% Y is the P x T matrix of responses at t_1,...,T_(T+1).
% X is the P x T matrix of predictors at t_0,...,t_T.
% prods is the number of products to look at in the design matrix. Takes
% values 'linear' or 'interaction'.
% � Chris Oates 2011.

function A = frequentist(Y,X,prods,lag)

[P,T] = size(X);

% enumerate parent sets
models = enumerate_parent_sets(P);
    
% maximised log likelihoods
Q = max_log_likelihoods(Y,X,models,prods,lag);
    
% AICc model selection
A = AICc_select(Q,models,T);

end

function A = AICc_select(Q,models,T)

[P,M] = size(models);

AICc = zeros(P,M);
for m = 1:M
    % number of parameters
    q = sum(models(:,m));
    
    % corrected AIC score
    AICc(:,m) = 2*(q + q*(q+1)/(T-q-1))*ones(P,1) - 2*Q(:,m);
end

% identify the best models
A = false(P);
for p = 1:P
    [~,m] = min(AICc(p,:));
    A(:,p) = models(:,m);
end
end

function Q = max_log_likelihoods(Y,X,models,prods,lag)

[~,T] = size(X);
[P,M] = size(models);
Q = zeros(P,M);

Y = standardise(Y,'rows');

for m = 1:M
    % compute the design matrix
    B = design_matrix(X,m,models,prods,lag);
    
    % evaluate the maximised log likelihood
    for g = 1:P
        if rcond(B'*B)<0.001
            y_hat = B*((B'*B+0.001*eye(length(B'*B)))\(B'*Y(g,:)'));
        else
            y_hat = B*(pinv(B)*Y(g,:)');
        end
        sigma_hat_2 = 1/T * norm(Y(g,:)'-y_hat)^2;
        Q(g,m) = -T/2 * log(sigma_hat_2);
    end 
end
end




    
    
    
    