% COMPUTE THE DESIGN MATRIX
% Inputs:
% X = P*n matrix of predictor variables.
% m = index of required model.
% Outputs:
% B = design matrix.
% © Chris Oates 2011.

function B = design_matrix(X,m,models,prods,lag)

% include lagged predictors
delay = 10; % delay about T/delay minutes
temp0 = X(models(:,m),:);
if lag == 1
    [~,T] = size(X);
    tau = ceil(T/delay); 
    temp1 = temp0(:,1:end-tau); 
    temp1 = [repmat(temp1(:,1),1,tau) temp1];
    temp = [temp0; temp1];
else
    temp = temp0;
end

% compute the design matrix
B = temp';
B = x2fx2(B,prods);
B(:,1) = [];
B = standardise(B,'columns');
end

function z = x2fx2(B,prod)
if strcmp(prod,'cubic')
    [Q,P] = size(B);
    B = [ones(Q,1) B];
    z = [];
    for i = 1:P+1
        for j = i:P+1
            for k = j:P+1
                z = [z B(:,i).*B(:,j).*B(:,k)];
            end
        end
    end
else 
    z = x2fx(B,prod);
end
end