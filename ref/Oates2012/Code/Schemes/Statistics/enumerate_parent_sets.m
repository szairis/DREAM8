% ENUMERATE PARENT SETS
% Input:
% P = number of variables.
% Output:
% models = P*M matrix containing all parent sets.
% © Chris Oates 2011.

function models = enumerate_parent_sets(P)

% consider parent sets of size at most max_in_degree
max_in_degree = min([P,2]);

% number of possible parent sets
n_models = 0;
for i=0:max_in_degree
    n_models = n_models + nchoosek(P,i);
end

% calculate all parent sets
models = false(P,n_models);
counter = 0;
for d = 0:max_in_degree
    nck = nchoosek(1:P,d);
    for i = 1:size(nck,1);
        counter = counter + 1;
        for j = 1:size(nck,2)
            models(nck(i,j),counter) = true;
        end
    end
end
end