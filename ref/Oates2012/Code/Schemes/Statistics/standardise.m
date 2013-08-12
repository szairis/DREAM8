% STANDARDISE MEANS AND STANDARD DEVIATIONS
% Inputs:
% v = vector or matrix to standardise.
% dimension = direction in which to standardise.
% Output:
% z = standardised vector or matrix.
% © Chris Oates 2011.

function z = standardise(v,dimension)

if strcmp(dimension,'columns')
    v = v';
end

% standardise along rows
[N,M] = size(v);
z = zeros(N,M);
mu = mean(v,2);
sigma = std(v,0,2);
for n = 1:N
    z(n,:) = (v(n,:)-mu(n)) / sigma(n);
end
    
if strcmp(dimension,'columns')
    z = z';
end

end