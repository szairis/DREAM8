function [A,C] = B_l_L_0(y,t)

C = [1,0,1,0];

X = []; Y = []; delta = [];

for d = 1:size(y,3)
    X = [X y(:,1:end-1,d)];
    delta = [delta diff(t)];
    Y = [Y diff(y(:,:,d),1,2)];
end

for j = 1:length(delta)
    Y(:,j) = Y(:,j) / delta(j);
end

prods = 'linear';
lag = 1;

cd('Statistics')
A = Bayesian(Y,X,prods,lag); 
cd('../')