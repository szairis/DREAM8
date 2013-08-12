function [A,C] = B_l_N_x(y,~)

C = [1,0,0,-1];

X = []; Y = [];

for d = 1:size(y,3)
    X = [X y(:,1:end-1,d)];
    Y = [Y y(:,2:end,d)];
end

prods = 'linear';
lag = 0;

cd('Statistics')
A = Bayesian(Y,X,prods,lag); 
cd('../')