function [A,C] = B_q_L_x(y,~)

C = [1,1,1,-1];

X = []; Y = [];

for d = 1:size(y,3)
    X = [X y(:,1:end-1,d)];
    Y = [Y y(:,2:end,d)];
end

prods = 'quad';
lag = 1;

cd('Statistics')
A = Bayesian(Y,X,prods,lag); 
cd('../')