function [A,C] = I_l_L_x(y,~)

C = [0,0,1,-1];

X = []; Y = []; 

for d = 1:size(y,3)
    X = [X y(:,1:end-1,d)];
    Y = [Y y(:,2:end,d)];
end

X = y(:,1:end-1);
Y = y(:,2:end);

prods = 'linear';
lag = 1;

cd('Statistics')
A = frequentist(Y,X,prods,lag); 
cd('../')