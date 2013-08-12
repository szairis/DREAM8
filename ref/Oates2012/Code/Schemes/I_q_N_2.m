function [A,C] = I_q_N_2(y,t)

C = [0,1,0,2];

X = []; Y = []; delta = [];

for d = 1:size(y,3)
    X = [X y(:,1:end-1,d)];
    delta = [delta diff(t)];
    Y = [Y diff(y(:,:,d),1,2)];
end

for j = 1:length(delta)
    Y(:,j) = Y(:,j) / delta(j);
end

h = delta.^(-2);
for j = 1:length(h)
    X(:,j) = X(:,j) / sqrt(h(j));
    Y(:,j) = Y(:,j) / sqrt(h(j));
end

prods = 'quad';
lag = 0;

cd('Statistics')
A = frequentist(Y,X,prods,lag); 
cd('../')




