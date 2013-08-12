% GENERATE LONGITUDINAL SINGLE CELL DATASETS USING CANTONE
% Inputs: refer to curves.m function.
% Outputs:
% y = time course data P*(n+1)*C.
% A = true adjacency matrix.
% genes = names of the model variables.
% © Chris Oates 2011.

function [y,A,genes] = Cantone_concat(t,snr_cell,snr_meas,~,varargin)

A = [0 1 0 0 0; 0 0 1 0 0; 1 0 0 1 1; 0 -1 0 0 0; -1 0 0 0 0];
genes = cell(5,1);
genes(1) = {'Cbf1'};
genes(2) = {'Gal4'};
genes(3) = {'Swi5'};
genes(4) = {'Gal80'};
genes(5) = {'Ash1'};

% simulate the experiment
C = 10;
P = length(A(:,1));
y = zeros(P,length(t),C);
for c = 1:C
    y(:,:,c) = SDDE(t,snr_cell);
    
    % add measurement error
    noise = mean(mean(y(:,:,c))) / snr_meas;
    y(:,:,c) = mvnrnd(y(:,:,c)',noise^2*eye(P))';
end

end

function sol = SDDE(t,snr_cell)

t = round(t);

% average initial condition
x_0 = [0.0045; 0.0324; 0.0075; 0.0221; 0.012];

% multiplicative cell noise
diffusion = @(X) abs(diag(X)) / (100000*snr_cell);

% initial condition
P = length(x_0);
lag = 100;
noise = mean(x_0) / snr_cell;
initial = mvnrnd(x_0',noise^2*eye(P))';
x(:,1:lag) = repmat(initial,1,lag);

% numerically solve
T = ceil(t(end));
for n = 1:T
    x(:,lag+n) = x(:,lag+n-1) + drift(x(:,lag+n-1),x(:,n));
    x(:,lag+n) = mvnrnd(x(:,lag+n)',diffusion(x(:,lag+n-1)))';
end

sol = zeros(P,length(t));
for i = 1:length(t)
    sol(:,i) = x(:,lag+t(i));
end
end

% CANTONE DRIFT FUNCTION 
% © Cantone et al, Cell, 2009

function x_dot = drift(x,x_lag)

% prevent blow-up
x = max(x,zeros(size(x)));
x_lag = max(x_lag,zeros(size(x_lag)));

% parameters
alpha = [0 0.000149 0.003 0.00074 0.00061];
v = [0.04 0.000882 0.0201 0.0147 0.0182];
k = [1 0.0356 0.0372 0.01 1.814 1.814];
h = [1 1 1 4 1 1];
d = [0.0222 0.0478 0.4217 0.098 0.05];
gamma = 0.6;

% drift function
x_dot = zeros(5,1);
x_dot(1) = alpha(1) + v(1) * (x_lag(3)^h(1)) / ( (k(1)^h(1)+x_lag(3)^h(1)) * (1+(x(5)/k(2))^h(2)) ) - d(1)*x(1);
x_dot(2) = alpha(2) + v(2) * (x(1)^h(3)) / (k(3)^h(3)+x(1)^h(3)) - d(2)*x(2);
x_dot(3) = alpha(3) + v(3) * (x(2)^h(4)) / (k(4)^h(4)+x(2)^h(4) * (1+(x(4)/gamma)^4)) - d(3)*x(3);
x_dot(4) = alpha(4) + v(4) * (x(3)^h(5)) / (k(5)^h(5)+x(3)^h(5)) - d(4)*x(4);
x_dot(5) = alpha(5) + v(5) * (x(3)^h(6)) / (k(6)^h(6)+x(3)^h(6)) - d(5)*x(5);
end
