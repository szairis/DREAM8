% GENERATE DATA FROM CANTONE MODEL
% Inputs: refer to curves.m function.
% varargin = index of a variable to inhibit.
% Outputs: 
% y = time course P*(n+1).
% A = true adjacency matrix.
% genes = names of the model variables.
% © Chris Oates 2011.

function [y,A,genes] = Cantone(t,snr_cell,snr_meas,mod_S,varargin)

A = [0 1 0 0 0; 0 0 1 0 0; 1 0 0 1 1; 0 -1 0 0 0; -1 0 0 0 0];
genes = cell(5,1);
genes(1) = {'Cbf1'};
genes(2) = {'Gal4'};
genes(3) = {'Swi5'};
genes(4) = {'Gal80'};
genes(5) = {'Ash1'};

% simulate the experiment
P = length(A(:,1));
y = zeros(P,length(t));
for s = 1:length(t)
    if nargin == 5
    y(:,s) = SDDE_ave(t(s),snr_cell,mod_S,varargin{1});
    else
        y(:,s) = SDDE_ave(t(s),snr_cell,mod_S);
    end
end    

% add measurement error
noise = mean(mean(y)) / snr_meas;
y = mvnrnd(y',noise^2*eye(P))';

end

% AVERAGE MANY CELLS DESTRUCTIVELY
% © Chris Oates 2011.

function y = SDDE_ave(t,snr_cell,mod_S,varargin)

% average initial condition
x_0 = [0.0045; 0.0324; 0.0075; 0.0221; 0.012];

% model a perfect inhibition
if nargin == 4
    x_0(varargin{1}) = 0.0001;
end

% empirical distribution threshold
mod_S_star = 30;

% multiplicative cell noise
diffusion = @(X) abs(diag(X)) / (100000*snr_cell);

% destructively sample cells at time t
P = length(x_0);
sample_size = min([mod_S mod_S_star]);
sample = zeros(P,sample_size);
lag = 100;
for cell = 1:sample_size
    noise = mean(x_0) / snr_cell;
    initial = mvnrnd(x_0',noise^2*eye(P))';
    if nargin == 4
        sample(:,cell) = SDDE(@drift,diffusion,lag,initial,t,varargin{1});
    else
        sample(:,cell) = SDDE(@drift,diffusion,lag,initial,t);
    end
end

% average mod_S cells
y = zeros(P,1);
for cell = 1:mod_S
    y = y + sample(:,unidrnd(sample_size));
end
y = y / mod_S;
end

% CUSTOM BUILT SDDE SOLVER
% Inputs:
% drift = drift function P*1.
% diffusion = diffusion matrix P*P.
% lag = lag in seconds.
% x_0 = initial state.
% t = SDDE is solved on [0 t].
% Outputs:
% sol = state vector at time t P*1.
% © Chris Oates 2011.

function sol = SDDE(drift,diffusion,lag,x_0,t,varargin)

res = 1;
t = res * round(t / res);
lag = res * round(lag / res);

% initial condition
x(:,1:lag) = repmat(x_0,1,lag);

% numerically solve
for n = 1:t
    if nargin == 6
        x(:,lag+n) = x(:,lag+n-1) + drift(x(:,lag+n-1),x(:,n),varargin{1})*res;
    else
        x(:,lag+n) = x(:,lag+n-1) + drift(x(:,lag+n-1),x(:,n))*res;
    end
    x(:,lag+n) = mvnrnd(x(:,lag+n)',diffusion(x(:,lag+n-1)))';
end
sol = x(:,lag+t);
end

% CANTONE DRIFT FUNCTION 
% © Cantone et al, Cell, 2009

function x_dot = drift(x,x_lag,varargin)

% prevent blow-up
x = max(x,zeros(size(x)));
x_lag = max(x_lag,zeros(size(x_lag)));

% model a perfect inhibition
if nargin == 3
    x(varargin{1}) = 0.0001;
end

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

% model a perfect inhibition
if nargin == 3
    x_dot(varargin{1}) = 0;
end
end
