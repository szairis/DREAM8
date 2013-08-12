% GENERATE DATA FROM SWAT MODEL
% Inputs: refer to curves.m function.
% varargin = index of a variable to inhibit.
% Outputs: 
% y = time course P*(n+1).
% A = true adjacency matrix.
% genes = names of the model variables.
% © Chris Oates 2011.

function [y,A,genes] = Swat(t,snr_cell,snr_meas,mod_S,varargin)

A = [1     1     1     0     1     1     0     1     0; ...
     1     1     1     0     1     0     0     1     0; ...
     0     0     1     1     0     0     0     0     0; ...
     1     0     1     1     0     1     0     0     0; ...
     0     0     1     0     1     0     0     0     0; ...
     1     1     1     0     1     1     1     1     0; ...
     0     0     0     0     0     1     1     0     0; ...
     0     0     0     0     0     0     0     1     1; ...
     0     0     0     0     0     1     1     1     1];
 
genes = cell(9,1);
genes(1) = {'pRB'};
genes(2) = {'E2F1'};
genes(3) = {'CycDi'};
genes(4) = {'CycDa'};
genes(5) = {'AP1'};
genes(6) = {'pRBp'};
genes(7) = {'pRBpp'};
genes(8) = {'CycEi'};
genes(9) = {'CycEa'};

% simulate the experiment
P = length(A(:,1));
y = zeros(P,length(t));
for s = 1:length(t)
    if nargin == 5
        y(:,s) = SDE_ave(t(s),snr_cell,mod_S,varargin{1});
    else
        y(:,s) = SDE_ave(t(s),snr_cell,mod_S);
    end
end    

% add measurement error
noise = mean(mean(y)) / snr_meas;
y = mvnrnd(y',noise^2*eye(P))';
end

% AVERAGE MANY CELLS DESTRUCTIVELY
% © Chris Oates 2011.

function y = SDE_ave(t,snr_cell,mod_S,varargin)

% average initial state
x_0 = [6; 1; 1; 1; 1; 1; 1; 1; 1];

% model a perfect inhibition
if nargin == 4
    x_0(varargin{1}) = 0.1;
end

% empirical distribution threshold
mod_S_star = 30;

% destructively sample cells at time t
P = length(x_0);
sample_size = min([mod_S mod_S_star]);
sample = zeros(P,sample_size);
for cell = 1:sample_size
    noise = mean(x_0) / snr_cell;
    initial = mvnrnd(x_0',noise^2*eye(P))';
    if nargin == 4
        sample(:,cell) = SDE_sol(initial,t,snr_cell,varargin{1});
    else
        sample(:,cell) = SDE_sol(initial,t,snr_cell);
    end
end

% average mod_S cells
y = zeros(P,1);
for cell = 1:mod_S
    y = y + sample(:,unidrnd(sample_size));
end
y = y / mod_S;
end

% SDE SOLVER
% Inputs:
% x_0 = initial state.
% t = SDE is solved on [0 t].
% Outputs:
% sol = state vector at time t P*1.
% © Chris Oates 2011.

function sol = SDE_sol(x_0,t,snr_cell,varargin)

t = round(t);
if t == 0
    sol = x_0;
else
    
    % multiplicative cell noise:
    diffusion = @(t,X) abs(diag(X)) / (400000*snr_cell);
    if nargin == 4
        D = @(t,X) drift(t,X,varargin{1});
    else
        D = @(t,X) drift(t,X);
    end
    s = sde(D,diffusion,'StartState',x_0);
    [paths,~,~] = s.simByEuler(t);
    sol = paths(end,:)';
end
end

% SWAT DRIFT FUNCTION
% © Swat et al, Bioinformatics, 2004

function z = drift(~,y,varargin)

% prevent blow-up
y = max(y,zeros(size(y)));

% model a perfect inhibition
if nargin == 3
    y(varargin{1}) = 0.1;
end

pRB = y(1);
E2F1 = y(2);
CycDi = y(3);
CycDa = y(4);
AP1 = y(5);
pRBp = y(6);
pRBpp = y(7);
CycEi = y(8);
CycEa = y(9);

% parameters
k1 = 1;
k2 = 1.6;
k3 = 0.05;
k16 = 0.4;
k34 = 0.04;
k43 = 0.01;
k61 = 0.3;
k67 = 0.7;
k76 = 0.1;
k23 = 0.3;
k25 = 0.9;
k28 = 0.06;
k89 = 0.07;
k98 = 0.01;
a = 0.04;
J11 = 0.5;
J12 = 5;
J15 = 0.001;
J18 = 0.6;
J61 = 5;
J62 = 8;
J65 = 6;
J68 = 7;
J13 = 0.002;
J63 = 2;
Km1 = 0.5;
Km2 = 4;
Km4 = 0.3;
Km9 = 0.005;
kp = 0.05;
opRB = 0.005;
oE2F1 = 0.1;
oCycDi = 0.023;
oCycDa = 0.03;
oAP1 = 0.01;
opRBp = 0.06;
opRBpp = 0.04;
oCycEi = 0.06;
oCycEa = 0.05;
Fm = 0.0065;

% drift function
z = zeros(9,1);
z(1) = k1 * E2F1/(Km1+E2F1) * J11/(J11+pRB) * J61/(J61+pRBp) - k16*pRB*CycDa + k61*pRBp - opRB*pRB;
z(2) = kp + k2 * (a^2+E2F1^2)/(Km2^2+E2F1^2) * J12/(J12+pRB) * J62/(J62+pRBp) - oE2F1*E2F1;
z(3) = k3*AP1 + k23*E2F1 * J13/(J13+pRB) * J63/(J63+pRBp) + k43*CycDa - k34*CycDi * CycDa/(Km4+CycDa) - oCycDi*CycDi;
z(4) = k34*CycDi * CycDa/(Km4+CycDa) - k43*CycDa - oCycDa*CycDa;
z(5) = Fm + k25*E2F1 * J15/(J15+pRB) * J65/(J65+pRBp) - oAP1*AP1;
z(6) = k16*pRB*CycDa - k61*pRBp - k67*pRBp*CycEa + k76*pRBpp - opRBp*pRBp;
z(7) = k67*pRBp*CycEa - k76*pRBpp - opRBpp*pRBpp;
z(8) = k28*E2F1 * J18/(J18+pRB) * J68/(J68+pRBp) + k98*CycEa - k89*CycEi * CycEa/(Km9+CycEa) - oCycEi*CycEi;
z(9) = k89*CycEi * CycEa/(Km9+CycEa) - k98*CycEa - oCycEa*CycEa;

% model a perfect inhibition
if nargin == 3
    z(varargin{1}) = 0;
end

% prevent blow-up
z = min(z,2*ones(size(z)));
end










