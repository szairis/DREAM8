% GENERATE MULTIPLE LONGITUDINAL SINGLE CELL DATASETS USING SWAT
% Inputs: refer to curves.m function.
% Outputs:
% y = time course data P*(n+1)*C.
% A = true adjacency matrix.
% genes = names of the model variables.
% © Chris Oates 2011.

function [y,A,genes] = Swat_concat(t,snr_cell,snr_meas,~,varargin)

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

% average initial state
x_0 = [6; 1; 1; 1; 1; 1; 1; 1; 1];

% simulate the experiment
C = 10;
P = length(A(:,1));
y = zeros(P,length(t),C);
for c = 1:C
    noise = mean(x_0) / snr_cell;
    initial = mvnrnd(x_0',noise^2*eye(P))';
    y(:,:,c) = SDE_sol(initial,t,snr_cell);
    
    % add measurement error
    noise = mean(mean(y(:,:,c))) / snr_meas;
    y(:,:,c) = mvnrnd(y(:,:,c)',noise^2*eye(P))';
end

end

% SDE SOLVER
% Inputs:
% x_0 = initial state.
% t = SDE is solved on [0 t].
% Outputs:
% sol = state vector at time t P*1.
% © Chris Oates 2011.

function sol = SDE_sol(x_0,t,snr_cell)

% multiplicative cell noise:
diffusion = @(t,X) abs(diag(X)) / (400000*snr_cell);
D = @(t,X) drift(t,X);
s = sde(D,diffusion,'StartState',x_0);
t = round(t);
[paths,~,~] = s.simByEuler(t(end));

P = length(x_0);
sol = zeros(P,length(t));
for i = 1:length(t)
    sol(:,i) = paths(t(i)+1,:)';
end
end

% SWAT DRIFT FUNCTION
% © Swat et al, Bioinformatics, 2004

function z = drift(~,y,varargin)

% prevent blow-up
y = max(y,zeros(size(y)));

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
end

