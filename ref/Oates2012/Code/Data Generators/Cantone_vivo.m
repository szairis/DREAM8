% CANTONE IN VIVO SWITCH ON DATA
% Outputs: 
% y = time course P*(n+1).
% A = true adjacency matrix.
% genes = names of the model variables.
% © Chris Oates 2011.

function [y,A,genes] = Cantone_vivo(~,~,~,~,varargin)

A = [0 1 0 0 0; 0 0 1 0 0; 1 0 0 1 1; 0 -1 0 0 0; -1 0 0 0 0];
genes = cell(5,1);
genes(1) = {'Cbf1'};
genes(2) = {'Gal4'};
genes(3) = {'Swi5'};
genes(4) = {'Gal80'};
genes(5) = {'Ash1'};

load switch_on
y = data';

end
