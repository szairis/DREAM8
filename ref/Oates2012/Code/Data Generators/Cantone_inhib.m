% GENERATE INHIBITION DATASETS USING CANTONE
% Inputs: refer to curves.m function.
% Outputs:
% y = time course data P*(n+1)*C.
% A = true adjacency matrix.
% genes = names of the model variables.
% © Chris Oates 2011.

function [y,A,genes] = Cantone_inhib(t,snr_cell,snr_meas,mod_S)

y = [];
for inhib_p = 1:5
    [yc,A,genes] = Cantone(t,snr_cell,snr_meas,mod_S,inhib_p);
    y(:,:,inhib_p) = yc;
end
end