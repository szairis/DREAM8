% INVESTIGATE CONSISTENCY OF THE INFERENCE SCHEMES
% © Chris Oates 2011.

function Experiment_3()

cd('Assessment')

num_datasets = 20;
snr_cell = 1000000000;
snr_meas = 1000000000;
num_cells = 1;

% sampling times
t_Cantone_even = [0:280/99:280];
t_Swat_even = [0:100/99:100];
t_Cantone_uneven = 0:log(280)/99:log(280);
t_Cantone_uneven = exp(t_Cantone_uneven);
t_Cantone_uneven(1) = 0;
t_Swat_uneven = 0:log(100)/99:log(100);
t_Swat_uneven = exp(t_Swat_uneven);
t_Swat_uneven(1) = 0;

% produce analysis for Cantone model, even sampling
Investigate('Cantone',num_datasets,t_Cantone_even,snr_cell,snr_meas,num_cells)

% produce analysis for Cantone model, even sampling
Investigate('Swat',num_datasets,t_Swat_even,snr_cell,snr_meas,1)

% produce analysis for Cantone model, uneven sampling
Investigate('Cantone',num_datasets,t_Cantone_uneven,snr_cell,snr_meas,1)

% produce analysis for Cantone model, uneven sampling
Investigate('Swat',num_datasets,t_Swat_uneven,snr_cell,snr_meas,1)

cd('../')
end