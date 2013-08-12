% INVESTIGATE PERFORMANCE UNDER UNEVEN SAMPLING
% © Chris Oates 2011.

function Experiment_2()

cd('Assessment')

num_datasets = 20;
num_cells = 10000;
snr_cell = 10;
snr_meas = 10;

% sampling times
t_Cantone = 0:log(280)/19:log(280);
t_Cantone = exp(t_Cantone);
t_Cantone(1) = 0;
t_Swat = 0:log(100)/19:log(100);
t_Swat = exp(t_Swat);
t_Swat(1) = 0;

% produce analysis for Cantone model
Investigate('Cantone',num_datasets,t_Cantone,snr_cell,snr_meas,num_cells)

% produce analysis for Swat model
Investigate('Swat',num_datasets,t_Swat,snr_cell,snr_meas,num_cells)

cd('../')
end