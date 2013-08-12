% INVESTIGATE THE BENEFIT OF INHIBITIONS
% © Chris Oates 2011.

function Experiment_5()

cd('Assessment')

num_datasets = 10;  
num_cells = 10000; 
snr_cell = 10;
snr_meas = 10;

% sampling times
t_Cantone = [0:280/19:280];
t_Swat = [0:100/19:100];

% produce analysis for Cantone model
Investigate('Cantone_inhib',num_datasets,t_Cantone,snr_cell,snr_meas,num_cells)

% produce analysis for Swat model
Investigate('Swat_inhib',num_datasets,t_Swat,snr_cell,snr_meas,num_cells)

cd('../')
end