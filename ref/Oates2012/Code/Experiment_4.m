% INVESTIGATE THE BENEFIT OF LONGITUDINAL SINGLE CELL DATA
% © Chris Oates 2011.

function Experiment_4()

cd('Assessment')

num_datasets = 10;
snr_cell = 10;
snr_meas = 10;
num_cells = 1;

% sampling times
t_Cantone = [0:280/19:280];
t_Swat = [0:100/19:100];

% produce analysis for Cantone model
Investigate('Cantone_concat',num_datasets,t_Cantone,snr_cell,snr_meas,num_cells)

% produce analysis for Swat model
Investigate('Swat_concat',num_datasets,t_Swat,snr_cell,snr_meas,num_cells)

cd('../')
end