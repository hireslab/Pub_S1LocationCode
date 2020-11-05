%% Load/build touch structures 
clear
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';
load([data_directory 'excitatory_all.mat']);
feature_list = {'pole','angle'};
touch_struct = CompileTStruct(data_directory, feature_list);
%% Main Figures 3
%Raster + PSTH one cell (A) 
trial_number = 29; %example trial in publication
plot_trial_raster(U, trial_number);

% firing rate X depth of recording (B)
plot_cell_depth(U, touch_struct)

% touch psth by quartiles of far, close and near (C)
units_to_plot = [13 ,5, 8];
plot_example_PSTH(U, touch_struct, units_to_plot);

% heatmap for object location tuned touch units (D)
feature = 'pole';
heatmap_builder(touch_struct,feature)

% modulation width (E)
modulation_width(touch_struct)
