%% Load touch and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_clean.mat') %L5b excitatory cells recorded by Jon and Phil
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';
feature_list = {'pole'};
touch_struct = CompileTStruct(data_directory, feature_list);
%% Main Figures 4
% accuracy during recording sessions (B)
session_accuracy(U)

% Comparison of naive vs expert proportion of touch/OL cells (C)
unit_proportions(U,touch_struct)

% shape of tuning (D)
modulation_width_comparison(U,touch_struct)

% Heatmap of tuning (E)
heatmap_comparison(U,touch_struct)

%% Supplemental Figure 4 
Fig4_supplemental(U)
