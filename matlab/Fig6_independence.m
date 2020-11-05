%% Load whisking and neural time series struct
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';
load([data_directory 'excitatory_clean.mat']);
feature_list = {'angle','phase'};
whisk_struct = CompileWStruct(data_directory, feature_list);
touch_struct = CompileTStruct(data_directory, feature_list);

%% Plots
% Plots for Figure 6 and Supporting Figures 5 and 6 
% For SF 6, change hilbert_var to 'phase' 

hilbert_var = 'angle'; 
% proportion of neurons that intersect pie (A)
pie_comparison(touch_struct,whisk_struct,hilbert_var)

% whisk and touch tuning curves (B)
plot_tuning_curves(touch_struct,whisk_struct,hilbert_var) 

% whisk and touch modulation depth (C)
population_modulation_comparison(touch_struct,whisk_struct,hilbert_var)

% shape correlation of tuning curves (D)
intersect_correlation_v2(touch_struct,whisk_struct,hilbert_var); 

% scatter of tuning preference of whisk and touch (E)
tuning_preference_scatter(touch_struct,whisk_struct,hilbert_var)

% heatmaps (SF)
population_heatmap_builder(touch_struct,whisk_struct,hilbert_var)
