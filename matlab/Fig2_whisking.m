%% Load/build whisking structures 
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';
load([data_directory 'excitatory_clean.mat']);
feature_list = {'angle','phase','midpoint','amplitude','velocity'};
whisk_struct = CompileWStruct(data_directory, feature_list);
%% Main Figures 2
% Figure 2A: whisking vs quiet comparison
whisk_quiet_comparison(U)

% Figure 2B: angle and phase tuning curves
plot_phase_angle_curves(whisk_struct)

% Figure 2C: population heatmap of phase and angle tuning
plot_phase_angle_heat(whisk_struct)

% Figure 2E: Proportion of units tuned to phase, angle, neither, or both
plot_phase_angle_pie(whisk_struct)

% Figure 2F/G: Phase and angle modulation comparison 
plot_phase_angle_modulation(whisk_struct)

% Figure 2H-J: Hilbert decomposition comparison 
whisk_hilbert_comparison(U,whisk_struct)

%% Supplemental Figure 2
Fig2_supplemental(U,whisk_struct)
