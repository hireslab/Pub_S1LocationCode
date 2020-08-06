%% Load/build whisking structures 
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';
load([data_directory 'Raw\excitatory_all.mat']);
whisk_struct = CompileWStruct(data_directory);
%% Main Figures 2
% Figure 2A: whisking vs quiet comparison
WhiskQuietComparison(U)

% Figure 2B: angle and phase tuning curves
PlotPhaseAngleCurves(whisk_struct)

% Figure 2C: population heatmap of phase and angle tuning
PlotPhaseAngleHeat(whisk_struct)

% Figure 2E: Proportion of units tuned to phase, angle, neither, or both
PlotPhaseAnglePie(whisk_struct)

% Figure 2F/G: Phase and angle modulation comparison 
PlotPhaseAngleModulation(whisk_struct)

% Figure 2H-J: Hilbert decomposition comparison 
WhiskHilbertComparison(U,whisk_struct)

%% Supplemental Figure 2
 Fig2Supplemental(U,whisk_struct)
