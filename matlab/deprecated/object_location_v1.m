%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

viewWindow = [-25:50]; %window for analyses around touch

%% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,viewWindow);

%% Plotter for feature tuning around touch window
gaussFilt = 1; %smoothing function for tuning plots
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
touchFeatureBinned_plotter(U,popV,selectedCells,fieldsList(1),whichTouches,viewWindow,gaussFilt)

%% Quantifying object location tuning
fieldsList = fields(popV{1}.allTouches);
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),whichTouches(1),viewWindow,'on');