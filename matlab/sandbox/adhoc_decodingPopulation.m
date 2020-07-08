%% Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
%load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%% Across the population, which timescales best capture location tuning? 
% prereq : need to normalize all neurons to all the same sampling space
% thus making all spaces close : far
% timescales to consider
% 1) spike train of all neurons across availability period
% 2) mean of firing rate for all neurons across availability period 
% 3) mean of all touches
% 4) spike train of all touches from 5:45ms post touch 

