% Randomized indices for train and test set
% This assumes that there is no bias in the data and that all outcomes are
% sampled evenly 
%
% Input vector of prediction values and proportion to be used for training
% Output vector for each train and test indices

function [trainIdx,testIdx] = randomizedIndices(DmatY,proportionTrain)

if nargin<2
    error('must input proportion to set aside for train (e.g. .7)')
end

numTrials = numel(DmatY );
shuffIdx = randperm(numTrials); 
trainIdx = shuffIdx(1:numTrials*proportionTrain); 
testIdx = setdiff(shuffIdx,trainIdx);


