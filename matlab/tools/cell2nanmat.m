function outputMat = cell2nanmat(inputCell) 

outputMat = nan(max(cellfun(@numel,inputCell)),numel(inputCell));

for i = 1:length(inputCell)
    outputMat(1:length(inputCell{i}),i) = inputCell{i};
end