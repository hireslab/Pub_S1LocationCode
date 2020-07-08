function [glmModel] = designMatrixBuilder_hilbert(glmModel,glmnetOpt,selectedFeatures,interpOption)

%selectedFeatures is a vector with values corresponding to the values in
%the DmatFields 

DmatFields = fields(glmModel{1}.io.components);

for i = 1:length(glmModel)
    
    DmatX = []; 
    dims = []; 
    for g = 1:length(selectedFeatures)
        DmatX = [DmatX glmModel{i}.io.components.(DmatFields{selectedFeatures(g)})];
        dims = [dims size(glmModel{i}.io.components.(DmatFields{selectedFeatures(g)}),2)];
    end

    if strcmp(interpOption,'on')
        DmatX = fillmissing(DmatX,'linear',1);
    elseif strcmp(interpOption,'off')
        disp('no interpolation of missing values') 
    end
    
    DmatY = glmModel{i}.io.DmatY; 
    
    [row,~] = find(isnan(DmatX));
    trialsToRemove = unique(ceil(row./length(glmnetOpt.buildIndices)));
    trialsToRemoveIdx = (length(glmnetOpt.buildIndices)*(trialsToRemove-1)) + (1:(length(glmnetOpt.buildIndices))) ;
    DmatX(trialsToRemoveIdx,:) = [] ;
    DmatY(trialsToRemoveIdx,:) = [] ;
    
    
    glmModel{i}.io.DmatX = DmatX;
    glmModel{i}.io.DmatY = DmatY; 
    glmModel{i}.io.DmatXNormalized = (DmatX - mean(DmatX)) ./ std(DmatX);
    glmModel{i}.io.selectedFeatures.name = DmatFields([selectedFeatures]); 
    glmModel{i}.io.selectedFeatures.dims = dims;
    
    if isfield(glmModel{i},'raw')
        glmModel{i}.raw.trimmedPole = glmModel{i}.raw.pole;
        glmModel{i}.raw.trimmedPole(trialsToRemove) = [];
    end
    
end