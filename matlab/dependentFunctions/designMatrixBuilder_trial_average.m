function [glmModel] = designMatrixBuilder_trial_average(glmModel,glmnetOpt,selectedFeatures)

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
    
    DmatY = glmModel{i}.io.DmatY; 
    
    [nan_remove,~] = find(sum(isnan(DmatX),2)>0);
    
    trialsToRemove = unique([nan_remove ]); 
    
    DmatX(trialsToRemove,:) = [] ;
    DmatY(trialsToRemove,:) = [] ;
    
    
    glmModel{i}.io.DmatX = DmatX;
    glmModel{i}.io.DmatY = DmatY; 
    glmModel{i}.io.DmatXNormalized = (DmatX - mean(DmatX)) ./ std(DmatX);
    
    if any(any(isnan(glmModel{i}.io.DmatXNormalized)))
        glmModel{i}.io.DmatXNormalized = [];
    end

    glmModel{i}.io.selectedFeatures.name = DmatFields([selectedFeatures]); 
    glmModel{i}.io.selectedFeatures.dims = dims;
    
    if isfield(glmModel{i},'raw')
        glmModel{i}.raw.trimmedPole = glmModel{i}.raw.pole;
        glmModel{i}.raw.trimmedPole(trialsToRemove) = [];
    end
    
end