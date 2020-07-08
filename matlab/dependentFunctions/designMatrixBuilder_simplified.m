function [glmModel] = designMatrixBuilder_simplified(glmModel,glmnetOpt,selectedFeatures,interpOption)

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
    
    [nan_remove,~] = find(isnan(DmatX));
    
    if strcmp(glmnetOpt.touchDirection,'protraction')
        dir_remove = find((glmModel{i}.io.components.phase < 0 & glmModel{i}.io.components.pt_velocity > 0) ==0);
    elseif strcmp(glmnetOpt.touchDirection,'retraction')
        dir_remove = find((glmModel{i}.io.components.phase > 0 & glmModel{i}.io.components.pt_velocity < 0) ==0);
    else
        dir_remove = [];
        disp('using all touches since no input')
    end

    trialsToRemove = unique([nan_remove ; dir_remove]); 
    
    DmatX(trialsToRemove,:) = [] ;
    DmatY(trialsToRemove,:) = [] ;
    
    
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