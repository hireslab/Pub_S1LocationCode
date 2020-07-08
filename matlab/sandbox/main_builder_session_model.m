%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
% U = defTouchResponse(U,.95,'off');
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tStruct = object_location_quantification(U,touchCells,'pole','off');

%% Builder for identifying hilbert components that generate tuning
locationUnits = find(cellfun(@(x) x.is_tuned==1,tStruct));
trained = find(cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U));

% selectedUnits = intersect(touchCells,trained);
selectedUnits = locationUnits;
selectedArray = U(selectedUnits); 

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 5;

glmnetOpt.interpOption = 'yes';
glmnetOpt.touchDirection = []; %leave blank to use all touches 

glmnetOpt.downsampling_rate = 20; %in ms, must be divisible by 4000 
glmnetOpt.shift = -1:1; 

%GLMdesign Matrix Set-up
fileName = 'glm_session';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_session(selectedArray,glmnetOpt,glmModel);
end

%GLMdesign Matrix Build
selectedFeatures = [1:8]; 
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_session(glmModel,glmnetOpt,selectedFeatures);

parfor i = 12:40
    disp(['iterating for neuron ' num2str(i) '/' num2str(length(selectedArray))])

    if size(glmModel{i}.io.DmatXNormalized,1)<20
        disp(['skipping neuron ' num2str(i) 'b/c less than 20 trials'])
    else
        glmModel{i} = poissonModel_session(glmModel{i},glmnetOpt);
        glmModel{i}.meta = selectedUnits(i);
        glmModel{i}.name = fileName;
    end
end

cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
save(fileName,'glmModel','-v7.3')
