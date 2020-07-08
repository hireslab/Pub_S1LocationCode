%Load whisking and neural time series struct 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
is_tuned = object_location_quantification(U,selectedCells,'pole');

%% Builder for identifying hilbert components that generate tuning
tunedIdx = find(is_tuned==1);
selectedArray = U(tunedIdx);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 5;

%BUILD parameters
glmnetOpt.buildIndices = [0:50]; %Indices around touch

%basis function and convolution using gaussian distribution
glmnetOpt.bf.bfwidth =5;
glmnetOpt.bf.bfstd = 3;
glmnetOpt.bf.bfspacing = 3;
basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);
glmnetOpt.bf.indicesToAdd  = [-41:glmnetOpt.bf.bfspacing:-4];

%GLMdesign Matrix Set-up
fileName = 'glm_all_units';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_v2(selectedArray,glmnetOpt,glmModel);
end

%GLMdesign Matrix Build
selectedFeatures = [3:9]; 
interpOption = 'on'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_hilbert(glmModel,glmnetOpt,selectedFeatures,interpOption);

%Build model 
toBuild = find(~cellfun(@(x) isfield(x,'predicted'),glmModel)); 
parfor i = 36:45
 disp(['iterating for neuron ' num2str(i) '/' num2str(length(tunedIdx))])
 glmModel{i} = binomialModel_hilbert(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});
 glmModel{i}.meta = tunedIdx(i);
 glmModel{i}.name = fileName;
end

cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
save(fileName,'glmModel','-v7.3')
