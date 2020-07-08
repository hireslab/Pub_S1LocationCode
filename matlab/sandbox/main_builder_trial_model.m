%Load whisking and neural time series struct 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
% U = defTouchResponse(U,.95,'off');
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tStruct = object_location_quantification(U,touchCells,'pole','off');

%% Poisson model for predicting spikes within pre-decision period 
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
glmnetOpt.numIterations = 25;

%BUILD parameters
glmnetOpt.touchDirection = 'pro';
glmnetOpt.touchOrder = 'all';

fileName = 'glm_trial_average';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_trial_average(selectedArray,glmnetOpt,glmModel);
end

selectedFeatures = [1 3:5 7];
interpOption = 'off'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_trial_average(glmModel,glmnetOpt,selectedFeatures);

for i = 1:length(glmModel)
    disp(['iterating for neuron ' num2str(i) '/' num2str(length(selectedArray))])

    if size(glmModel{i}.io.DmatXNormalized,1)<20
        disp(['skipping neuron ' num2str(i) 'b/c less than 20 trials'])
    else
%         glmModel{i} = poissonModel(glmModel{i}.io.DmatXNormalized, glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i}); 
        glmModel{i} = gaussianModel(glmModel{i}.io.DmatXNormalized, glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i}); 
        glmModel{i}.meta = selectedUnits(i);
        glmModel{i}.name = fileName;
    end
end



%% heatmap of input features and output predictions
builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));
meangof = cellfun(@(x) mean(x.gof.devExplained),glmModel(builtUnits));
wellfitUnits = intersect(builtUnits,find(meangof>.1)); 

figure(658);clf
for cellNum = wellfitUnits(3)
    %predictors
    [s_pole,idx] = sort(normalize_var(glmModel{cellNum}.raw.trimmedPole',1,-1));
    features = normalize_var(glmModel{cellNum}.io.DmatXNormalized(idx,:),-1,1);
    
    %colormap 
     stretch_resolution = 100 ;
    stretch_map = redbluecmap;
    new_map = nan(stretch_resolution,3);
    for j = 1:3
        new_map(:,j) = interp1(linspace(1,stretch_resolution,length(stretch_map)),stretch_map(:,j),1:stretch_resolution);
    end
    
    figure(658);subplot(4,1,1:3)
    imagesc([s_pole features]')
    set(gca,'yticklabel',['pole' ; glmModel{1}.io.selectedFeatures.name ; 'spikes'])
    colorbar
    colormap(new_map)
    title('inputs')
     
    %targets
    pred = cell2mat(glmModel{cellNum}.predicted.spikeProb');
    raw = cell2mat(glmModel{cellNum}.predicted.spikeTestRaw');
    [s_pole_mdl,s_idx] = sort(normalize_var(cell2mat(glmModel{cellNum}.predicted.pole),1,-1));
    pred = pred(s_idx);
    raw = raw(s_idx); 
    [~,sa,sb] = unique(s_pole_mdl);
    start_end_indices = [sa [sa(2:end)-1;length(s_pole_mdl)]];
    raw_responses = nan(length(start_end_indices),1);
    pred_responses = nan(length(start_end_indices),1);
    for g = 1:length(start_end_indices)
        raw_responses(g) = mean(raw(start_end_indices(g,1):start_end_indices(g,2)));
        pred_responses(g) = mean(pred(start_end_indices(g,1):start_end_indices(g,2)));
    end
    
 
    figure(658);subplot(4,1,4)
    imagesc(normalize_var([raw_responses  pred_responses],0,1)')
    set(gca,'ytick',1:2,'yticklabel',[{'true spikes'} ; {'predicted spikes'}])
    colorbar
    colormap(new_map)
    title(['outputs for gof de ' num2str(meangof(cellNum))])
    xlabel('trials sorted by pole position')
    
end

%% goodness of fit metrics 
builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));

meangof = cellfun(@(x) nanmean(x.gof.devExplained),glmModel(builtUnits));

figure(250);clf
[s_gof,idx] = sort(meangof);
plotUnits = builtUnits(fliplr(idx));
gof = fliplr(s_gof);

rc=numSubplots(length(plotUnits));
for g = 1:length(plotUnits)
    rec = plotUnits(g);

    subplot(rc(1),rc(2),g)

    if gof(g) > .1
        scatter(cell2mat(glmModel{rec}.predicted.spikeTestRaw'),cell2mat(glmModel{rec}.predicted.spikeProb'),'r.')
    else
        scatter(cell2mat(glmModel{rec}.predicted.spikeTestRaw'),cell2mat(glmModel{rec}.predicted.spikeProb'),'k.')
    end
    title(num2str(gof(g))); 
end
suplabel('true spike count')
suplabel('predicted spike count','y')


meangof(meangof<0) = 0; 
figure(230);
histogram(meangof,0:.05:1)
xlabel('gof: deviance explained');
ylabel('number of units'); 

%% deviance explained leave one out feature importance
builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));
meangof = cellfun(@(x) mean(x.gof.devExplained),glmModel(builtUnits));
wellfitUnits = intersect(builtUnits,find(meangof>.1)); 

glmModel = poisson_leaveOneOut(glmModel,glmnetOpt,wellfitUnits);
loo_importance = cell2mat(cellfun(@(x) nanmean(x.gof.LOO_importance,2),glmModel(wellfitUnits),'uniformoutput',0));

imp = cell2mat(cellfun(@(x) mean(x.gof.LOO_importance,2),glmModel(wellfitUnits),'uniformoutput',0));

imp(imp<0) = 0;
ranking = imp ./ nansum(imp);

figure(249);clf
subplot(1,2,1); 
shadedErrorBar(1:size(imp,1)-1,nanmean(ranking(2:end,:),2),nanstd(ranking(2:end,:),[],2)./numel(wellfitUnits))
set(gca,'xtick',1:size(imp,1)-1,'xticklabel',glmModel{1}.io.selectedFeatures.name)
title('leave one out feature ranking')

figure(249);subplot(1,2,2)
bestFeatures = cell2mat(cellfun(@(x) mean(x.coeffs.raw(2:end,:),2),glmModel(wellfitUnits),'uniformoutput',0));
bfNorm = normalize_var(abs(bestFeatures),0,1);

shadedErrorBar(1:size(bestFeatures,1),nanmean(bfNorm,2),nanstd(bfNorm,[],2)./sum(~isnan(sum(bfNorm))))
set(gca,'xtick',1:length(glmModel{1}.io.selectedFeatures.name),'xticklabel',glmModel{1}.io.selectedFeatures.name); 
ylabel('normalized feature weight')
title(['n = ' num2str(sum(~isnan(sum(bfNorm)))) ' well fit units'])


%%




