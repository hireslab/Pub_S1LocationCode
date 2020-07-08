%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tStruct = object_location_quantification(U,selectedCells,'pole');

%%
tunedIdx = find( cellfun(@(x) x.is_tuned,tStruct)==1);
selectedArray = U(tunedIdx);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 20;

%BUILD parameters
glmnetOpt.buildIndices = [-25:50]; %Indices around touch
glmnetOpt.touchDirection = 'protraction';

fileName = 'glm_simplified';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_simplified(selectedArray,glmnetOpt,glmModel);
end
    
%GLMdesign Matrix Build
selectedFeatures = [1 5:8]; 
interpOption = 'off'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_simplified(glmModel,glmnetOpt,selectedFeatures,interpOption);

for i = 1:length(glmModel)
    disp(['iterating for neuron ' num2str(i) '/' num2str(length(tunedIdx))])
    glmModel{i} = poissonModel(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});
    glmModel{i}.meta = tunedIdx(i);
    glmModel{i}.name = fileName;
end
%% heatmap of input features and output predictions
builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));
meangof_de = cellfun(@(x) mean(x.gof.devExplained),glmModel(builtUnits));
wellfitUnits = intersect(builtUnits,find(meangof_de>.1)); 

figure(658);clf
for cellNum = wellfitUnits(5)
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
    imagesc([s_pole' features]')
    set(gca,'yticklabel',['pole' ; glmModel{cellNum}.io.selectedFeatures.name ; 'spikes'])
    colorbar
    colormap(new_map)
    title('inputs')
     
    %targets
    pred = cell2mat(glmModel{cellNum}.predicted.spikeProb');
    raw = cell2mat(glmModel{cellNum}.predicted.spikeTestRaw');
    pred_responses = pred(idx);
    raw_responses = raw(idx); 
   
 
    figure(658);subplot(4,1,4)
    imagesc(normalize_var([raw_responses  pred_responses],0,1)')
    set(gca,'ytick',1:2,'yticklabel',[{'true spikes'} ; {'predicted spikes'}])
    colorbar
    colormap(new_map)
    title(['outputs for gof de ' num2str(meangof_de(cellNum))])
    xlabel('touches sorted by pole position')
    
end

figure(348);clf
imagesc(corr([s_pole' features]))
axis square 
caxis([0 1])
colorbar
colormap bone
set(gca,'yticklabel',['pole' ; glmModel{cellNum}.io.selectedFeatures.name],...
    'xticklabel',['pole' ; glmModel{cellNum}.io.selectedFeatures.name])
axis square

%% goodness of fit of model 
samplesPerBin = 250; 
rc = numSubplots(length(glmModel)); 

figure(32);clf

% [~,idx] = sort(gof_tuning);
% plotIdx = fliplr(idx); 

for rec = 1:length(glmModel)
    
%     i = plotIdx(rec); 
    i = rec; 
    poles = normalize_var(cell2mat(glmModel{i}.predicted.pole'),-1,1);
    rawResponses = cell2mat(glmModel{i}.predicted.spikeTestRaw');
    predResponses = cell2mat(glmModel{i}.predicted.spikeProb');
    
    numBins = round(numel(poles)./samplesPerBin); 
    [raw_sorted,raw_sortedBy] = binslin(poles,rawResponses,'equalN',numBins);
    [pred_sorted,pred_sortedBy] = binslin(poles,predResponses,'equalN',numBins);
    
    figure(32);subplot(rc(1),rc(2),rec)
    hold on; plot(cellfun(@median, raw_sortedBy),normalize_var(smooth(cellfun(@mean, raw_sorted)),0,1),'b')
    hold on; plot(cellfun(@median, pred_sortedBy),normalize_var(smooth(cellfun(@mean, pred_sorted)),0,1),'r')
    
    gof_tuning(rec) = corr(smooth(cellfun(@mean, raw_sorted)),smooth(cellfun(@mean, pred_sorted)));
    mae(rec) = mean(smooth(cellfun(@mean, raw_sorted)) - smooth(cellfun(@mean, pred_sorted)));
    
    title(num2str(gof_tuning(i))); 
end

figure(80);clf
scatter(gof_tuning,gof_de,'filled','k')
title(['corr = ' num2str(corr(gof_tuning',gof_de'))])
ylabel('dev explained')
xlabel('tuning goodness of fit')

%% deviance explained leave one out feature importance
builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));
meangof_de = cellfun(@(x) mean(x.gof.devExplained),glmModel(builtUnits));
wellfitUnits = intersect(builtUnits,find(gof_tuning>.5)); 

glmModel = poisson_leaveOneOut(glmModel,glmnetOpt,wellfitUnits);
loo_importance = cell2mat(cellfun(@(x) nanmean(x.gof.LOO_importance,2),glmModel(wellfitUnits),'uniformoutput',0));

imp = cell2mat(cellfun(@(x) mean(x.gof.LOO_importance,2),glmModel(wellfitUnits),'uniformoutput',0));

imp(imp<0) = 0;
ranking = imp ./ nansum(imp);

figure(249);clf
subplot(1,2,1); 
hold on; plot(ranking(2:end,:),'color',[.9 .9 .9])
hold on;shadedErrorBar(1:size(imp,1)-1,nanmean(ranking(2:end,:),2),nanstd(ranking(2:end,:),[],2)./numel(wellfitUnits))

set(gca,'xtick',1:size(imp,1)-1,'xticklabel',glmModel{1}.io.selectedFeatures.name,'ytick',0:.5:1)
title('leave one out feature ranking')

figure(249);subplot(1,2,2)
bestFeatures = cell2mat(cellfun(@(x) mean(x.coeffs.raw(2:end,:),2),glmModel(wellfitUnits),'uniformoutput',0));
bfNorm = normalize_var(abs(bestFeatures),0,1);

hold on;plot(bfNorm,'color',[.9 .9 .9])
hold on;shadedErrorBar(1:size(bestFeatures,1),nanmean(bfNorm,2),nanstd(bfNorm,[],2)./sum(~isnan(sum(bfNorm))))
% hold on;shadedErrorBar(1:size(bestFeatures,1),nanmean(bfNorm,2),nanstd(bfNorm,[],2))

set(gca,'xtick',1:length(glmModel{1}.io.selectedFeatures.name),'xticklabel',glmModel{1}.io.selectedFeatures.name,'ytick',0:.5:1); 
ylabel('normalized feature weight')
title('absolute feature weight ranking')
suptitle(['n = ' num2str(sum(~isnan(sum(bfNorm)))) ' well fit units'])

%best feature in predicting choice 
[~,midx] = max(ranking(1:end,:));
keptIndices = wellfitUnits((midx-1)~=0); 
featureIndices = midx((midx-1)~=0)-1;

[sort_featureIdx,sort_idx] = sort(featureIndices);

topFeature = glmModel{1}.io.selectedFeatures.name(sort_featureIdx);

plotUnits = keptIndices(sort_idx);

figure(33);clf
rc = numSubplots(length(plotUnits)); 

samplesPerBin = 500;
for rec = 1:length(plotUnits)
    
    i = plotUnits(rec); 
    poles = normalize_var(cell2mat(glmModel{i}.predicted.pole'),-1,1);
    rawResponses = cell2mat(glmModel{i}.predicted.spikeTestRaw');
    predResponses = cell2mat(glmModel{i}.predicted.spikeProb');
    
    numBins = round(numel(poles)./samplesPerBin); 
    [raw_sorted,raw_sortedBy] = binslin(poles,rawResponses,'equalN',numBins);
    [pred_sorted,pred_sortedBy] = binslin(poles,predResponses,'equalN',numBins);
    
    subplot(rc(1),rc(2),rec)
    hold on; plot(cellfun(@median, raw_sortedBy),normalize_var(smooth(cellfun(@mean, raw_sorted)),0,1),'b')
    hold on; plot(cellfun(@median, pred_sortedBy),normalize_var(smooth(cellfun(@mean, pred_sorted)),0,1),'r')
    
    gof_tuning_indiv= corr(smooth(cellfun(@mean, raw_sorted)),smooth(cellfun(@mean, pred_sorted)));

    title([topFeature{rec} ' ' num2str(gof_tuning_indiv)]); 
end



%% random metrics 
%percent of touches responsive
pct_responsive = cellfun(@(x) mean(x.io.DmatY>0),glmModel);
figure(4);clf
histogram(pct_responsive,0:.1:1)

    
    
