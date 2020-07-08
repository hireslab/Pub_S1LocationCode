function [glmModel] = designMatrixBlocks_trial_average(selectedArray,glmnetOpt,glmModel)


touchIndices.all = cellfun(@(x) ~isnan(squeeze(x.S_ctk(9,:,:))) + ~isnan(squeeze(x.S_ctk(12,:,:))),selectedArray,'uniformoutput',0)';
sampled_mask.all = cellfun(@(x) ones(size(x)),touchIndices.all,'uniformoutput',0);

[touchIndices.preDecision, sampled_mask.preDecision] = preDecisionTouchMat(selectedArray);

sampling_mask = cellfun(@(x) zeros(size(x)),touchIndices.all,'uniformoutput',0);
for d = 1:length(sampling_mask)
    sampling_period = round(mean(selectedArray{d}.meta.poleOnset)*1000):round(mean(selectedArray{d}.meta.poleOnset)*1000)+750; 
    sampling_mask{d}(sampling_period,:) = 1;
    sampled_mask.sampling_period{d} = sampling_mask{d};
    touchIndices.sampling_period{d} = touchIndices.all{d} .* sampling_mask{d};
end

varNames = {'theta','velocity','amplitude','midpoint','phase','kappa'};

for g = 1:length(selectedArray)
    array = selectedArray{g};
    varx=[1:5 17];
    
    %Input variables 
    for variableNumber=1:length(varx)
        
        [numTouches, meanFeature] = preDecisionTouchFeatures(array,touchIndices.sampling_period{g},varx(variableNumber),glmnetOpt.touchDirection,glmnetOpt.touchOrder);
        
        glmModel{g}.io.components.(varNames{variableNumber}) = nanmean(meanFeature)';
    end
    glmModel{g}.io.components.touchNumber = numTouches';
    
    %response variables
    glmModel{g}.io.DmatY = nansum(sampled_mask.sampling_period{g} .* squeeze(array.R_ntk))';
    
    %raw variables
    glmModel{g}.raw.pole = array.meta.motorPosition;
    
end

