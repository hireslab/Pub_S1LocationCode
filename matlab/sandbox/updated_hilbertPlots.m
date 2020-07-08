%% pearson correlation of predicted FR and true FR
gauss_std = [1 2 4 8 16 32];
pearson_corr = zeros(length(gauss_std),length(glmModel));
for i = 1:length(gauss_std)
    
    buildIndices = glmModel{1}.modelParams.buildIndices; 
    shiftedResponseWindows = cellfun(@(x) find(U{x.meta}.meta.responseWindow(1) == buildIndices) : find(U{x.meta}.meta.responseWindow(2) == buildIndices) ,glmModel,'uniformoutput',0);
    
    selected_rawFR = cellfun(@(x,y) x.predicted.spikeTestRaw(:,y),glmModel,shiftedResponseWindows,'uniformoutput',0);
    rawFR = cellfun(@(x) imgaussfilt(x(:),gauss_std(i)),selected_rawFR,'uniformoutput',0);
    
    selected_predictedFR = cellfun(@(x,y) x.predicted.spikeProb(:,y),glmModel,shiftedResponseWindows,'uniformoutput',0);
    predictedFR = cellfun(@(x) x(:),selected_predictedFR,'uniformoutput',0);

    pearson_corr(i,:) = cellfun(@(x,y) corr(x,y),rawFR,predictedFR);
end
%% plotting of the above corr values
figure(123);clf
subplot(2,1,1)
xs = repmat((1:length(gauss_std))',1,length(glmModel)); 
scatter(xs(:),pearson_corr(:),'kx')
hold on; boxplot(pearson_corr')
set(gca,'xtick',1:length(gauss_std),'xticklabel',gauss_std)
xlabel('Gaussian sigma');
ylabel('Pearson correlation')

subplot(2,1,2)
[~,idx] = sort(pearson_corr(end,:));
imagesc(pearson_corr(:,idx))
set(gca,'ytick',1:length(gauss_std),'yticklabel',gauss_std,'xdir','normal','ydir','normal')
colorbar
xlabel('cell number');
ylabel('Gaussian sigma')

%% plotting of touch psth responses raw vs predicted
gauss_std = [1 2 4 8 16 32];
rc = numSubplots(length(glmModel));

for i = 3
    buildIndices = glmModel{1}.modelParams.buildIndices; 
    shiftedResponseWindows = cellfun(@(x) find(U{x.meta}.meta.responseWindow(1) == buildIndices) : find(U{x.meta}.meta.responseWindow(2) == buildIndices) ,glmModel,'uniformoutput',0);

    rawFR = cellfun(@(x) imgaussfilt(mean(x.predicted.spikeTestRaw)*1000,gauss_std(i)),glmModel,'uniformoutput',0);
    predictedFR = cellfun(@(x) mean(x.predicted.spikeProb)*1000,glmModel,'uniformoutput',0);
    
    [~,idx] = sort(pearson_corr(i,:));
    idx = fliplr(idx);
    
    figure(67+i);clf
    for b = 1:length(rawFR)
        cellNum = idx(b); %set to b to plot in build order
        subplot(rc(1),rc(2),b)
        plot(buildIndices,rawFR{cellNum},'k')
        hold on; plot(buildIndices,predictedFR{cellNum},'r')
        title(num2str(pearson_corr(i,cellNum)))
        set(gca,'xtick',0:25:50,'xlim',[min(buildIndices) max(buildIndices)])
    end
        
end
suptitle(['gaussian sigma ' num2str(gauss_std(i))])



%% deviance explained
DE = cellfun(@(x) mean(x.gof.devExplained),glmModel);
figure(5221);clf
for i = 3
hold on;scatter(DE,pearson_corr(i,:),'filled')
end

%% fitting to OL tuning
% tunedCells = cellfun(@(x) x.meta,glmModel);
% [~,tcxy_real] = object_location_quantification(U,tunedCells,'angle');

numTouchesPerBin = 75; %number of touches to assign in each bin for quantification.
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying
gauss_filt = .5;

rc = numSubplots(numel(tunedCells));


figure(22);clf
figure(24); clf

tcxy_predicted = cell(1,length(glmModel)); 

for rec = 1:length(glmModel)
    array = U{glmModel{rec}.meta};
    buildIndices = glmModel{rec}.modelParams.buildIndices;
    touchResponse = array.meta.responseWindow;
    trColumns = find(buildIndices==touchResponse(1)) : find(buildIndices==touchResponse(2));
    
    predResponseRaw = glmModel{rec}.predicted.spikeProb;
    predResponse = mean(glmModel{rec}.predicted.spikeProb(:,trColumns),2) * 1000;
    
    stimulus = glmModel{rec}.predicted.angles;
    numBins = round(numel(stimulus)./numTouchesPerBin);
    
    %% Heatmap
    
    [sorted_heat, sortedBy_heat] = binslin(stimulus,predResponseRaw,'equalN',numBins);
    
    
    if numel(sortedBy_heat)>min_bins
        figure(22);subplot(rc(1),rc(2),rec)
        heat_resp = cell2mat(cellfun(@(x) mean(x,1),sorted_heat,'uniformoutput',0));
        smoothed_heat_resp = imgaussfilt(heat_resp,gauss_filt,'padding','replicate');
        imagesc(smoothed_heat_resp)
        hold on; plot([find(buildIndices==0) find(buildIndices==0)],[1 length(sortedBy_heat)],'-.w')
%         caxis([0 prctile(smoothed_heat_resp(:),100)])
        
        set(gca,'ydir','normal','ytick',1:length(sortedBy_heat),'yticklabel',round(cellfun(@median,sortedBy_heat)),...
            'xtick',1:25:length(buildIndices),'xticklabel',0:25:50)
        
    else
        disp(['Not plotting cell ' num2str(selectedCells(rec)) ' b/c less than 5 sampled bins'])
    end
    
    %% model OL tuning
    
    [sorted, sortedBy] = binslin(stimulus,predResponse,'equalN',numBins);
    quant_ol_p(rec) = anova1(cell2nanmat(sorted),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    if numel(sortedBy)>min_bins
        figure(24);subplot(rc(1),rc(2),rec)
        shadedErrorBar(cellfun(@median, sortedBy), smooth(cellfun(@mean,sorted),smoothing_param),smooth(CI,smoothing_param),'k')
        
        if quant_ol_p(rec) < alpha_value
            %plot smoothed response first
            smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
            if strcmp(array.meta.responseType,'excited')
                [maxResponse,idx] = max(smooth_response);
                
                %plot scatter of first sig diff from max
                sd_p = nan(length(sorted),1);
                pThresh = .05;
                for g = 1:numel(sorted)
                    [~,sd_p(g)] = ttest2(sorted{idx},sorted{g});
                end
                all_idx = find(sd_p < pThresh); %all points sig diff from max
                [~,sd_idx_tmp] = min(abs(idx - all_idx));
                sd_idx = all_idx(sd_idx_tmp);
                
                if ~isempty(sd_idx) && ~isempty(idx)
                    hold on; scatter(median(sortedBy{idx}),maxResponse,'b','filled');
                    hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'b','filled');
                end
                
            elseif strcmp(array.meta.responseType,'inhibited')
                [minResponse,idx] = min(smooth_response);
                
                %plot scatter of first sig diff from min
                sd_p = nan(length(sorted),1);
                pThresh = .05;
                for g = 1:numel(sorted)
                    [~,sd_p(g)] = ttest2(sorted{idx},sorted{g});
                end
                all_idx = find(sd_p < pThresh); %all points sig diff from min
                [~,sd_idx_tmp] = min(abs(idx - all_idx));
                sd_idx = all_idx(sd_idx_tmp);
                
                if ~isempty(sd_idx) && ~isempty(idx)
                    hold on; scatter(median(sortedBy{idx}),minResponse,'r','filled');
                    hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'r','filled');
                end
            end

        end
        
        tcxy_predicted{rec} = [cellfun(@median, sortedBy) smooth(cellfun(@mean,sorted),smoothing_param)];

        set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
    end
end

%% compare the above 

%find intersecting bins
bins_predicted = cellfun(@(x) round(min(x(:,1))):round(max(x(:,1))),tcxy_predicted,'uniformoutput',0);
bins_real = cellfun(@(x) round(min(x(:,1))):round(max(x(:,1))),tcxy_real,'uniformoutput',0);

ix_bins = cellfun(@(x,y) intersect(x,y),bins_predicted,bins_real,'uniformoutput',0);

%stretching intepolations to match 
pred_response = cellfun(@(x,y) interp1(x(:,1),x(:,2),y)',tcxy_predicted,ix_bins,'uniformoutput',0);
real_response = cellfun(@(x,y) interp1(x(:,1),x(:,2),y)',tcxy_real,ix_bins,'uniformoutput',0);

tuning_correlation = cellfun(@(x,y) corr(x,y,'rows','complete'),pred_response,real_response);


best_mdled = intersect(find(pearson_corr(3,:)>.2),find(tuning_correlation>.5));


figure(19);clf
scatter(pearson_corr(3,:),tuning_correlation)
ylabel('tuning correlation')
xlabel('pearson correlation (sigma 4)')
hold on; scatter(pearson_corr(3,best_mdled),tuning_correlation(best_mdled),'filled','r')
title(['"well" modeled units = ' num2str(numel(best_mdled)) '/' num2str(numel(tuning_correlation))])

[~,idx] = sort(tuning_correlation);

idx = fliplr(idx);

figure(21);clf
for i = 1:length(pred_response)
    cellNum = idx(i); %set to i to 
    subplot(rc(1),rc(2),i)
    plot(1:length(pred_response{cellNum}),normalize_var(pred_response{cellNum},0,1),'r')
    hold on; plot(1:length(real_response{cellNum}),normalize_var(real_response{cellNum},0,1),'k')
    title(num2str(tuning_correlation(cellNum)))
end


%scatter of gof metrics x firing rate of cell
figure(580);clf
cellfiringRate = cellfun(@(x) nanmean(U{x.meta}.R_ntk(:)),glmModel);
scatter(cellfun(@(x) mean(x.predicted.spikeTestRaw(:)),glmModel) * 1000, tuning_correlation)
hold on; scatter(cellfun(@(x) mean(x.predicted.spikeTestRaw(:)),glmModel) * 1000, pearson_corr(3,:))
set(gca,'xtick',0:25:100)
legend({'tuning correlation','touch psth correlation'})
xlabel('touch response firing rate (Hz)')
ylabel('correlation')

corr((cellfun(@(x) mean(x.predicted.spikeTestRaw(:)),glmModel) * 1000)', tuning_correlation')
corr((cellfun(@(x) mean(x.predicted.spikeTestRaw(:)),glmModel) * 1000)', pearson_corr(3,:)')

%scatter of gof metrics x SNR
modeled_cells = cellfun(@(x) x.meta,glmModel);
[~,SNR] = defTouchResponse(U(modeled_cells),.95,'on');
figure(320);clf
scatter(SNR,tuning_correlation,'b')
hold on; scatter(SNR, pearson_corr(3,:),'r')
legend({'tuning correlation','touch psth correlation'})
xlabel('SNR of touch response')
ylabel('correlation')

corr(SNR', tuning_correlation')
corr(SNR', pearson_corr(3,:)')



%% kernels
best_mdled = intersect(find(pearson_corr(3,:)>.2),find(tuning_correlation>.5));
figure(128);clf

rc = numSubplots(numel(best_mdled));


for i = 1:length(best_mdled)
    selCell = best_mdled(i); 
    subplot(rc(1),rc(2),i)
    BI = glmModel{selCell}.modelParams.buildIndices;
%     coeffsToPlot = fields(glmModel{i}.coeffs);
    
    coeffsToPlot= glmModel{selCell}.io.selectedFeatures.name;
    for u = 1:length(coeffsToPlot)
        coeffsToPlotName = coeffsToPlot{u};
        if strcmp(coeffsToPlotName,'touch') || strcmp(coeffsToPlotName,'touchDur')
            hold on; plot(BI,sum(glmModel{selCell}.coeffs.(coeffsToPlotName)'.*glmModel{selCell}.basisFunctions.touch,2))
        else
            hold on; plot(BI,sum(glmModel{selCell}.coeffs.(coeffsToPlotName)'.*glmModel{selCell}.basisFunctions.features,2))
        end
    end
    set(gca,'xtick',[0:25:50])
end
%     ylabel('kernel weight')
%     xlabel('time from touch onset (ms)')
legend(coeffsToPlot(1:end))
suptitle('kernel weights in spike prediction')