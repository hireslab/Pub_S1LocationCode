load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\glm_session.mat')

builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));

for i = 1:length(builtUnits)

    trial_length = glmModel{builtUnits(i)}.modelParams.trial_length;
    full_x = trial_length + sum(glmModel{i}.modelParams.shift~=0);
    x_values = linspace(0,4000,full_x);
        
        
        true = cellfun(@(x) reshape(x,trial_length,numel(x)./trial_length), glmModel{builtUnits(i)}.predicted.spikeTestRaw,'uniformoutput',0);
        predicted = cellfun(@(x) reshape(x,trial_length,numel(x)./trial_length), glmModel{builtUnits(i)}.predicted.spikeProb,'uniformoutput',0);
        polePositions = cell2mat(glmModel{builtUnits(i)}.predicted.pole);
    
        true_psth = mean(cell2mat(true)',2);
        predicted_psth =  mean(cell2mat(predicted)',2);
    
        [sort_values,idx] = sort(polePositions);
        [~,ia] = unique(sort_values)
        
        full_true = cell2mat(true)'; 
        full_predict = cell2mat(predicted)';
        sorted_raster_true = full_true(idx(ia),:);
%         sorted_raster_predict = poissrnd(full_predict(idx(ia),:));
         sorted_raster_predict = full_predict(idx,:);
    
        ceil_value = prctile(sorted_raster_true(:),99);

        figure(40);clf
        subplot(2,1,1)
        imagesc(sorted_raster_true)
        set(gca,'ytick',[])
        caxis([0 ceil_value])
        ylabel('close - far')
        subplot(2,1,2);
        imagesc(sorted_raster_predict)
        caxis([0 ceil_value])
        newcmap = flipud(bone)
        colormap(newcmap)
        set(gca,'ytick',[])
        ylabel('close - far')

        
        %psth
        figure(51);clf
        plot(smooth(nanmean(sorted_raster_true)),'k')
        hold on; 
        plot(smooth(nanmean(sorted_raster_predict)),'r')
        

    meanGOF(i) = mean(glmModel{builtUnits(i)}.gof.devExplained);
end

