function sample_justification(uberarray,hilbert_feature)

disp('Iterating BARS fitting over cells and different bins. This may take a while...')

whisk_num_bins = [200, 100, 50, 20, 10];
touches_per_bin = [10 25 50 75 100 125 150 175 200];

whisk = [];
touch = [];

for rec = 1:numel(uberarray)
    disp([num2str(rec) '/' num2str(numel(uberarray))])
    % WHISKING
    array = uberarray{rec};
    spikes = squeeze(array.R_ntk);
    
    %whisking mask
    timePostTouchToTrim = 30;
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touchOnIdx = touchOnIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchOffIdx = touchOffIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(touchOnIdx)
        touchEx_mask(touchOnIdx(i):touchOffIdx(i)+timePostTouchToTrim) = NaN; %added 30 to add time from touch offset
    end
    touchEx_mask(1:100,:) = 1; %since bleedover from end of trials before, tihs ensure we keep end
    amplitude = squeeze(array.S_ctk(3,:,:));
    whisking = nan(size(squeeze(array.S_ctk(1,:,:))));
    whisking(amplitude>5)=1;
    whisking_mask = whisking .* touchEx_mask;
    
    if strcmp(hilbert_feature,'angle')
        conversion_feature = squeeze(array.S_ctk(1,:,:));
    elseif strcmp(hilbert_feature,'amplitude')
        conversion_feature = squeeze(array.S_ctk(3,:,:));
    elseif strcmp(hilbert_feature,'midpoint')
        conversion_feature = squeeze(array.S_ctk(4,:,:));
    elseif strcmp(hilbert_feature,'phase')
        conversion_feature = squeeze(array.S_ctk(5,:,:));
    end
    
    current_feature = conversion_feature(whisking_mask==1);
    filtered_spikes =spikes(whisking_mask==1);
    
    
    for b = 1:numel(whisk_num_bins)
        [sorted, sortedBy] = binslin(current_feature,filtered_spikes*1000,'equalN',whisk_num_bins(b));
        nanmat = cell2nanmat(sorted);
        [whisk.p_value{rec}(b),~,~] = anova1(nanmat,[],'off');
        
%         BARS fitting
        x = cellfun(@median,sortedBy);
        y = cellfun(@mean, sorted);
        if strcmp(hilbert_feature,'phase') || strcmp(hilbert_feature,'pole')
            xq = min(x):.1:max(x);
        else
            xq = min(x):1:max(x);
        end
        yq = interp1(x,y,xq);

        try 
            barsFit = [];
            barsFit = barsP(yq,[min(xq) max(xq)],round(mean(cellfun(@numel,sorted))));
            barsFit.x = xq;
            smooth_response = barsFit.mean(2:end-1);
            smooth_stimulus = xq(2:end-1);
            [maxResponse,maxidx] = max(smooth_response);
            [minResponse] = min(smooth_response);
            
            whisk.mod_idx{rec}(b) = (maxResponse - minResponse) ./ (maxResponse + minResponse);
            whisk.tune_peak{rec}(b) = smooth_stimulus(maxidx);
        catch
            disp(['skipping ' num2str(rec) ' due to ill fitting of bars'])
            whisk.mod_idx{rec}(b) = nan;
            whisk.tune_peak{rec}(b) = nan;
        end
        
        
    end
    
    %% TOUCH
    %stimulus and response variables definitions
    if isfield(array.meta.touchProperties,'responseWindow')
        touch_rw = array.meta.touchProperties.responseWindow;
        [tVar] = atTouch_sorter(array,touch_rw(1):touch_rw(2));
        
        if ~isempty(hilbert_feature)
            if strcmp(hilbert_feature,'angle')
                selected_feature = tVar.allTouches.S_ctk(:,1);
            elseif strcmp(hilbert_feature,'amplitude')
                selected_feature = tVar.allTouches.S_ctk(:,3);
            elseif strcmp(hilbert_feature,'midpoint')
                selected_feature = tVar.allTouches.S_ctk(:,4);
            elseif strcmp(hilbert_feature,'phase')
                selected_feature = tVar.allTouches.S_ctk(:,5);
            elseif strcmp(hilbert_feature,'curvature')
                selected_feature = tVar.allTouches.S_ctk(:,6);
            elseif strcmp(hilbert_feature,'pole')
                selected_feature = normalize_var(tVar.allTouches.S_ctk(:,7),-1,1);
            else
                error('select features of "angle", "amplitude", "midpoint", "phase", "curvature", or "pole"')
            end
        else
            error('select features of "angle", "amplitude", "midpoint", "phase", "curvature", or "pole"')
        end
        
        response = mean(tVar.allTouches.R_ntk,2) * 1000;
        
        for b = 1:numel(touches_per_bin)
            numBins = round(numel(selected_feature)./touches_per_bin(b));
            if numBins > 1
                [sorted_touch, sorted_by_touch] = binslin(selected_feature,response,'equalN',numBins);
                nanmat = cell2nanmat(sorted_touch);
                [touch.p_value{rec}(b),~,~] = anova1(nanmat,[],'off');
                
                %BARS fitting
                x = cellfun(@median,sorted_by_touch);
                y = cellfun(@mean, sorted_touch);
                if strcmp(hilbert_feature,'phase') || strcmp(hilbert_feature,'pole')
                    xq = min(x):.1:max(x);
                else
                    xq = min(x):1:max(x);
                end
                yq = interp1(x,y,xq);
                
                try
                    barsFit = [];
                    barsFit = barsP(yq,[min(xq) max(xq)],round(mean(cellfun(@numel,sorted))));
                    barsFit.x = xq;
                    smooth_response = barsFit.mean(2:end-1);
                    smooth_stimulus = xq(2:end-1);
                    [maxResponse,maxidx] = max(smooth_response);
                    [minResponse] = min(smooth_response);
                    
                    touch.mod_idx{rec}(b) = (maxResponse - minResponse) ./ (maxResponse + minResponse);
                    touch.tune_peak{rec}(b) = smooth_stimulus(maxidx);
                catch
                    disp(['skipping ' num2str(rec) ' due to ill fitting of bars'])
                    
                    touch.mod_idx{rec}(b) = nan;
                    touch.tune_peak{rec}(b) = nan;
                end
            else
                touch_ol_p{rec}(b) = nan;
            end
            nanmat = nan(1,numel(touches_per_bin));
            touch.mod_idx{rec} = nanmat;
            touch.tune_peak{rec} = nanmat;
            touch_ol_p{rec} = nanmat;
            
        end
    end
    
end

%% plotting p-value
selected_whisk = cellfun(@(x) x(3),whisk.p_value)<0.05;
figure(3840);clf

subplot(1,2,1);
p_mat = cell2mat(whisk.p_value((selected_whisk))')';
plot(whisk_num_bins,p_mat,'color',[.8 .8 .8])
hold on; plot(whisk_num_bins,nanmean(p_mat,2),'r')
hold on; plot(whisk_num_bins,nanmedian(p_mat,2),'r--')
set(gca,'ylim',[0.0001 .5],'xtick',fliplr(whisk_num_bins),'xticklabel',fliplr(1./whisk_num_bins)*100)
axis square
ylabel('ANOVA p-value');
xlabel('whisk data per bin (%)')

subplot(1,2,2);
t_mat = cell2mat(touch.p_value')';
plot(touches_per_bin,t_mat,'color',[.8 .8 .8])
hold on; plot(touches_per_bin,nanmean(t_mat,2),'r')
hold on; plot(touches_per_bin,nanmedian(t_mat,2),'r--')
set(gca,'ylim',[0.0001 .1],'xtick',touches_per_bin,'xdir','reverse')
axis square
ylabel('ANOVA p-value');
xlabel('number of touches per bin')

% modulation index plotting
mods_whisk = cell2mat(whisk.mod_idx(selected_whisk)')';
figure(5180);clf
subplot(2,2,1);
plot(whisk_num_bins,mods_whisk,'color',[.8 .8 .8])
hold on; plot(whisk_num_bins,nanmean(mods_whisk,2),'r')
hold on; plot(whisk_num_bins,nanmedian(mods_whisk,2),'r--')
set(gca,'xtick',fliplr(whisk_num_bins),'xticklabel',fliplr(1./whisk_num_bins)*100)
axis square
ylabel('modulation depth');
xlabel('whisk data per bin (%)')

subplot(2,2,2);
t_mat = cell2nanmat(touch.mod_idx');
plot(touches_per_bin,t_mat,'color',[.8 .8 .8])
hold on; plot(touches_per_bin,nanmean(t_mat,2),'r')
hold on; plot(touches_per_bin,nanmedian(t_mat,2),'r--')
set(gca,'xtick',touches_per_bin,'xdir','reverse')
axis square
ylabel('modulation depth');
xlabel('number of touches per bin')

% modulation tune peak plotting
mean_norm_tune = cellfun(@(x) x-nanmean(x), whisk.tune_peak((selected_whisk)),'uniformoutput',0);
whisk_peak = cell2mat(mean_norm_tune')';
subplot(2,2,3);
plot(whisk_num_bins,whisk_peak,'color',[.8 .8 .8])
hold on; plot(whisk_num_bins,nanmean(whisk_peak,2),'r')
hold on; plot(whisk_num_bins,nanmedian(whisk_peak,2),'r--')
set(gca,'xtick',fliplr(whisk_num_bins),'xticklabel',fliplr(1./whisk_num_bins)*100,'ylim',[-20 20])
axis square
ylabel('mean normalized tune peak');
xlabel('whisk data per bin (%)')

subplot(2,2,4);
mean_norm_tune = cellfun(@(x) x-nanmean(x), touch.tune_peak,'uniformoutput',0);
touch_peak = cell2nanmat(mean_norm_tune);
plot(touches_per_bin,touch_peak,'color',[.8 .8 .8])
hold on; plot(touches_per_bin,nanmean(touch_peak,2),'r')
hold on; plot(touches_per_bin,nanmedian(touch_peak,2),'r--')
set(gca,'xtick',touches_per_bin,'xdir','reverse','ylim',[-5 5])
axis square
ylabel('mean normalized tune peak');
xlabel('number of touches per bin')





