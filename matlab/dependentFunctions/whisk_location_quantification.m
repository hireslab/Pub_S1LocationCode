function [tuneStruct] = whisk_location_quantification(uberarray,selectedCells,hilbert_feature,displayOpt,capture_window)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function is used to plot a heat map of whisking
% tuning for selectedCells. Specifically, use touch cells.
%
% inputs:
% uberarray - packaged uber array with all recorded units
% selectedCells - indices of units in uberarray that are touch cells
%
% outputs:
% whisk location tuning as defined from
% tuneStruct = struct with calculations of tuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 4), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

rc = numSubplots(numel(selectedCells));

%function parameters
alpha_value = .01; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 10; %smoothing parameter for smooth f(x) in shadedErrorBar
binslin_bins = 20; 

%populating struct for tuning quantification
tuneStruct = cell(1,length(uberarray));
for i = 1:length(uberarray)
    tuneStruct{i}.is_tuned = nan;
    tuneStruct{i}.calculations = [];
end
cell_counter = 0 
for rec = 1:length(selectedCells)
    array = uberarray{selectedCells(rec)};
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
    elseif strcmp(hilbert_feature,'velocity')
        conversion_feature = squeeze(array.S_ctk(2,:,:));
    elseif strcmp(hilbert_feature,'amplitude')
        conversion_feature = squeeze(array.S_ctk(3,:,:));
    elseif strcmp(hilbert_feature,'midpoint')
        conversion_feature = squeeze(array.S_ctk(4,:,:));
    elseif strcmp(hilbert_feature,'phase')
        conversion_feature = squeeze(array.S_ctk(5,:,:));
    elseif strcmp(hilbert_feature,'curvature')
        conversion_feature = squeeze(array.S_ctk(6,:,:));
    elseif strcmp(hilbert_feature,'pole') %conversion of angle to pole using touch positions
        viewWindow = -25:50;
        [tVar] = atTouch_sorter(array,viewWindow);
        pole_at_touch = normalize_var(tVar.allTouches.S_ctk(:,end),-1,1);
        angle_at_touch = tVar.allTouches.S_ctk(:,1);
        
        [s,sby] = binslin(pole_at_touch,angle_at_touch,'equalN',12); %bin responses based on pole positions
        cleaned = cell2mat(cellfun(@(x,y) rmoutliers([x y]),s,sby,'uniformoutput',0)); %remove outliers using median
        cleaned(sum(isnan(cleaned),2)>0,:) = [];
        clean_aat = cleaned(:,1);
        clean_pat = cleaned(:,2);
        p = polyfit(clean_aat,clean_pat,2); %use polyfit to convert angle to pole position
        
        angle = squeeze(array.S_ctk(1,:,:));
        conversion_feature = polyval(p,angle);
        
    end
    
    current_feature = conversion_feature(whisking_mask==1);
    whiskIdx = find(whisking_mask==1); 
    if strcmp(capture_window,'lag_window')
        try
            touch_wind = array.meta.touchProperties.responseWindow;
        catch
            touch_cells = cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),uberarray);
            all_windows = cell2mat(cellfun(@(x) x.meta.touchProperties.responseWindow,uberarray(touch_cells),'uniformoutput',0)');
            touch_wind = median(all_windows); %filling empty touch windows w/ median touch windows
        end
        response_idx = repmat(whiskIdx,1,length(touch_wind(1):touch_wind(2))) + repmat(touch_wind(1):touch_wind(2),numel(whiskIdx),1);
    elseif strcmp(capture_window,'lag')
        try
            touch_wind = array.meta.touchProperties.peakResponse;
        catch
            touch_cells = cellfun(@(x) isfield(x.meta.touchProperties,'peakResponse'),uberarray);
            all_windows = cell2mat(cellfun(@(x) x.meta.touchProperties.peakResponse,uberarray(touch_cells),'uniformoutput',0)');
            touch_wind = round(median(all_windows)); %filling empty touch windows w/ median touch windows
        end
        response_idx = whiskIdx+touch_wind;
    elseif strcmp(capture_window,'instant')
        response_idx = whiskIdx;
    end
    
    filt_spikes = spikes .* whisking_mask;
    filt_spikes = [filt_spikes nan(size(filt_spikes,1),1)]; %padding w/ nans at end for indexing
    
    stretch_touch_mask = [touchEx_mask nan(size(touchEx_mask,1),1)];
    whisks_with_touch = find(any(isnan(stretch_touch_mask(response_idx)),2));
    keep_whisks = 1:length(response_idx);
    filtered_spikes = nanmean(filt_spikes(response_idx),2);
    keep_whisks(unique([whisks_with_touch;find(isnan(filtered_spikes))])) = [];%only keep whisk examples that does not have any touch contamination or with at least 1 timepoint of spiking.
    
%     disp([capture_window ' : ' num2str(numel(keep_whisks) ./ length(response_idx).*100) '% of whisking timepoints w/ no touch contamination and 1 timepoint of spiking'])
    current_feature = current_feature(keep_whisks); 
    filtered_spikes = filtered_spikes(keep_whisks);
    
    %% Tuning in whisk windows
    %     equalN_numBins = round(sum(~isnan(current_feature(:)))./numWhiskSamplesPerBin);
    %     if strcmp(hilbert_feature,'phase')
    %         [sorted, sortedBy] = binslin(current_feature,filtered_spikes*1000,'equalE',equalE_bins,-pi,pi);
    %     else
    %         [sorted, sortedBy] = binslin(current_feature,filtered_spikes*1000,'equalE',equalE_bins,min(current_feature),max(current_feature));
    [sorted, sortedBy] = binslin(current_feature,filtered_spikes*1000,'equalN',binslin_bins);
    %     end
    
    nanmat = cell2nanmat(sorted);
    [quant_ol_p,tbl,stats] = anova1(nanmat,[],'off');
    real_f = tbl{2,5};
    
%     %shuffle method
    shuff_num = 100;
    nanmat = cell2nanmat(sorted);
    f_stats = zeros(1,shuff_num); 
    for i = 1:shuff_num
        shuff = randperm(numel(nanmat));
        [~,tbl] = anova1(reshape(nanmat(shuff),size(cell2nanmat(sorted))),[],'off');
       f_stats(i) = tbl{2,5};
    end
    %percentile of real_f compared to shuffled f distribution
    real_pctile = double(find(sort(f_stats) < real_f,1,'last') ./ shuff_num);
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    clear barsFit
    %BARSP FOR MOD IDX CALCULATIONS
    try
        x = cellfun(@median,sortedBy);
        y = cellfun(@mean, sorted);
        if strcmp(hilbert_feature,'phase') || strcmp(hilbert_feature,'pole')
            xq = min(x):.1:max(x);
        elseif strcmp(hilbert_feature,'velocity')
            xq = min(x): 200 : max(x);
        else
            xq = min(x):1:max(x);
        end
        yq = interp1(x,y,xq);
        
        numSamples = round(mean(cellfun(@numel,sorted))); 
        barsFit = barsP(yq,[min(xq) max(xq)],numSamples);
        barsFit.x = xq;
        
%         figure(9);subplot(rc(1),rc(2),rec)
%         bar(x,y,'k')
%         if sig_by_chance < 0.05 && quant_ol_p < 0.05 && nonzero_bins > 10
%             hold on; shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'g')
%         else
%             hold on; shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'k')
%         end
        
        smooth_response = barsFit.mean(2:end-1);
        smooth_stimulus = xq(2:end-1);
        
        [maxResponse,maxidx] = max(smooth_response);
        [minResponse,~] = min(smooth_response);
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = smooth_stimulus(maxidx);
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = (maxResponse - minResponse);
    catch
        barsFit = [];
        disp(['skipping ' num2str(rec) ' due to ill fitting of bars'])
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = 0;
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = nan;
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = 0;
    end
    
    %for table
    smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
    [~,maxidx] = max(smooth_response);
    [~,minidx] = min(smooth_response);
    tuneStruct{selectedCells(rec)}.calculations.responses_at_peak = sorted{maxidx};
    tuneStruct{selectedCells(rec)}.calculations.responses_at_trough = sorted{minidx};
    
    
    if willdisplay
        figure(23);subplot(rc(1),rc(2),rec)
        %             shadedErrorBar(smooth_stimulus,smooth_response,smooth(CI,smoothing_param),'k')
        %             shadedErrorBar( xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'k')
        
    end
    
    
%    if  quant_ol_p < alpha_value && ~isempty(barsFit) && mean(cellfun(@mean,sorted))>2 
   if  quant_ol_p < alpha_value && ~isempty(barsFit) && real_pctile > .95 && mean(cellfun(@mean,sorted))>2 
        %right and left tuning by idx
        compare_table = multcompare(stats,'Display','off');
        max_compares = compare_table(any(compare_table(:,[1 2]) == maxidx,2),:);
        sig_max_compares = max_compares(max_compares(:,end) < alpha_value,:);
        compare_idx = sig_max_compares(:,[1 2]);
        other_idx = compare_idx(compare_idx ~= maxidx);
        
        left_tuning_idx = other_idx(find(other_idx<maxidx,1,'last'));
        right_tuning_idx = other_idx(find(other_idx>maxidx,1,'first'));
        
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = median(sortedBy{maxidx});
        
        %finding tuning width before peak
        if ~isempty(left_tuning_idx)
            tuneStruct{selectedCells(rec)}.calculations.tune_left_width = median(sortedBy{maxidx})-median(sortedBy{left_tuning_idx});
            if willdisplay
                hold on; scatter(median(sortedBy{maxidx}),maxResponse,'g','filled');
                hold on; scatter(median(sortedBy{left_tuning_idx}),smooth_response(left_tuning_idx),'r','filled');
            end
        else
            tuneStruct{selectedCells(rec)}.calculations.tune_left_width = nan;
        end
        
        %finding tuning width after peak
        if ~isempty(right_tuning_idx)
            tuneStruct{selectedCells(rec)}.calculations.tune_right_width = median(sortedBy{right_tuning_idx}) - median(sortedBy{maxidx});
            if willdisplay
                hold on; scatter(median(sortedBy{maxidx}),maxResponse,'g','filled');
                hold on; scatter(median(sortedBy{right_tuning_idx}),smooth_response(right_tuning_idx),'r','filled');
            end
        else
            tuneStruct{selectedCells(rec)}.calculations.tune_right_width = nan;
        end
        
        if ~isnan(tuneStruct{selectedCells(rec)}.calculations.tune_right_width) | ~isnan(tuneStruct{selectedCells(rec)}.calculations.tune_left_width)
            tuneStruct{selectedCells(rec)}.is_tuned = 1;
        else
            tuneStruct{selectedCells(rec)}.is_tuned = 0;
        end
        tuneStruct{selectedCells(rec)}.calculations.responses_at_peak = sorted{maxidx};
        tuneStruct{selectedCells(rec)}.calculations.responses_at_trough = sorted{minidx};
    else
        tuneStruct{selectedCells(rec)}.is_tuned = 0; 
    end
    
    if willdisplay
        if strcmp(hilbert_feature,'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
        elseif strcmp(hilbert_feature,'pole')
            set(gca,'xtick',-5:1:5,'xdir','reverse')
        else
            set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
        end
    end
    
    if tuneStruct{selectedCells(rec)}.is_tuned == 1
        cell_counter = cell_counter+1;
%         [quant_ol_p, real_pctile]
%         figure(280);clf;scatter(cellfun(@median,sortedBy),cellfun(@mean,sorted))
        disp(['num tuned cells = ' num2str(cell_counter) '/' num2str(rec)])
%         pause
    end
    
    tuneStruct{selectedCells(rec)}.stim_response.varNames = {'median S_ctk','mean R_ntk','std R_ntk','95CI R_ntk'};
    tuneStruct{selectedCells(rec)}.stim_response.values = [cellfun(@nanmedian, sortedBy) smooth(cellfun(@nanmean,sorted),smoothing_param) smooth(cellfun(@nanstd,sorted),smoothing_param) smooth(CI,smoothing_param)];
    tuneStruct{selectedCells(rec)}.stim_response.raw_stim = sortedBy;
    tuneStruct{selectedCells(rec)}.stim_response.raw_response = sorted;
    if ~isempty(barsFit)
        tuneStruct{selectedCells(rec)}.stim_response.bars_fit= barsFit;
        tuneStruct{selectedCells(rec)}.stim_response.bars_stim = xq;
    end
    
end
