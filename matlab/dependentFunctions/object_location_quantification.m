
function [tuneStruct] = object_location_quantification(uberarray,selectedCells,hilbert_feature,displayOpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function is used to plot a heat map of location (e.g. angle at touch)
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


%function parameters
viewWindow = [-25:50]; %viewing window around touch for heatmap mainly
proportionDataPerBin = .05; %binning param
alpha_value = .01; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_touches_per_bin = 0; %minimum number of angle bins to consider quantifying
gauss_filt = .5; %2D smoothing filter for heatmap

%populating struct for tuning quantification
tuneStruct = cell(1,length(uberarray));
for i = 1:length(uberarray)
    tuneStruct{i}.calculations = [];
    tuneStruct{i}.is_tuned = nan;
end

rc = numSubplots(numel(selectedCells));
if willdisplay
    figure(22);clf
    figure(23);clf
end

for rec = 1:length(selectedCells)
    %stimulus and response variables definitions
    array = uberarray{selectedCells(rec)};
    [tVar] = atTouch_sorter(array,viewWindow);
    
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
    try
        rw = find(viewWindow == array.meta.touchProperties.responseWindow(1)) : find(viewWindow == array.meta.touchProperties.responseWindow(2));
    catch
        touch_cells = cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),uberarray);
        all_windows = cell2mat(cellfun(@(x) x.meta.touchProperties.responseWindow,uberarray(touch_cells),'uniformoutput',0)');
        touch_wind = median(all_windows); %filling empty touch windows w/ median touch windows
        rw = find(viewWindow == touch_wind(1)) : find(viewWindow == touch_wind(2));
    end
    
    response = mean(tVar.allTouches.R_ntk(:,rw),2) * 1000;
    numTouchesPerBin = round(numel(selected_feature).* proportionDataPerBin);   
    numBins = round(numel(selected_feature)./numTouchesPerBin);

    
    %% Heatmap of responses around touch. Window around touch defined as view_window
    if strcmp(hilbert_feature,'phase')
        [sorted_heat, sortedBy_heat] = binslin(selected_feature,tVar.allTouches.R_ntk,'equalE',13,-pi,pi);
    else
        [sorted_heat, sortedBy_heat] = binslin(selected_feature,tVar.allTouches.R_ntk,'equalN',numBins);
    end
    if willdisplay
        if numel(sortedBy_heat)>min_bins
            heat_resp = cell2mat(cellfun(@(x) mean(x,1),sorted_heat,'uniformoutput',0));
            smoothed_heat_resp = imgaussfilt(heat_resp,gauss_filt,'padding','replicate');
            figure(22);subplot(rc(1),rc(2),rec)
            imagesc(smoothed_heat_resp)
            hold on; plot([find(viewWindow==0) find(viewWindow==0)],[1 length(sortedBy_heat)],'-.w')
            caxis([0 prctile(smoothed_heat_resp(:),99)])
            if strcmp(hilbert_feature,'phase')
                set(gca,'ydir','normal','ytick',[1 12],'yticklabel',{'-\pi','\pi'},...
                    'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
            elseif strcmp(hilbert_feature,'pole')
                [~,zero_idx] = min(abs(cellfun(@median,sortedBy_heat)));
                set(gca,'ytick',zero_idx,'yticklabel',0,...
                    'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
            else
                set(gca,'ydir','normal','ytick',1:length(sortedBy_heat),'yticklabel',round(cellfun(@median,sortedBy_heat)),...
                    'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
            end
        else
            disp(['Not plotting cell ' num2str(selectedCells(rec)) ' b/c less than 5 sampled bins'])
        end
    end
    %% Tuning in touch response window
%     if strcmp(hilbert_feature,'phase') %commenting out for modulation
%         %index calculation. Need all variables to have same amount of bins
%         [sorted, sortedBy] = binslin(selected_feature,response,'equalE',13,-pi,pi);
%     else
        [sorted, sortedBy] = binslin(selected_feature,response,'equalN',numBins);
%     end
    disp(['num touches per bin is ' num2str(round(mean(cellfun(@numel,sorted)))) ' touches'])
    [quant_ol_p,tbl,stats] = anova1(cell2nanmat(sorted),[],'off');
    real_f = tbl{2,5};
    
%     %shuffle method     
    shuff_num = 1000;
    nanmat = cell2nanmat(sorted);
    f_stats = zeros(1,shuff_num); 
    for i = 1:shuff_num
        shuff = randperm(numel(nanmat));
        [~,tbl] = anova1(reshape(nanmat(shuff),size(cell2nanmat(sorted))),[],'off');
       f_stats(i) = tbl{2,5};
    end
    %percentile of real_f compared to shuffled f distribution
    comp_pctile = find(sort(f_stats) < real_f,1,'last') ./ 1000 * 100;
   
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
        else
            xq = min(x):1:max(x);
        end
        yq = interp1(x,y,xq);
        
        barsFit = barsP(yq,[min(xq) max(xq)],round(mean(cellfun(@numel,sorted))));
        barsFit.x = xq;
        
%         figure(9);subplot(rc(1),rc(2),rec)
%         bar(x,y,'k')
%         if sig_by_chance < 0.05 && quant_ol_p < 0.05 && numel(sortedBy)>min_bins
%             hold on; shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'g')
%         else
%             hold on; shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'k')
%         end
        
        smooth_response = barsFit.mean(2:end-1);
        smooth_stimulus = xq(2:end-1);
        
        [maxResponse,maxidx] = max(smooth_response);
        [minResponse,~] = min(smooth_response);
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
        tuneStruct{selectedCells(rec)}.calculations.mod_depth = (maxResponse - minResponse) ./ mean(barsFit.mean(2:end-1));
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = smooth_stimulus(maxidx);
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = (maxResponse - minResponse);
        tuneStruct{selectedCells(rec)}.calculations.anova_p = quant_ol_p;
    catch
        barsFit = [];
        disp(['skipping ' num2str(rec) ' due to ill fitting of bars'])
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = 0;
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = nan;
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = 0;
        tuneStruct{selectedCells(rec)}.calculations.mod_depth = 0;
        tuneStruct{selectedCells(rec)}.calculations.anova_p = quant_ol_p;
    end
    

    %for table 
    smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
    [~,maxidx] = max(smooth_response);
    [~,minidx] = min(smooth_response);
    tuneStruct{selectedCells(rec)}.calculations.responses_at_peak = sorted{maxidx};
    tuneStruct{selectedCells(rec)}.calculations.responses_at_trough = sorted{minidx};
    
    % making sure we've sampled enough bins before plotting.
    % min_bins defined as a global param above.
    if numel(sortedBy{1}) >= min_touches_per_bin
        
        if willdisplay
            figure(23);subplot(rc(1),rc(2),rec)
            try
                shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'k')
            catch
                shadedErrorBar(cellfun(@median, sortedBy), smooth(cellfun(@mean,sorted),smoothing_param),smooth(CI,smoothing_param),'r')
            end
        end
        
%         if sig_by_chance < 0.05 && quant_ol_p < 0.05
%          if quant_ol_p<alpha_value && mean(cellfun(@mean,sorted))>2
         if quant_ol_p < alpha_value && comp_pctile > 95
%          if sig_by_chance<=alpha_value 
            %plot smoothed response first
            smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
            %             if strcmp(array.meta.touchProperties.responseType,'excited')
            [maxResponse,maxidx] = max(smooth_response);
            [~, minidx] = min(smooth_response);
            
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
            
            tuneStruct{selectedCells(rec)}.is_tuned = 1;
            
        end
        
        if willdisplay
            if strcmp(hilbert_feature,'phase')
                set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
            elseif strcmp(hilbert_feature,'pole')
                set(gca,'xlim',[-1 1],'xtick',-1:1:1,'xdir','reverse')
            else
                set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
            end
        end
        
        tuneStruct{selectedCells(rec)}.stim_response.varNames = {'median S_ctk','mean R_ntk','std R_ntk','95CI R_ntk'};
        tuneStruct{selectedCells(rec)}.stim_response.values = [cellfun(@nanmedian, sortedBy) smooth(cellfun(@nanmean,sorted),smoothing_param) smooth(cellfun(@nanstd,sorted),smoothing_param) CI];
        tuneStruct{selectedCells(rec)}.stim_response.raw_stim = sortedBy;
        tuneStruct{selectedCells(rec)}.stim_response.raw_response = sorted;
    else
        tuneStruct{selectedCells(rec)}.is_tuned  = .5;  %not enough samples
    end
    
    if ~isempty(barsFit)
        tuneStruct{selectedCells(rec)}.stim_response.bars_fit= barsFit;
        tuneStruct{selectedCells(rec)}.stim_response.raw_stim = sortedBy;
        tuneStruct{selectedCells(rec)}.stim_response.raw_response = sorted;
    end
    
    
    
end