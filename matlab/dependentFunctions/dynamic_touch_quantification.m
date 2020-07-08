function [tuneStruct] = dynamic_touch_quantification(uberarray,selectedCells,dynamic_touch_variable,displayOpt)

if (nargin < 4), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    rc = numSubplots(numel(selectedCells));
    figure(25);clf
end

dt_var_list = {'dkappa','dtheta'};
if ~any(contains(dt_var_list,dynamic_touch_variable))
    error('use dynamic touch variables of "dkappa" or "dtheta"');
end

numTouchesPerBin = 75; %number of touches to assign in each bin for quantification.
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying

tuneStruct = cell(1,numel(uberarray)); %dynamic touch tuning structure
dt = cell(1,numel(uberarray)); %dynamic touch variable structure
%% BUILD dynamic touch data structure

for rec = 1:length(selectedCells)
    array = uberarray{selectedCells(rec)};
    spikes = squeeze(array.R_ntk);
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touch_response_window = array.meta.touchProperties.responseWindow(1): array.meta.touchProperties.responseWindow(2);
    phase = squeeze(array.S_ctk(5,:,:)); 
    touch_phase = phase(touchOnIdx); 
    
    for i = 1:length(touchOnIdx)
        touch_window_idx = touchOnIdx(i):touchOffIdx(i);
        
        pt_velocity = mean(array.S_ctk(2,touchOnIdx(i)-5:touchOnIdx(i)-1));
        kwin=array.S_ctk(19,touch_window_idx); %get values in touch window
        if pt_velocity>0 && touch_phase(i)<0 %if protraction touch...
            [minValue ,val_idx] = min(kwin);
            dt{selectedCells(rec)}.stim.dkappa(i)=minValue;
        elseif pt_velocity<0 && touch_phase(i)>0 %if retraction touch...
            [maxValue ,val_idx] = max(kwin);
            dt{selectedCells(rec)}.stim.dkappa(i)=maxValue;
        else
            [~ ,val_idx] = max(abs(kwin)); %find idx of max kappa within each touch window, neg or pos
             dt{selectedCells(rec)}.stim.dkappa(i)=kwin(val_idx); %use idx to pull max kappa
        end
        dk_response_window = touch_window_idx(val_idx) + touch_response_window ;
        dt{selectedCells(rec)}.response.dkappa(i,:) = spikes(dk_response_window);
        
        dtheta=array.S_ctk(18,touchOnIdx(i):touchOffIdx(i)); %get values in touch window
        [~ ,maxidx] = max(abs(dtheta)); %find idx of max theta within each touch window, neg or pos
        dt{selectedCells(rec)}.stim.dtheta(i)=dtheta(maxidx); %use idx to pull max theta
        dt_response_window = touch_window_idx(maxidx) + touch_response_window ;
        dt{selectedCells(rec)}.response.dtheta(i,:) = spikes(dt_response_window);
    end
end

%% QUANTIFY dynamic touch data structure

for rec = 1:length(selectedCells)
    
    curr_dt = dt{selectedCells(rec)};
    numBins = round(size(curr_dt.stim.(dynamic_touch_variable),2) ./ numTouchesPerBin);
    
    [sorted, sortedBy] = binslin(curr_dt.stim.(dynamic_touch_variable),mean(curr_dt.response.(dynamic_touch_variable),2) * 1000,'equalN',numBins);
    
    if numBins >= 2
        x = cellfun(@nanmedian,sortedBy);
        y = smooth(cellfun(@nanmean,sorted),smoothing_param);
        err = smooth(cellfun(@(x)nanstd(x)./sqrt(numel(x)),sorted),smoothing_param);
        if willdisplay
            subplot(rc(1),rc(2),rec)
            shadedErrorBar(x,y,err)
            title(num2str((max(y) - min(y)) ./ (max(y) + min(y))))
        end
        
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (max(y) - min(y)) ./ (max(y) + min(y));
    else
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = 0;
    end
    
    tuneStruct{selectedCells(rec)}.dynamic_touch_raw = curr_dt;
end


