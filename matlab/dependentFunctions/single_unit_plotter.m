%% plot for each cell
hilbert_feature = 'angle'

fileName = ['whisk_' hilbert_feature '_tune'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    wStruct = whisking_location_quantification(U,1:numel(U),hilbert_feature,'off');
end
%%
figure(8);clf
% for rec = datasample(1:length(U),1)
for rec = 32
    clearvars -except rec U wStruct hilbert_feature
    
    [tVar] = atTouch_sorter(U{rec},-25:50);
    touch_response = nanmean(tVar.allTouches.R_ntk) * 1000;
    
    subplot(1,3,1)
    bar(-25:50,touch_response,'k')
    title('touch response')
    xlabel('time from touch onset (ms)')
    ylabel('firing rate (spks/s)')
    
    whisk_xy = [wStruct{rec}.stim_response.bars_fit.x' wStruct{rec}.stim_response.bars_fit.mean...
        wStruct{rec}.stim_response.bars_fit.confBands - wStruct{rec}.stim_response.bars_fit.mean];
    norm_whisk =[wStruct{rec}.stim_response.bars_fit.x' normalize_var(whisk_xy(:,2),0,1)];
    
    if wStruct{rec}.is_tuned==1
        subplot(1,3,2)
        hold on; shadedErrorBar(whisk_xy(:,1),whisk_xy(:,2),whisk_xy(:,3),'c')
        subplot(1,3,3)
        hold on; plot(norm_whisk(:,1),norm_whisk(:,2),'c')
    else
        subplot(1,3,2)
        hold on; shadedErrorBar(whisk_xy(:,1),whisk_xy(:,2),whisk_xy(:,3),'c--')
        subplot(1,3,3)
        hold on; plot(norm_whisk(:,1),norm_whisk(:,2),'c--')
    end
    
    %if touch cell...
    if strcmp(U{rec}.meta.touchProperties.responseType,'excited')
        responseWindow = U{rec}.meta.touchProperties.responseWindow;
        peak_y = ceil(max(touch_response)/10)*10;
        
        touch_tune = object_location_quantification(U,rec,hilbert_feature,'off');
        loc_xy = [touch_tune{rec}.stim_response.bars_fit.x' touch_tune{rec}.stim_response.bars_fit.mean...
            touch_tune{rec}.stim_response.bars_fit.confBands - touch_tune{rec}.stim_response.bars_fit.mean];
        norm_loc = [loc_xy(:,1) normalize_var(loc_xy(:,2),0,1)];
        
        subplot(1,3,1)
        hold on; plot(responseWindow,[peak_y peak_y],'r')
        set(gca,'ylim',[0 peak_y * 1.2])
        subplot(1,3,2)
        if touch_tune{rec}.is_tuned==1
            shadedErrorBar(loc_xy(:,1),loc_xy(:,2),loc_xy(:,3),'b')
        else 
            shadedErrorBar(loc_xy(:,1),loc_xy(:,2),loc_xy(:,3),'b--')
        end
        subplot(1,3,3)
        plot(norm_loc(:,1),norm_loc(:,2),'b')
    end
    
    title_labels = {'raw','norm'};
    for d = 2:3
        subplot(1,3,d)
        legend('whisk','touch')
        title([title_labels{d-1} ' whisk and touch tuning curve'])
        xlabel(hilbert_feature)
        ylabel('firing rate (spks/s)')
    end
    
    suptitle(['cell number : ' num2str(rec)])
end