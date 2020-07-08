%% Traces of features and examples
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));

cellNumber = 32; % good SNR cells 55, 57

[tVar] = atTouch_sorter(U{cellNumber},-25:50);
touch_response = nanmean(tVar.allTouches.R_ntk) * 1000;

try
    figure(480);clf
    bar(-25:50,touch_response,'k')
    title('touch response')
    xlabel('time from touch onset (ms)')
    ylabel('firing rate (spks/s)')
    
    responseWindow = U{cellNumber}.meta.touchProperties.responseWindow;
    peak_y = ceil(max(touch_response)/10)*10;
    hold on; plot(responseWindow,[peak_y peak_y],'r')
    set(gca,'ylim',[0 peak_y * 1.2])
catch 
    bar(-25:50,touch_response,'k')
    title('NOT TOUCH CELL!')
end


uTrialNumber = datasample(1:U{cellNumber}.k,1); %57/92/52/83/34 on cell 52 is good
trialNumber = U{cellNumber}.meta.usedTrialNums(uTrialNumber);

for uTrialNumber = 59
    trialNumber = U{cellNumber}.meta.usedTrialNums(uTrialNumber);
    
    
    touchOn = [find(U{cellNumber}.S_ctk(9,:,trialNumber)==1) find(U{cellNumber}.S_ctk(12,:,trialNumber)==1)];
    touchOff = [find(U{cellNumber}.S_ctk(10,:,trialNumber)==1) find(U{cellNumber}.S_ctk(13,:,trialNumber)==1)];
    
    angle = U{cellNumber}.S_ctk(1,:,trialNumber);
    amp = U{cellNumber}.S_ctk(3,:,trialNumber);
    midpoint = U{cellNumber}.S_ctk(4,:,trialNumber);
    phase = U{cellNumber}.S_ctk(5,:,trialNumber);
    phase(smooth(amp,50)<5) = nan;
    licks= find(U{cellNumber}.S_ctk(16,:,trialNumber)==1);
    sweep_name = ['sweepArray_' U{cellNumber}.meta.mouseName '_' U{cellNumber}.meta.sessionName '_' U{cellNumber}.meta.cellName '_' U{cellNumber}.meta.cellCode ];
    tarray_name = ['trial_array_' U{cellNumber}.meta.mouseName '_' U{cellNumber}.meta.sessionName '_' U{cellNumber}.meta.cellName '_' U{cellNumber}.meta.cellCode ];
    try
        load(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays\' sweep_name]);
    catch
        cd('C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays_Phil')
        file_names = dir('*.mat');
        for g = 1:length(file_names)
            matched_idx(g) = ~isempty(strfind(file_names(g).name, sweep_name));
        end
        load(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays_Phil\' file_names(find(matched_idx,1)).name]);
    end
    
    try
        load(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\TrialArrayBuilders\' tarray_name]);
    catch
        cd('C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\TrialArrayBuilders_Phil')
        file_names = dir('*.mat');
        for g = 1:length(file_names)
            matched_idx(g) = ~isempty(strfind(file_names(g).name, tarray_name));
        end
        load(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\TrialArrayBuilders_Phil\' file_names(find(matched_idx,1)).name]);
    end
    
    spikes = s.sweeps{trialNumber}.rawSignal - mean(s.sweeps{trialNumber}.rawSignal);
    spikeTimes = round(T.trials{trialNumber}.spikesTrial.spikeTimes./10)-150;
    uTimes = find(U{cellNumber}.R_ntk(:,:,uTrialNumber))';
    

    figure(8);clf
    subplot(6,1,1)
    plot(1:U{cellNumber}.t,angle);
    ylabel('angle')
    set(gca,'xtick',[])
    subplot(6,1,2)
    plot(1:U{cellNumber}.t,amp);
    ylabel('amp')
    set(gca,'xtick',[])
    subplot(6,1,3)
    plot(1:U{cellNumber}.t,midpoint);
    ylabel('midpoint')
    set(gca,'xtick',[])
    subplot(6,1,4)
    scatter(1:U{cellNumber}.t,phase,'.');
    ylabel('phase')
    set(gca,'xtick',[])
    subplot(6,1,5)
    if ~isempty(touchOn)
        for b = 1:length(touchOn)
            hold on; plot([touchOn(b) touchOff(b)],[1 1],'k')
        end
    end
    if ~isempty(spikeTimes)
        hold on; scatter(spikeTimes,ones(1,numel(spikeTimes)).*2,'k.')
    end
    
    if ~isempty(uTimes)
        hold on; scatter(uTimes,ones(1,numel(uTimes)).*3,'k.')
    end
    
     if ~isempty(uTimes)
        hold on; scatter(licks,ones(1,numel(licks)).*4,'c.')
    end
    
    set(gca,'xlim',[0 4000],'ylim',[0 4],'ytick',[1 2 3 4 ],'yticklabel',{'touch','spike times','Uspikes','licks'})
    subplot(6,1,6)
    
    plot(spikes)
    if sum(spikeTimes(end)-uTimes(end)) ~= 0
        spikeTimes = round(T.trials{trialNumber}.spikesTrial.spikeTimes./10)-200;
        uTimes = find(U{cellNumber}.R_ntk(:,:,uTrialNumber))';
        set(gca,'xtick',2000:10000:42000,'xticklabel',0:1000:4000,'xlim',[2000 42000])
    else
        set(gca,'xtick',1500:10000:41500,'xticklabel',0:1000:4000,'xlim',[1500 41500])
    end

    suptitle(['cell num:' num2str(cellNumber) ' trialNum:' num2str(uTrialNumber)]);
    pause
end

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
fn = 'example_trial_signals.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% firing rate X depth of recording
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,touchCells,'pole','off');
location_cells = touchCells(cellfun(@(x) x.is_tuned==1,pole_tuned(touchCells)));

jc_silent_cell = [766 819 895 631 776 815 910 871 844 902   941   840   888   748   732   940   686   944   950   933]; %Phils  from 902

figure(480);clf
subplot(3,1,[1 2])
scatter( cellfun(@(y) y.meta.depth,U),cellfun(@(y) mean(y.R_ntk(:))*1000, U),'k','filled')
hold on;scatter( cellfun(@(y) y.meta.depth,U(location_cells)),cellfun(@(y) mean(y.R_ntk(:))*1000, U(location_cells)),'g','filled')
hold on; scatter(jc_silent_cell,ones(length(jc_silent_cell),1)*-1,[],'c','filled');
hold on; plot([700 700],[0 30],'--k')
hold on; plot([900 900],[0 30],'--k')
title(['n = response (' num2str(numel(U)) ') + silent (' num2str(numel(jc_silent_cell)) ')' ])
set(gca,'ytick',0:10:30,'xtick',600:100:1000,'xlim',[600 1000])
ylabel('firing rate (hz)');
xlabel('depth (um)')


hold on; subplot(3,1,3)
histogram(cellfun(@(y) y.meta.depth,U(setdiff(1:length(U),location_cells))),600:25:1000,'facecolor','k')
hold on; histogram(cellfun(@(y) y.meta.depth,U(location_cells)),600:25:1000,'facecolor','g')
hold on;histogram(jc_silent_cell,600:25:1000,'facecolor','c')
set(gca,'xtick',600:100:1000)

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
fn = 'scatter_depth_firingrate.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Raster subplot(1) is unsorted and subplot(2) is sorted
for i = 29
    spikes = squeeze(U{i}.R_ntk);
    motors = U{i}.meta.motorPosition;
    [~,sidx] = sort(motors);
    sidx = fliplr(sidx);
    
    figure(38);clf
    for g = 1:length(motors)
        %         subplot(1,2,1)
        %         spikeIdx = find(spikes(:,g));
        %         hold on; scatter(spikeIdx,ones(numel(spikeIdx),1).*g,'k.')
        %         set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
        %         axis square
        %         subplot(1,2,2)
        %         spikeIdx = find(spikes(:,sidx(g)));
        %         hold on; scatter(spikeIdx,ones(numel(spikeIdx),1).*g,'k.')
    end
    %     set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
    %     axis square
    
    touchOn = [find(U{i}.S_ctk(9,:,:)==1)  ;find(U{i}.S_ctk(12,:,:)==1)];
    touchOff = [find(U{i}.S_ctk(10,:,:)==1)  ;find(U{i}.S_ctk(13,:,:)==1)];
    whisking = find(U{i}.S_ctk(3,:,:)>5);
    touch_matrix = nan(size(spikes));
    whisk_matrix = nan(size(spikes));
    for g = 1:length(touchOn)
        touch_matrix(touchOn(g):touchOff(g)) = 1;
    end
    whisk_matrix(setdiff(whisking,find(touch_matrix==1))) = 1;
    
    subplot(1,2,1)
    pcolor(whisk_matrix')
    set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000],'ylim',[48 94])
    axis square
    %SORTED TOUCH OVERLAY
    %     subplot(1,2,2)
    %     pcolor(touch_matrix(:,sidx)')
    %set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
    %
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
    fn = 'example_raster_whisk_2.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
end
%% PSTH sorted by go and nogo
for i = 28
    smoothTime = 100;
    
    spikes = squeeze(U{i}.R_ntk);
    go = find(U{i}.meta.trialType==1);
    nogo = find(U{i}.meta.trialType==0);
    
    figure(80);clf
    shadedErrorBar(1:U{i}.t,smooth(nanmean(spikes(:,go) ,2)*1000,smoothTime),smooth( (nanstd(spikes(:,go),[] ,2)*1000) ./ numel(go) ,smoothTime),'b')
    hold on;shadedErrorBar(1:U{i}.t,smooth(nanmean(spikes(:,nogo) ,2)*1000,smoothTime),smooth( (nanstd(spikes(:,nogo),[] ,2)*1000) ./ numel(nogo),smoothTime),'r')
    
    set(gca,'xtick',0:1000:4000,'ytick',0:10:40)
end



