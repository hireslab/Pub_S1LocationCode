function whisking = whisking_general(array,displayOpt)
% 
% This function outputs a general whisking analysis of:
% 1) comparison of spiking based on whisking vs quiescent periods
% 2) comparison of spiking in protraction vs retraction in whisking cells
% 3) tuning across the population as whisking reaches the peak or trough of whisking cycles
% all of this is completed using a alpha value of 0.01
%
% results are packaged in a struct of with value names and matrix
%%

if (nargin < 2), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));


%TOP level params
pthresh = .01;

%% 1) comparison of spiking based on whisking vs quiescent periods
for i = 1:length(array)
    
    [objmask]= maskBuilder(array{i});
    spikes = squeeze(array{i}.R_ntk(:,:,:));
    filteredSpikes = spikes.*objmask.touch;
    whiskSpikes = filteredSpikes .* objmask.whisking;
    quietSpikes = filteredSpikes .* objmask.quiet;
    
    [~,p(i)] = ttest2(quietSpikes(:),whiskSpikes(:));
    
    nwfr(i) = nansum(quietSpikes(:))./numel(quietSpikes)*1000;
    wfr(i) = nansum(whiskSpikes(:))./numel(whiskSpikes)*1000;
    
end

if willdisplay
    figure(9);subplot(1,3,1)
end

nwc = nan(length(array),1);
wc = nan(length(array),1);
whiskingDefinition = nan(1,length(array));

for g = 1:length(array)
    if p(g)<pthresh
        if nwfr(g) > wfr(g)
            if willdisplay
                hold on;scatter(nwfr(g),wfr(g),[],'filled','r')
            end
            whiskingDefinition(g) = -1;
        elseif nwfr(g)<wfr(g)
            if willdisplay
                hold on;scatter(nwfr(g),wfr(g),[],'filled','b')
            end
            whiskingDefinition(g) = 1;
        end
    else
        if willdisplay
            hold on;scatter(nwfr(g),wfr(g),[],'filled','markerfacecolor',[.7 .7 .7])
        end
        whiskingDefinition(g) = 0;
    end
end

if willdisplay
    set(gca,'xlim',[0 30],'ylim',[0 30],'xtick',0:10:30,'ytick',0:10:30)
    hold on; plot([0 30],[0 30],'-.k')
    xlabel('no whisk spks/s');ylabel('whisk spks/s')
    axis square
end



%histogram of the above
% vals = [nwfr' wfr'];
% nwcells = vals(find(whiskingDefinition==-1),:);
% wcells = vals(find(whiskingDefinition==1),:);
% figure(8);clf;
% subplot(2,1,1);histogram(wcells(:,1),0:2.5:30,'facealpha',1,'facecolor','b')
% hold on; histogram(nwcells(:,1),0:2.5:30,'facealpha',1,'facecolor','r')
% set(gca,'ylim',[0 20],'ytick',0:20:20)
% ylabel('nw spks')
% subplot(2,1,2);histogram(wcells(:,2),0:2.5:30,'facealpha',1,'facecolor','b')
% hold on; histogram(nwcells(:,2),0:2.5:30,'facealpha',1,'facecolor','r')
% set(gca,'ylim',[0 20],'ytick',0:20:20)
% ylabel('w spks')


%prop tuned blocking out touch periods.
sum(abs(whiskingDefinition))./numel(whiskingDefinition);


%% 2) comparison of spiking in protraction vs retraction in whisking cells
wccells = find(whiskingDefinition==1);
rVSp = zeros(length(wccells),2);
directionalityTuning = nan(1,length(array));

if willdisplay
    figure(9);subplot(1,3,2)
end


for d = 1:length(wccells)
    curr = array{wccells(d)};
    masks = assist_touchmasks(curr);
    
    amps = squeeze(curr.S_ctk(3,:,:)).*masks.touch;
    phase = squeeze(curr.S_ctk(5,:,:)).*masks.touch;
    spks = squeeze(curr.R_ntk(1,:,:)).*masks.touch;
    
    retWidx = intersect(find(amps>5),find(phase>0));
    proWidx = intersect(find(amps>5),find(phase<0));
    
    [~,rppval(d)] = ttest2(spks(retWidx),spks(proWidx));
    
    rVSp(d,:) = [mean(spks(retWidx))*1000 mean(spks(proWidx))*1000];
    
    [~,rORp]  = max(rVSp(d,:));
    
    if rppval(d)<=pthresh
        if rORp ==2
            if willdisplay
                hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','c')
            end
            directionalityTuning(wccells(d)) = 1;
        elseif rORp ==1
            if willdisplay
                hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','m')
            end
            directionalityTuning(wccells(d)) = -1;
        end
    elseif rppval(d)>pthresh
        if willdisplay
            hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','markerfacecolor',[.7 .7 .7])
        end
        directionalityTuning(wccells(d)) = 0;
    end
    
end
if willdisplay
    hold on; plot([0 30],[0 30],'-.k')
    axis square
    xlabel('retraction spks/s'); ylabel('protraction spks/s')
    set(gca,'xtick',0:10:30,'ytick',0:10:30)
end

propWhiskDtuned = sum(rppval<=pthresh)./numel(rppval);
dirCells = (rppval<=pthresh);
rVSp(dirCells,:);

%% 3) tuning across the population as whisking reaches the peak or trough of whisking cycles
viewWindow = [-25:50];
amplitudeThreshold = 5;
[rc] = numSubplots(length(array));

for i = 1:length(array)
    
    [objmask]= maskBuilder(array{i});
    P = findPeaksTroughs(array{i},amplitudeThreshold);
    spikes = squeeze(array{i}.R_ntk(:,:,:));
    filteredSpikes = spikes.*objmask.touch;
    whiskSpikes = filteredSpikes .* objmask.whisking;
    quietSpikes = filteredSpikes .* objmask.quiet;
    
    peakWindow = find(P.peakIdx==1) + viewWindow;
    troughWindow = find(P.troughIdx==1) + viewWindow;
    
    
    %trimming for edge cases
    peakWindow(sum(peakWindow > array{i}.t * array{i}.k,2)>0,:) = [];
    peakWindow(sum(peakWindow <= 0,2)>0,:) = [];
    troughWindow(sum(troughWindow > array{i}.t * array{i}.k,2)>0,:) = [];
    troughWindow(sum(troughWindow <= 0)>0,:) = [];
    
    %building spike mats
    peakSpikes = filteredSpikes(peakWindow);
    troughSpikes = filteredSpikes(troughWindow);
    
    %filtering spike mats with nan values which are touch windows
    %(don't want to just nanmean because it'll give false view of decrease firing)
    peakSpikes(sum(isnan(peakSpikes),2)>0,:) = [];
    troughSpikes(sum(isnan(troughSpikes),2)>0,:) = [];
    
    %     figure(21);subplot(rc(1),rc(2),i)
    %     plot(viewWindow,nanmean(peakSpikes,1)*1000,'b')
    %     hold on; plot(viewWindow,nanmean(troughSpikes,1)*1000,'c')
    %
    % % %     shadedErrorBar(viewWindow,nanmean(peakSpikes,1)*1000,nanstd(peakSpikes,1)*1000./sqrt(size(peakWindow,1)),'b')
    % % %     hold on;shadedErrorBar(viewWindow,nanmean(troughSpikes,1)*1000,nanstd(troughSpikes,1)*1000./sqrt(size(troughWindow,1)),'c')
    %     hold on; plot([0 0],[0 max(nanmean(peakSpikes,1)*1000)],'-.k')
    % %
    
    peaks(:,i) = mean(peakSpikes);
    troughs(:,i) = mean(troughSpikes);
    
end

if willdisplay
    figure(9);subplot(1,3,3)
    % plot(viewWindow,normalize_var(peaks,0,1),'color',[.8 .8 .8])
    hold on; shadedErrorBar(viewWindow,nanmean(peaks(:,wccells),2)*1000,nanstd(peaks,[],2)*1000,'b')
    hold on; shadedErrorBar(viewWindow,nanmean(troughs(:,wccells),2)*1000,nanstd(troughs,[],2)*1000,'c')
    set(gca,'xlim',[min(viewWindow) max(viewWindow)],'xtick',-25:25:50)
    hold on; plot([0 0],[0 25],'-.k')
    title('spks/s from whisking peak(blue) and trough(cyan)')
    axis square
    xlabel('time from peak or trough')
    ylabel('spks/s')
end



whisking.rowFeatNames = {'on/off tuning', 'pro/ret tuning'};
whisking.matrix = [whiskingDefinition;directionalityTuning];
