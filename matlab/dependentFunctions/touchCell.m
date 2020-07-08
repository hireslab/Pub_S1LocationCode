%Input U array and std threshold above baseline to classify as touch cell
%Finds whether 5:30ms post touch is greater than the mean 50ms pre touch.
%Future edit could involve just grabbing peak post touch within that window
%and using that to define

function [touchORnaw] = touchCell(U,displayOpt)

if (nargin < 2), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

for rec=1:length(U)
    window = [-50:150];
    touchIdx= [find(U{rec}.S_ctk(9,:,:)==1) ; find(U{rec}.S_ctk(12,:,:)==1)];
    spks = squeeze(U{rec}.R_ntk);
    
    within_range = ~ (logical(sum((touchIdx + window) > numel(spks),2)) | logical(sum((touchIdx + window) < 0, 2)));
    touchIdx = touchIdx(within_range); 
    
    touchSpks = spks(repmat(touchIdx,1,length(window))+repmat(window,size(touchIdx,1),1));
    sTrains(rec,:) = mean(touchSpks);
    blclose = mean(touchSpks(:,36:50),2);
    prclose = mean(touchSpks(:,52:66),2);
    
    blmid = mean(touchSpks(:,26:50),2);
    prmid = mean(touchSpks(:,52:76),2);
    
    blfar =  mean(touchSpks(:,1:25),2);
    prfar = mean(touchSpks(:,77:101),2);
    
    
    blpr(1,:) = [mean(blclose) mean(prclose)];
    blpr(2,:)=[mean(blmid) mean(prmid)];
    blpr(3,:)=[mean(blfar) mean(prfar)];
    
    trClose = ttest2(blclose,prclose,'alpha',.01); 
    trMid = ttest2(blmid,prmid,'alpha',.01);
    trFar = ttest2(blfar,prfar,'alpha',.01);
    trtest = [trClose trMid trFar];
    trtest(isnan(trtest)) = 0;
    ttanalysis = logical(trtest');
    
    if ismember(1,ttanalysis)
        sigs = find(ttanalysis==1);
        diffs = blpr(sigs,2)-blpr(sigs,1) ;
        
        if mean(diffs>0)
            touchORnaw(rec) = 1;% TOUCH ON
        elseif  mean(diffs<0)
            touchORnaw(rec) = -1; %TOUCH OFF
        else
            touchORnaw(rec) = 0; %NONTOUCH
        end
    else
        touchORnaw(rec) = 0; 
    end
    
end

%Display option for histogram and heatmpa plotting
if (willdisplay)
    figure(40);clf;
    rowColumns = numSubplots(length(U));
    plotrow=rowColumns(1);
    plotcol=rowColumns(2);
    
    for k = 1:size(sTrains,1)
        subplot(plotrow,plotcol,k);
        avgSpikes = sTrains(k,:);
        if touchORnaw(k)==1
            bar(-50:150,avgSpikes,'b');
        elseif touchORnaw(k) == -1
            bar(-50:150,avgSpikes,'r');
        else
            b=bar(-50:150,avgSpikes,'k');
            b.FaceColor = [.5 .5 .5];
        end
        hold on; plot([0 0],[0 max(sTrains(k,:))*1.5],'-.k','linewidth',2)
        set(gca,'xlim',[-25 50],'xtick',0,'ytick',[]);
        
    end
    
    %% HEATMAP PLOTTING
    tc=find(touchORnaw==1);
    heatSpks = nan(length(tc),76);
    for g = 1:sum(touchORnaw==1)
        avgSpikes = sTrains(tc(g),:);
        heatSpks(g,:) = normalize_var(avgSpikes(25:100),0,1);
        excitRaw(g,:) = avgSpikes(25:100);
    end
    [~,idx]  = max(heatSpks');
    [~,ftidx] = sort(idx);
    
    ntc = find(touchORnaw==0);
    for g = 1:sum(touchORnaw==0)
        avgSpikes = sTrains(ntc(g),:);
        ntheatSpks(g,:) = normalize_var(avgSpikes(25:100),0,1);
        ntRaw(g,:) = avgSpikes(25:100);
    end
    
    INHIBtc=find(touchORnaw==-1);
    INHIBheatSpks = nan(length(INHIBtc),76);
    for g = 1:sum(touchORnaw==-1)
        avgSpikes = sTrains(INHIBtc(g),:);
        INHIBheatSpks(g,:) = normalize_var(avgSpikes(25:100),0,1);
        inhibRaw(g,:) = avgSpikes(25:100);
    end
    [~,idx]  = min(INHIBheatSpks(:,25:end)');
    [~,INHIBftidx] = sort(idx);
    
    if ~isempty(tc)
        figure(48);clf;
        subplot(4,1,[1 2])
        imagesc(heatSpks(ftidx,:))
        hold on; plot([26 26],[0 length(tc)],'-.w','linewidth',3)
        colorbar
        colormap(jet)
        caxis([0 1])
        set(gca,'xtick',1:25:76,'xticklabel',-25:25:75,'ytick',[])
        
        figure(580);clf
        subplot(3,1,1)
        bar(-25:50,mean(excitRaw)*1000,'facecolor',[.5 .5 .5])
        set(gca,'xlim',[-25 50],'xtick',-25:25:50)
    end
    
    if ~isempty(INHIBtc)
        figure(48);hold on;
        subplot(4,1,[3])
        imagesc(INHIBheatSpks(INHIBftidx,:))
        hold on; plot([26 26],[0 length(tc)],'-.w','linewidth',3)
        colorbar
        colormap(jet)
        caxis([0 1])
        set(gca,'xtick',125:76,'xticklabel',-25:25:75,'ytick',[])
        
        figure(580);subplot(3,1,2)
        bar(-25:50,mean(inhibRaw)*1000,'facecolor',[.5 .5 .5])
        set(gca,'xlim',[-25 50],'xtick',-25:25:50)
    end
    
    if ~isempty(ntc)
        subplot(4,1,[4])
        imagesc(ntheatSpks)
        hold on; plot([26 26],[0 length(ntc)],'-.w','linewidth',3)
        colorbar
        colormap(jet)
        caxis([0 1])
        set(gca,'xtick',1:25:76,'xticklabel',-25:25:75,'ytick',[])
        
        
        figure(580);subplot(3,1,3)
        bar(-25:50,mean(ntRaw)*1000,'facecolor',[.5 .5 .5])
        set(gca,'xlim',[-25 50],'xtick',-25:25:50)
    end
end





