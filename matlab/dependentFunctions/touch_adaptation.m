%% Figure 5 Tuning of neuronal responses to touch order, pretouch velocity
% and maximum curvature. It requires that variables created with the 
% S1datasetWrapper.m and intermediateAnalyasisStructures.mat are in the 
% main workspace.
% Tested on Matlab 2018b
% Adapted from Samuel Andrew Hires 2016-07-14
% Edited by Jonathan Cheung 2019-07-11

U = defTouchResponse(U,.99,'on');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)==1);

[rc] = numSubplots(length(selectedCells));
%%
for rec = 1:length(selectedCells)
    array=U{selectedCells(rec)}
    array.onsetLatency=8;
    trange=1:30;
    k_range = [1:20];
    % Adaptation across touches
    % Adaptation by ITI
    % Phase at touch
    % Velocity at touch
    % Peak force at touch
    % Direction of touch
    
    v(rec).spikes = squeeze(array.R_ntk);
    v(rec).allTouchIdx = find(nansum([array.S_ctk(9,:,:);array.S_ctk(12,:,:)]));
    v(rec).firstTouchIdx = find(array.S_ctk(9,:,:)==1);
    v(rec).lateTouchIdx = find(array.S_ctk(12,:,:)==1);
    v(rec).allTouchITI = diff([0; v(rec).allTouchIdx]);
    
    vel = [zeros(1,array.k); squeeze(diff( array.S_ctk(1,:,:)))];
    theta=squeeze(array.S_ctk(1,:,:));
    amp = squeeze(array.S_ctk(3,:,:));
    setp = squeeze(array.S_ctk(4,:,:));
    phase = squeeze(array.S_ctk(5,:,:));
    dk = squeeze(array.S_ctk(6,:,:));
    M0Adj = squeeze(array.S_ctk(7,:,:));
    Fax = squeeze(array.S_ctk(8,:,:));
    
    %trange = [1:resL45(rec).idxT];
%     tspikesIdx = repmat(v(rec).allTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(v(rec).allTouchIdx),1); %index for touch + trange
    tspikesIdx = v(rec).allTouchIdx + (array.meta.responseWindow(1):array.meta.responseWindow(end));
    
    if isnan(tspikesIdx)
        v(rec).sc = 0;
        v(rec).lsc = 0;
    else
        v(rec).sc = sum(v(rec).spikes(tspikesIdx),2);  
        lspikesIdx = repmat(v(rec).lateTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(v(rec).lateTouchIdx),1);
        v(rec).lsc = sum(v(rec).spikes(lspikesIdx),2);
    end
    
    % Touch Adaptation
    v(rec).touchNumber = [];
    for i = 1:length(v(rec).firstTouchIdx)-1;
        v(rec).touchNumber = cat(2,v(rec).touchNumber,1:sum(v(rec).allTouchIdx >= v(rec).firstTouchIdx(i) & v(rec).allTouchIdx < v(rec).firstTouchIdx(i+1)));
    end
    v(rec).touchNumber = cat(2,v(rec).touchNumber,1:sum(v(rec).allTouchIdx >= v(rec).firstTouchIdx(end)));
    v(rec).lateTouchNumber = v(rec).touchNumber(v(rec).touchNumber>1);
    v(rec).adaptation.lh = [];
    v(rec).adaptation.lci = [];
 
    [lh, lci] = poissfit(v(rec).sc(v(rec).touchNumber == 1));
    v(rec).adaptation.lh(1) = lh;
    v(rec).adaptation.lci(1,:) = lci;
    [v(rec).adaptation.sorted v(rec).adaptation.sortedBy v(rec).adaptation.binBounds]=binslin(v(rec).lateTouchNumber, v(rec).lsc, 'equalN',9);
    
    for num = 1:length(v(rec).adaptation.sorted)
        [lh, lci] = poissfit(v(rec).adaptation.sorted{num});
        v(rec).adaptation.lh(num+1) = lh;
        v(rec).adaptation.lci(num+1,:) = lci;
    end
    


    %% ADAPTATION
    [sort_adapt,idx_adapt] = sortrows ([v(rec).allTouchIdx, v(rec).touchNumber'],2);
    adapt_raster = zeros(length(idx_adapt),151);
    adapt_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_adapt),1,151)+repmat([-50:100],length(idx_adapt),1));
    
    figure(52);
    subplot(rc(1),rc(2),rec)
%     subplot(3,4,[1 3])
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(adapt_raster'),151), ceil(find(adapt_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    set(gca,'xtick',[],'ytick',[])
    colormap([1 1 1;0 0 0])
    box off
%     set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
    ax1 = gca;
    ax1;hold on
    ax1_pos = get(ax1,'Position');
    axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'yaxislocation','right','Color','none');hold on
    plot(mat2gray(sort_adapt(:,2)),1:length(sort_adapt),'r')
    axis tight
    norm_adapt_range = mat2gray([sort_adapt(1,2) sort_adapt(end,2)]);
    plot([norm_adapt_range(2)]*[1;1],[1 size(sort_adapt,1)],'--','color',[1 .5 .5])
%     set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_adapt(1,2) sort_adapt(end,2)]),'xticklabel',round([sort_adapt(1,2) (sort_adapt(end,2))]))
set(gca,'xtick',[],'ytick',[])
%     xlabel('Touch order in trial','color','r')
  
    figure(78300);
%     subplot(3,4,4)
    subplot(rc(1),rc(2),rec)
    plot(v(rec).adaptation.lh, 1:size(v(rec).adaptation.lh,2),'k-')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).adaptation.lci,1:size(v(rec).adaptation.lh,2), 'color',[.5 .5 .5])
    axes(ax1)
    set(gca,'color','none')
    %  print(gcf,'-depsc2',[printdir 'Raster_Adaptation_Cell_' num2str(rec)])
    
  % ADAPTATION BY ITI 
%    v(rec).allTouchITI;
%    v(rec).sc;
%    
%    [sorted,sortedBy,edges] = binslin(v(rec).allTouchITI,v(rec).sc,'equalE',22,0,500);
%    figure(480)
%    subplot(rc(1),rc(2),rec)
%    errorbar(linspace(0,500,21),cellfun(@mean,sorted),cellfun(@std,sorted),'o')
%    set(gca,'xlim',[-5 205],'xtick',[0:50:200],'ytick',[0:1:5])
%    
%     figure(45)
%     subplot(rc(1),rc(2),rec)
%     histogram(v(rec).allTouchITI,[0:25:500],'normalization','probability')
%     set(gca,'xtick',[0:100:500],'ytick',[0:.1:1])
%     xlabel('ITI');ylabel('Proportion of touches')
end



