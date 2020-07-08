function [adaptation] = adaptation_quantification(U,selectedCells,displayOpt)

if (nargin < 3), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    figure(80);clf
end

rc = numSubplots(numel(U));
%%
adaptation = cell(1,numel(U));
g_vy = cell(1,numel(selectedCells));
ITI_data = [];
for rec = 1:length(selectedCells)
    array=U{selectedCells(rec)};
    array.onsetLatency=array.meta.touchProperties.responseWindow(1);
    trange=1:diff(array.meta.touchProperties.responseWindow);
    % Adaptation across touches
    % Adaptation by ITI
    
    spikes = squeeze(array.R_ntk);
    allTouchIdx = find(nansum([array.S_ctk(9,:,:);array.S_ctk(12,:,:)]));
    firstTouchIdx = find(array.S_ctk(9,:,:)==1);
    lateTouchIdx = find(array.S_ctk(12,:,:)==1);
    allTouchITI = diff([0; allTouchIdx]);
    
    tspikesIdx = repmat(allTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(allTouchIdx),1); %index for touch + trange
    
    if isnan(tspikesIdx)
        sc = 0;
        lsc = 0;
    else
        sc = sum(spikes(tspikesIdx),2);
        lspikesIdx = repmat(lateTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(lateTouchIdx),1);
        lsc = sum(spikes(lspikesIdx),2);
    end
    
    %% Touch Adaptation by touch order
    touchNumber = [];
    for i = 1:length(firstTouchIdx)-1
        touchNumber = cat(2,touchNumber,1:sum(allTouchIdx >= firstTouchIdx(i) & allTouchIdx < firstTouchIdx(i+1)));
    end
    touchNumber = cat(2,touchNumber,1:sum(allTouchIdx >= firstTouchIdx(end)));
    lateTouchNumber = touchNumber(touchNumber>1);
    adaptation{selectedCells(rec)}.lh = [];
    adaptation{selectedCells(rec)}.lci = [];
    
    [lh, lci] = poissfit(sc(touchNumber == 1));
    adaptation{selectedCells(rec)}.lh(1) = lh;
    adaptation{selectedCells(rec)}.lci(1,:) = lci;
    [adaptation{selectedCells(rec)}.sorted, adaptation{selectedCells(rec)}.sortedBy, adaptation{selectedCells(rec)}.binBounds]=binslin(lateTouchNumber, lsc, 'equalN',9);
    
    for num = 1:length(adaptation{selectedCells(rec)}.sorted)
        [lh, lci] = poissfit(adaptation{selectedCells(rec)}.sorted{num});
        adaptation{selectedCells(rec)}.lh(num+1) = lh;
        adaptation{selectedCells(rec)}.lci(num+1,:) = lci;
    end
    
    [adaptation{selectedCells(rec)}.calculations.p,~,~]= anova1(cell2nanmat(adaptation{selectedCells(rec)}.sorted),[],'off');
    
    %multcompare(stats,[],'off')
    if willdisplay
        if adaptation{selectedCells(rec)}.calculations.p<.05
            
            subplot(rc(1),rc(2),selectedCells(rec))
            shadedErrorBar(1:size(adaptation{selectedCells(rec)}.lh,2),adaptation{selectedCells(rec)}.lh,adaptation{selectedCells(rec)}.lci(:,2) - adaptation{selectedCells(rec)}.lh','b')
        else
            subplot(rc(1),rc(2),selectedCells(rec))
            shadedErrorBar(1:size(adaptation{selectedCells(rec)}.lh,2),adaptation{selectedCells(rec)}.lh,adaptation{selectedCells(rec)}.lci(:,2) - adaptation{selectedCells(rec)}.lh','k')
        end
    end
    
%     maxResponse = max(adaptation{selectedCells(rec)}.lh);
%     minResponse = min(adaptation{selectedCells(rec)}.lh);
%     adaptation{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
    adaptation{selectedCells(rec)}.calculations.mod_idx_relative = mean(adaptation{selectedCells(rec)}.lh(2:end)) ./ adaptation{selectedCells(rec)}.lh(1);
    adaptation{selectedCells(rec)}.calculations.adapt_ratio =  (adaptation{selectedCells(rec)}.lh) ./ adaptation{selectedCells(rec)}.lh(1);
    adaptation{selectedCells(rec)}.sampleNumber_x_y = [1:10 ; histcounts(touchNumber,.5:10.5)];
    title(num2str(mean(adaptation{selectedCells(rec)}.lh(2:end)) ./ adaptation{selectedCells(rec)}.lh(1)));
    %% Touch Adaptation by ITI
    
    time_interp = 20:10:4000; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE THIS BEORE RUNNING. BUILT FOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING
    %     rec=datasample(1:length(selectedCells),1);

    array = U{(selectedCells(rec))};
    smooth_param = 5;
    num_touches_per_bin = 75;
%     
    spks = squeeze(array.R_ntk(:,:,:));
    response_window = array.meta.touchProperties.responseWindow(1):array.meta.touchProperties.responseWindow(2);
    touch_times = [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    vis_touch_response = spks(touch_times + (-25:50));
    calc_touch_response = spks(touch_times + response_window);
    
    numBins = round(numel(touch_times)./num_touches_per_bin);
    
    if numBins>=2
        ITI = abs(touch_times  - [nan(1,1) ; touch_times(1:end-1)]);
        ITI(1) = max(ITI);
        ITI(ITI>4000) = 4000;
        
        ITI_data.histogram_x = time_interp;
        ITI_data.histogram(:,rec) = histcounts(ITI,time_interp);
        
        [sorted_vis,sortedBy_vis] = binslin(ITI,vis_touch_response*1000,'equalN',numBins);
        [sorted_calc,~] = binslin(ITI,calc_touch_response*1000,'equalN',numBins);
        
        minx = cellfun(@min,sortedBy_vis);
        maxx = cellfun(@max,sortedBy_vis);
        medianx = cellfun(@median,sortedBy_vis);
        y = cellfun(@nanmean,sorted_vis,'uniformoutput',0) ;
        err = cellfun(@(x) nanstd(x) ./ sqrt(size(sorted_vis,1)),sorted_vis,'uniformoutput',0) ;
        
        mean_responses = cellfun(@(x) mean(sum(x,2)./numel(response_window)),sorted_calc);
        sem_responses = cellfun(@(x) std(sum(x,2)./numel(response_window)) ./ sqrt(size(x,1)),sorted_calc);
        
        
        [g_x,unique_idx] = unique(medianx);
        g_y = smooth(mean_responses(unique_idx),smooth_param);
        if numel(g_x>2)
            g_vy{rec} = interp1(g_x,g_y,time_interp);
        end
        
        
        
%             figure(8);clf
%             plot_elements = numBins-1;
%             for g = 1:plot_elements
%                 subplot(3,plot_elements+1,g+1)
%                 shadedErrorBar(-25:50,smooth(y{g},smooth_param),smooth(err{g},smooth_param))
%                 title([num2str(minx(g)) '-' num2str(maxx(g))])
%                 set(gca,'ylim',[0 max(cell2mat(y),[],'all')],'xlim',[min(response_window) max(response_window)])
%                 
%                 if g == 1
%                     subplot(3,plot_elements+1,1)
%                     shadedErrorBar(-25:50,smooth(y{end},smooth_param),smooth(err{end},smooth_param))
%                     title('"first touch"')
%                     set(gca,'ylim',[0 max(cell2mat(y),[],'all')],'xlim',[min(response_window) max(response_window)])
%                 end
%                 
%             end
%             
%             subplot(3,plot_elements+1,[(plot_elements+2) : (2*(plot_elements+1))])
%             errorbar(medianx,smooth(mean_responses,smooth_param),smooth(sem_responses,smooth_param),'vertical','ko-')
%             set(gca,'xlim',[0 max(medianx)],'xscale','log')
%             xlabel('ITI')
%             ylabel('response window fr')
%             suptitle(['cell num = ' num2str(selectedCells(rec))])
%             
%             subplot(3,plot_elements+1,[(plot_elements+plot_elements+3) : (3*(plot_elements+1))])
%             f = fit(ITI,sum(calc_touch_response,2),'smoothingspline','SmoothingParam',.000001);
%             hold on; plot(f,ITI,sum(calc_touch_response,2))
%             set(gca,'xlim',[0 max(medianx)],'xscale','log')
      

    end
    
    
    
    %     firing_rate = sum(spks(touch_times + response_window),2) ./ numel(response_window) .* 1000;
    %
    %     if willdisplay
    %         figure(30);
    %         subplot(rc(1),rc(2),rec);
    %         scatter(ITI,firing_rate,'k.')
    %         set(gca,'xlim',[0 prctile(ITI,99)])
    %     end
    %
    
    
    
    %% ADAPTATION plotting raster by touch order
    %     [sort_adapt,idx_adapt] = sortrows ([allTouchIdx, touchNumber'],2);
    %     adapt_raster = zeros(length(idx_adapt),151);
    %     adapt_raster = spikes(repmat(allTouchIdx(idx_adapt),1,151)+repmat([-50:100],length(idx_adapt),1));
    %
    %
    
    %         figure(52);subplot(3,4,[1 3])
    %     plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    %     hold on
    %     plot(mod(find(adapt_raster'),151), ceil(find(adapt_raster')/151) ,'k.','markersize',4)
    %     axis([1 151 1 size(adapt_raster,1)]);
    %     colormap([1 1 1;0 0 0])
    %     box off
    %     set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
    %     ax1 = gca;
    %     ax1;hold on
    %     ax1_pos = get(ax1,'Position');
    %     axes('Position',ax1_pos,...
    %         'XAxisLocation','top',...
    %         'yaxislocation','right','Color','none');hold on
    %     plot(mat2gray(sort_adapt(:,2)),1:length(sort_adapt),'r')
    %     axis tight
    %     norm_adapt_range = mat2gray([sort_adapt(1,2) sort_adapt(end,2)]);
    %     plot([norm_adapt_range(2)]*[1;1],[1 size(sort_adapt,1)],'--','color',[1 .5 .5])
    %     set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_adapt(1,2) sort_adapt(end,2)]),'xticklabel',round([sort_adapt(1,2) (sort_adapt(end,2))]))
    %     xlabel('Touch order in trial','color','r')
    %
    %     plot(adaptation.lh, 1:size(adaptation.lh,2),'k-')%(velocity.binBounds(1:end-1)+velocity.binBounds(2:end))/2, 'color','b')
    %     hold on
    %     box off
    %     set(gca,'ylim',[0.5 10.5],'ytick',[])
    %     plot(adaptation.lci,1:size(adaptation.lh,2), 'color',[.5 .5 .5])
    %     axes(ax1)
    %     set(gca,'color','none')
    
end


%POPULATION PLOTTING

if willdisplay
    %% TOUCH ORDER PLOTTING
    figure(9);clf
    subplot(2,1,1)
    sample_mat = cell2mat(cellfun(@(x) x.sampleNumber_x_y(2,:),adaptation(selectedCells),'uniformoutput',0)')';
    sampled_proportions = sample_mat./sum(sample_mat);
    cdf_proportions = cumsum(sample_mat)./sum(sample_mat);
    yyaxis left
    shadedErrorBar(1:10,mean(sampled_proportions,2),std(sampled_proportions,[],2) ./ sqrt(numel(selectedCells)),'b')
    ylabel('proportion of touches')
    yyaxis right
    shadedErrorBar(1:10,mean(cdf_proportions,2),std(cdf_proportions,[],2) ./ sqrt(numel(selectedCells)),'r')
    ylabel('CDF of touches')
    title('touch order sampling')
    set(gca,'xlim',[1 10])
    
    
    subplot(2,1,2)
    adapt_mat = cell2mat(cellfun(@(x) x.calculations.adapt_ratio,adaptation(selectedCells),'uniformoutput',0)')';
    adapt = find(mean(adapt_mat(2:end,:))<1);
    facil = find(mean(adapt_mat(2:end,:))>1);
    plot(1:10,adapt_mat,'color',[.8 .8 .8])
    hold on; shadedErrorBar(1:10,mean(adapt_mat(:,adapt),2),std(adapt_mat(:,adapt),[],2)./sqrt(numel(adapt)),'r')
    hold on; shadedErrorBar(1:10,mean(adapt_mat(:,facil),2),std(adapt_mat(:,facil),[],2)./sqrt(numel(facil)),'b')
    title(['facil(blue n=' num2str(numel(facil)) ') and adapt (red n=' num2str(numel(adapt)) ')'])
    ylabel('adaptation ratio (x/first touch response)')
    hold on; plot([1 10],[1 1],'--k')
    set(gca,'ylim',[0 3],'xlim',[1 10],'ytick',0:1:3)
    xlabel('touch order')
    
    %% ITI PLOTTING 
    
    ITI_mat = normalize_var(cell2mat(g_vy')',0,1);
    x = linspace(25,3995,398);
    xticks = 1:10:length(time_interp);
    xticklabel = time_interp(xticks);

    [~,maxidx] = max(ITI_mat,[],1);
    [~,idx] = sort(maxidx);
%     [~,idx] = sort(nansum(ITI_mat(1:20,:)));
%     [~,idx] = sort(ITI_mat(end,:)); 
   
    figure(8324);clf
    subplot(3,10,[2:10])
    %CDF OF SAMPLES
    cdf = cumsum(ITI_data.histogram) ./ sum(ITI_data.histogram);
    yyaxis right
    hold on; shadedErrorBar(x,nanmean(cdf,2),nanstd(cdf,[],2)./ sqrt(sum(~isnan(cdf),2)),'r-')
    set(gca,'xscale','log','xlim',[0 4000],'xtick',[25 50 100 500 1000 4000],'ylim',[0 1])
    xlabel('ITI')
    ylabel('CDF of touches')
    %HIST OF SAMPLES
    mean_histo_counts = nanmean(ITI_data.histogram ./ sum(ITI_data.histogram),2);
    sem_histo_counts = nanstd(ITI_data.histogram ./ sum(ITI_data.histogram),[],2) ./ sqrt(size(ITI_data.histogram,2));
    yyaxis left
    shadedErrorBar(x,mean_histo_counts,sem_histo_counts,'b')
    set(gca,'xscale','log','xlim',[0 4000],'xtick',[25 50 100 500 1000 4000],'ylim',[0 .1])
    title('touch ITI sampling')
    xlabel('ITI')
    ylabel('proportion of touches')
    
    %HEATMAP OF ITI RESPONSE
    subplot(3,10,11)
    imagesc(flipud(ITI_mat(end,idx)')); %had to flipUD to keep it consistent with pcolor which doesnt reverse
    ylabel('neurons sorted by time to peak response')
    set(gca,'xtick',1,'xticklabel','1st')
    caxis([0 1])
    subplot(3,10,[12:20])
    pcolor(ITI_mat(:,idx)')
    set(gca,'xtick',xticks,'xticklabel',xticklabel,'xlim',[1 50],'ytick',[])
    xlabel('touch ITI')
    caxis([0 1])
    colorbar
    
    %MEAN ITI RESPONSE x GROUPING
    subplot(3,10,21:30)
    facil = idx(1:19);
    adapt = idx(20:end); 
    shadedErrorBar(time_interp,nanmean(ITI_mat,2),(nanstd(ITI_mat,[],2)./size(ITI_mat,2)),'k')
    hold on; shadedErrorBar(time_interp,nanmean(ITI_mat(:,facil),2),(nanstd(ITI_mat(:,facil),[],2)./ sum(~isnan(ITI_mat(:,facil)),2) ),'b')
    hold on; shadedErrorBar(time_interp,nanmean(ITI_mat(:,adapt),2),(nanstd(ITI_mat(:,adapt),[],2)./sum(~isnan(ITI_mat(:,adapt)),2)),'r')
    title('gray:all, blue:last 19, red:first 41')
    set(gca,'xscale','log','xlim',[25 4000],'xtick',[25 50 100 500 1000 4000])
    xlabel('ITI')
    ylabel('normalized touch response')
    
    %look at this for cutoff (this is the time to peak response)
    %sort(maxidx)*10, cut off arbitrarily set at 250ms
    
    
    
    circ_mat = [ITI_mat ; ITI_mat];
    final_adapt = nan(size(mat_view,1),50);
    for g = 1:size(mat_view,1)
        adapt_raw = circ_mat(maxidx(g)+(0:399),g);
        filt_adapt = adapt_raw(~isnan(adapt_raw));
        final_adapt(g,1:50) = filt_adapt(1:50);
    end
    figure(9);clf
    imagesc(final_adapt(idx,:))
    caxis([.5 1])
    
    
end
