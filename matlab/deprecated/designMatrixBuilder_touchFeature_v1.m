function [DmatX,DmatYnorm,DmatYshuff] = designMatrixBuilder_touchFeature(U,popV,selectedCells,variableFields,touchOrderFields,viewWindow,numInterpPts)

%% Interpolating all neurons to within specific bounds so we can use for population decoding
preDecisionTouches = preDecisionTouchMat(U);

for g = 1:numel(touchOrderFields)
    
    modelYtouch = nan(numInterpPts,length(selectedCells));
    modelYtouchstd = nan(numInterpPts,length(selectedCells));
    modelYavail = nan(numInterpPts,length(selectedCells));
    modelYavailstd = nan(numInterpPts,length(selectedCells));
    
    for d = 1:length(selectedCells)
        
        currCell = selectedCells(d);
        
        %Response in pole availWindow
        masks = maskBuilder(U{currCell});
        motors = normalize_var(U{currCell}.meta.motorPosition,-1,1);
        spikes = squeeze(U{currCell}.R_ntk);
        
        inside = double(isnan(masks.avail)); 
        inside(inside==0) = nan; 
%       availResponses = nanmean(spikes.*masks.availtoMeanOffset) * 1000;
        availResponses = nanmean(spikes.*inside) * 1000;
        [sortedAvail,~,~] = binslin(motors',availResponses','equalE',11,-1,1);
        
        availMean = cellfun(@nanmean,sortedAvail); 
        availSTD = cellfun(@nanstd,sortedAvail);
     
        modelYavail(:,d) = interp1(linspace(-.9,.9,10),availMean,linspace(-.9,.9,numInterpPts));
        modelYavailstd(:,d) = interp1(linspace(-.9,.9,10),availSTD,linspace(-.9,.9,numInterpPts));
 
        %null response outside pole availWindow
        pOnset = round(mean(U{currCell}.meta.poleOnset)*1000);
%         outside = masks.avail;
%         outsideResponses = nanmean(spikes.*outside) * 1000;
        outsideResponses = nanmean(spikes(1:pOnset,:))*1000; 
        [sortedOutside,~,~] = binslin(motors',outsideResponses','equalE',11,-1,1);
        
        outMean = cellfun(@nanmean,sortedOutside); 
        outSTD = cellfun(@nanstd,sortedOutside);
     
        modelYout(:,d) = interp1(linspace(-.9,.9,10),outMean,linspace(-.9,.9,numInterpPts));
        modelYoutstd(:,d) = interp1(linspace(-.9,.9,10),outSTD,linspace(-.9,.9,numInterpPts));

        [po(d)] = anova1(cell2nanmat(sortedOutside),[],'off');
        
%         %Touch ZSCORE response
%         counts = popV{currCell}.(touchOrderFields{g}).(variableFields).counts;
%         ranges = popV{currCell}.(touchOrderFields{g}).theta.range ;
%         thetavect = [];
%         for k = 1:length(ranges)
%             thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
%         end
%         
%         allspikes = cell2mat(popV{currCell}.(touchOrderFields{g}).theta.raw);
%         tOnset = find(viewWindow==0);
%         responseWindow = U{currCell}.meta.responseWindow;
%         meanbase = mean(mean(allspikes(:,1:tOnset),2));
%         stdbase = std(mean(allspikes(:,1:tOnset),2));
%         postresponses = allspikes(:,tOnset+responseWindow(1):tOnset+responseWindow(2));
%         
%         xresp = mean(postresponses,2);
%         zscore = (xresp - meanbase) ./ stdbase;
%         
%         [sorted,~,~] = binslin(thetavect,zscore,'equalE',numel(bounds),bounds(1),bounds(end));
%         
%         samps = cellfun(@numel,sorted);
%         selBins = samps>sum(samps)./100;
%         
%         zraw = cell2nanmat(sorted)';
%         zraw = zraw(selBins,:);
%         [~,~,stats]=anova1(zraw',[],'off');
% 
%         x = str2double(stats.gnames);
%         x = bounds(x)';
%         y = stats.means;
%         zstd = nanstd(zraw,[],2);
%         ystd = zstd(~isnan(zstd))';
%         
%         modelYtouch(:,d) = interp1(normalize_var(x,0,1),y,linspace(0,1,numInterpPts));
%         modelYtouchstd(:,d) = interp1(normalize_var(x,0,1),ystd,linspace(0,1,numInterpPts));
%         
        %Each touch raw touch response      
        tV = atTouch_sorter(U{currCell},viewWindow,preDecisionTouches{currCell});
        
        spikeTrainResponse = tV.allTouches.R_ntk(:,find(viewWindow==0)+U{currCell}.meta.responseWindow(1) : find(viewWindow==0)+U{currCell}.meta.responseWindow(2)); %capture responses within tuch response window
%       spikeTrainResponse = tV.allTouches.spikeMat(:,find(viewWindow==0)+5: find(viewWindow==0)+35); %use all responses from 5:35ms post touch
        meanResponse = nanmean(spikeTrainResponse,2)*1000;
        motors = normalize_var(tV.allTouches.S_ctk(:,7),-1,1);
        
        [sortedTouch,~,~] = binslin(motors,meanResponse,'equalE',11,-1,1);
       
        touchMean = cellfun(@nanmean,sortedTouch); 
        touchSTD = cellfun(@nanstd,sortedTouch);
       
        modelYrawTouch(:,d) = interp1(linspace(-.9,.9,10),touchMean,linspace(-.9,.9,numInterpPts));
        modelYrawTouchstd(:,d) = interp1(linspace(-.9,.9,10),touchSTD,linspace(-.9,.9,numInterpPts));
        
        [pT(d)] = anova1(cell2nanmat(sortedTouch),[],'off');
        
        %Mean of each all touches response 
        Umotors = unique(motors);
        allTouchMean = nan(length(Umotors),1);
        for q = 1:length(Umotors)
            allTouchMean(q) = mean(meanResponse(find(motors == Umotors(q))));
        end
        
        [sortedTouchMean,~,~] = binslin(Umotors,allTouchMean,'equalE',11,-1,1);
       
        alltouchMean = cellfun(@nanmean,sortedTouchMean); 
        alltouchSTD = cellfun(@nanstd,sortedTouchMean);
       
        modelYrawAllTouch(:,d) = interp1(linspace(-.9,.9,10),alltouchMean,linspace(-.9,.9,numInterpPts));
        modelYrawAllTouchstd(:,d) = interp1(linspace(-.9,.9,10),alltouchSTD,linspace(-.9,.9,numInterpPts));
       
        [pTM(d)] = anova1(cell2nanmat(sortedTouchMean),[],'off');
    end
    
%     xNames = {'outside_pole','pole_avail','touch_responseZ','touch_responseR','trialMean_touch_responseR'};
%     Xs = {modelYout,modelYavail,modelYtouch,modelYrawTouch,modelYrawAllTouch};
%     XsSTD = {modelYoutstd,modelYavailstd,modelYtouchstd,modelYrawTouchstd,modelYrawAllTouchstd}; 
%     
    xNames = {'outside_pole','pole_avail','touch_responseR','trialMean_touch_responseR'};
    Xs = {modelYout,modelYavail,modelYrawTouch,modelYrawAllTouch};
    XsSTD = {modelYoutstd,modelYavailstd,modelYrawTouchstd,modelYrawAllTouchstd}; 
    
    for f = 1:length(Xs)
        resampNum = 500;
        resampX = nan(size(Xs{f},1)*resampNum,size(Xs{f},2));
        for i = 1:resampNum
%             resampX(size(Xs{f},1)*(i-1)+1:size(Xs{f},1)*(i-1)+size(Xs{f},1),:) =  normrnd(Xs{f},XsSTD{f});
                resampX(size(Xs{f},1)*(i-1)+1:size(Xs{f},1)*(i-1)+size(Xs{f},1),:) =  poissrnd(Xs{f});
        end
        
        DmatX.(xNames{f}) = resampX;
    end
    
    DmatYnorm = repmat([1:size(modelYtouch,1)]',resampNum,1);
    randshuff = randperm(length(DmatYnorm));
    DmatYshuff = DmatYnorm(randshuff);
end


% colors = jet(length(Xs));
% for b = 1:length(Xs)
%     raw = (Xs{b}-mean(Xs{b})) ./ std(Xs{b});
%     for g = 1:size(raw,2)
%     figure(23);subplot(4,8,g)
%     hold on;plot(linspace(-1,1,20),raw(:,g),'color',colors(b,:))
%     end
% end
% 
% legend('out','in','perTouch','trialTouches')




%     %     for d = 1:length(DmatY)
%     for d = 1
%         rando = randperm(length(DmatX));
%         tmpDmatX=DmatX(rando,:);
%         tmpDmatY=DmatY{d}(rando,:);
%         txp = [];
%
%         cvglm
%
%
%         [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY{d})),0);
%         [pred,opt_thresh,prob]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
%             ,tmpDmatY(end*.7:end,:),'Max');
%         txp = [txp ; tmpDmatY(end*.7:end) pred];
%
%         cmat = confusionmat(txp(:,1),txp(:,2));
%         ncmat = cmat./sum(cmat);
%         figure(9+d);clf;
%         imagesc(ncmat);
%         caxis([0 .40])
%         axis square
%         colorbar
%         title(['chance = ' num2str(1/size(modely,1))])
%         set(gca,'xtick',[],'ytick',[])
%         xlabel('truth');ylabel('predicted')
%
%         predDiff(:,d) = txp(:,2) - txp(:,1);
%     end
%     model{b}=histcounts(predDiff(:,1),-16.5:16.5,'normalization','probability');
%     stats(b,:) = [model{b}(17) std(predDiff)];
%     %     shuff{b}=histcounts(predDiff(:,2),-16.5:16.5,'normalization','probability');
%     figure(40);clf
%     plot(-16:16,model{b},'k')
%     %     hold on; plot(-16:16,shuff{b},'color',[.7 .7 .7])
%
%     set(gca,'xtick',-16:8:16,'ytick',0:.1:.5,'ylim',[0 .35])
%     %     decodeResolution = std(predDiff)* (mean(rangeExplored)./size(modely,1))
