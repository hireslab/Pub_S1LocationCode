function [DmatX,DmatYnorm] = designMatrixBuilder_touchFeature(U,selectedCells,touchOrderFields,viewWindow,numInterpPts)

% v2.0 of decoding design matrix. Instead of interpolating, this function
% bootstraps data in each sampled motor bin then choosing a resampNum to
% grab from each bin to use as the input design matrix.
preDecisionTouches = preDecisionTouchMat(U);

trialNumThreshold = 90; %min number of trials to consider using that cell 
numTperCell = cellfun(@(x) x.k,U(selectedCells));
resampNum = round(mean(numTperCell(numTperCell>trialNumThreshold)./numInterpPts));

for g = 1:numel(touchOrderFields)
    
    for d = 1:length(selectedCells)
        
        currCell = selectedCells(d);
        
        if U{currCell}.k>trialNumThreshold
            
            masks = maskBuilder(U{currCell});
            motors = normalize_var(U{currCell}.meta.motorPosition,-1,1);
            spikes = squeeze(U{currCell}.R_ntk);
            %         spikes = reshape(spikes(randperm(numel(spikes))),size(spikes));
            %         spikes = zeros(size(spikes));
            %         percTot = .01;
            %         spikes(datasample(1:numel(spikes), round(percTot*numel(spikes)))) = 1;
            
            %null response outside pole availWindow
            pOnset = round(mean(U{currCell}.meta.poleOnset)*1000);
            outsideResponses = nanmean(spikes(1:pOnset,:))*1000;
            [~,bootsam] = bootstrp(50,@corr,outsideResponses,motors);
            %         [~,bootsam] = bootstrp(50,@mean,outsideResponses);
            resampMotors = motors(bootsam);
            resampOutside = outsideResponses(bootsam);
            [sortedOutside,~,~] = binslin(resampMotors(:),resampOutside(:),'equalE',numInterpPts+1,-1,1);
            sortedOutside = cellfun(@(x) x(randperm(numel(x))),sortedOutside,'uniformoutput',0); %shuffling data in case of resampling bias
            resampedOutside = cell2nanmat(sortedOutside);
            selOutside = resampedOutside(1:resampNum,:);
            
            DmatX.outside_pole(:,d) = selOutside(:);
            DmatX.outside_poleShuffled(:,d) = selOutside(randperm(length(selOutside(:))))';
            
            %Response in pole availWindow
            inside = double(isnan(masks.avail));
            inside(inside==0) = nan;
            insideResponses = nanmean(spikes.*inside) * 1000;
            [~,bootsam] = bootstrp(50,@corr,insideResponses,motors);
            resampMotors = motors(bootsam);
            resampInside = insideResponses(bootsam);
            [sortedInside,~,~] = binslin(resampMotors(:),resampInside(:),'equalE',numInterpPts+1,-1,1);
            
            resampedInside = cell2nanmat(sortedInside);
            selInside = resampedInside(1:resampNum,:);
            
            DmatX.pole_avail(:,d) = selInside(:);
            
            [sortedInsideP,~,~] = binslin(motors(:),insideResponses(:),'equalE',numInterpPts+1,-1,1);
            resampedInsideP = cell2nanmat(sortedInsideP);
            [pI(d)] = anova1(resampedInsideP,[],'off');
            
            %Each touch raw touch response
            tV = atTouch_sorter(U{currCell},viewWindow,preDecisionTouches{currCell});
            
            if isfield(U{currCell}.meta,'responseWindow')
                spikeTrainResponse = tV.allTouches.R_ntk(:,find(viewWindow==0)+U{currCell}.meta.responseWindow(1) : find(viewWindow==0)+U{currCell}.meta.responseWindow(2)); %capture responses within tuch response window
            else
                spikeTrainResponse = tV.allTouches.R_ntk(:,5 : 35); %capture responses within tuch response window
            end
            
            meanResponse = nanmean(spikeTrainResponse,2)*1000;
            motors = normalize_var(tV.allTouches.S_ctk(:,7),-1,1);
            
            [~,bootsam] = bootstrp(50,@corr,meanResponse,motors);
            resampMotors = motors(bootsam);
            resampTouchResponse = meanResponse(bootsam);
            [sortedTouch,~,~] = binslin(resampMotors(:),resampTouchResponse(:),'equalE',numInterpPts+1,-1,1);
            
            resampedTouch = cell2nanmat(sortedTouch);
            selTouches = resampedTouch(1:resampNum,:);
            
            DmatX.touch_responseR(:,d) = selTouches(:);
            
            %Mean of each all touches response
            Umotors = unique(motors);
            allTouchMean = nan(length(Umotors),1);
            allTouchSum = nan(length(Umotors),1);
            for q = 1:length(Umotors)
                allTouchMean(q) = nanmean(meanResponse(find(motors == Umotors(q))));
                allTouchSum(q) =  nansum(spikeTrainResponse(find(motors == Umotors(q)),:),'all');
            end
            [~,bootsam] = bootstrp(50,@corr,allTouchMean,Umotors);
            resampMotors = Umotors(bootsam);
            resampTouches = allTouchMean(bootsam);
            resampSum = allTouchSum(bootsam); 
            [sortedTouchMean,~,~] = binslin(resampMotors(:),resampTouches(:),'equalE',numInterpPts+1,-1,1);
            [sortedTouchSum,~,~] = binslin(resampMotors(:),resampSum(:),'equalE',numInterpPts+1,-1,1);
            
            resampedTouchMean = cell2nanmat(sortedTouchMean);
            selTouchesTrial = resampedTouchMean(1:resampNum,:);
            resampedTouchAll = cell2nanmat(sortedTouchSum);
            selTouchesAllTrial = resampedTouchAll(1:resampNum,:);
            
            DmatX.trialMean_touch_responseR(:,d) = selTouchesTrial(:);
            DmatX.trialSum_touch_responseR(:,d) = selTouchesAllTrial(:);
            
        else
            disp(['skipping cell ' num2str(currCell) ' because num trials <' num2str(trialNumThreshold)])
            
            DmatX.outside_pole(:,d) = nan;
            DmatX.outside_poleShuffled(:,d)= nan;
            DmatX.pole_avail(:,d) = nan;
            DmatX.touch_responseR(:,d)= nan;
            DmatX.trialMean_touch_responseR(:,d)= nan;
        end
        
    end
    
    Ymat = repmat(1:numInterpPts,resampNum,1);
    DmatYnorm = Ymat(:);
end

%
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


