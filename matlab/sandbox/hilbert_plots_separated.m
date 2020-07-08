touchWindow = [-25:50]; %window for analyses around touch
popV = touchFeatureBinned(U,touchWindow);
touchOrderFields = fields(popV{1});
whichTouches = 1;

% Defining touch response
U = defTouchResponse(U,.95,'off');

%%
for i = 9
    array = U{glmModel{i}.meta};
    
    if ~isfield(array.meta,'responseWindow')
        array.meta.responseWindow = [8 30];
    end
    
    BI = glmModel{i}.modelParams.buildIndices;
    touchResp = array.meta.responseWindow;
    spikes = squeeze(array.R_ntk(:,:,:));
    cellFR(i) = nanmean(spikes(:))*1000;
    
    gaussFilt = 2;
    sampledAngles = glmModel{i}.heat.angle;
    anglePlotSpacing = 3;
    tuningPrctile = 99;
    countsThresh = 10;
    
    
    selThetaBins = popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.range(popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.counts>=countsThresh);
    selectedBins = ismember(glmModel{i}.heat.angle,selThetaBins)    ;
    selRawMat = glmModel{i}.heat.matrixRaw([min(find(selectedBins)):max(find(selectedBins))],:);
    selMdlMat = glmModel{i}.heat.matrix([min(find(selectedBins)):max(find(selectedBins))],:);
    
    rawThetaIdx = ismember(popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.range,selThetaBins);
    rawTheta = popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.range(rawThetaIdx);
    rawMap= popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.spikes;
    figure(123)
    subplot(4,5,i)
    imagesc(imgaussfilt(rawMap(rawThetaIdx,:),gaussFilt,'padding','replicate'));
    hold on;plot([find(BI==0)-1 find(BI==0)-1],[size(glmModel{i}.heat.matrixRaw,1) 0],'w:')
    set(gca,'ytick',1:anglePlotSpacing:length(rawTheta),'yticklabel', rawTheta(1:anglePlotSpacing:length(rawTheta)),'ydir','normal',...
        'xtick',(0:25:length(BI)),'xticklabel',[min(BI):25:max(BI)],'xlim',[0 length(BI)])
    ylabel('angle at touch')
    
    hold on; plot([find(BI==touchResp(1)) find(BI==touchResp(1))],[0 max(selThetaBins)+1],'r-.')
    hold on; plot([find(BI==touchResp(2)) find(BI==touchResp(2))],[0 max(selThetaBins)+1],'r-.')
    axis square
    figure(124)
    subplot(4,5,i)
    imagesc(imgaussfilt(selRawMat,gaussFilt,'padding','replicate'));
    %     caxis([0 prctile(selRawMat(:),tuningPrctile)])
    hold on;plot([find(BI==0)-1 find(BI==0)-1],[size(glmModel{i}.heat.matrixRaw,1) 0],'w:')
    set(gca,'ytick',1:anglePlotSpacing:length(selThetaBins),'yticklabel', selThetaBins(1:anglePlotSpacing:length(selThetaBins)),'ydir','normal',...
        'xtick',(0:25:length(BI)),'xticklabel',[min(BI):25:max(BI)],'xlim',[0 length(BI)])
    ylabel('angle at touch')
    axis square
    
    hold on; plot([find(BI==touchResp(1)) find(BI==touchResp(1))],[0 max(selThetaBins)+1],'r-.')
    hold on; plot([find(BI==touchResp(2)) find(BI==touchResp(2))],[0 max(selThetaBins)+1],'r-.')
    
    figure(125)
    subplot(4,5,i)
    imagesc(imgaussfilt(selMdlMat,gaussFilt,'padding','replicate'));
    %     caxis([0 prctile(selMdlMat(:),tuningPrctile)]) %cap limit at 99th percentile in case of outliers in predictions
    hold on;plot([find(BI==0)-1 find(BI==0)-1],[size(glmModel{i}.heat.matrixRaw,1) 0],'w:')
    set(gca,'ytick',1:anglePlotSpacing:length(selThetaBins),'yticklabel', selThetaBins(1:anglePlotSpacing:length(selThetaBins)),'ydir','normal',...
        'xtick',(0:25:length(BI)),'xticklabel',[min(BI):25:max(BI)],'xlim',[0 length(BI)])
    %     xlabel('time from touch onset (ms)')
    ylabel('angle at touch')
    axis square
    hold on; plot([find(BI==touchResp(1)) find(BI==touchResp(1))],[0 max(selThetaBins)+1],'r-.')
    hold on; plot([find(BI==touchResp(2)) find(BI==touchResp(2))],[0 max(selThetaBins)+1],'r-.')
    
    
end

figure(123); suptitle('original tuning curve')
figure(124);   suptitle('training tuning curve')
figure(125);   suptitle('modeled tuning curve')

%%
comparison = struct();
%OBJECT LOCATION TUNING
figure(126);clf
otune=nan(size(glmModel,1));
tuningSharpness = cell(3,1);
thetabounds = [-100:2:100];

for i = 9
    array = U{glmModel{i}.meta};
    
    if ~isfield(array.meta,'responseWindow')
        array.meta.responseWindow = [8 30];
    end
    
    BI = glmModel{i}.modelParams.buildIndices;
    touchResp = array.meta.responseWindow;
    spikes = squeeze(array.R_ntk(:,:,:));
    cellFR(i) = nanmean(spikes(:))*1000;
    
    gaussFilt = 2;
    sampledAngles = glmModel{i}.heat.angle;
    anglePlotSpacing = 3;
    tuningPrctile = 99;
    countsThresh = 10;
    
    selThetaBins = popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.range(popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.counts>=countsThresh);
    
    colors = {'b','k','r'};
    tangles = [];
    for g = 1:length(popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.range)
        tangles = [tangles ; repmat(popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.range(g),popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.counts(g),1)];
    end
    
    spikesGroup = {glmModel{i}.raw.spikes, cell2mat(popV{glmModel{i}.meta}.(touchOrderFields{whichTouches}).theta.raw) , glmModel{i}.predicted.spikeProb};
    angleGroup = {glmModel{i}.raw.trimmedAngle, tangles, glmModel{i}.predicted.angles};
    namesGroup = {'original','trainingData','predicted'};
    
    figure(126);
    subplot(4,5,i)
    for d = 1:length(spikesGroup)
        
        allspks = spikesGroup{d};
        angleatT = angleGroup{d};
        
        bl = allspks(:,1:find(BI==0));
        response = allspks(:,find(BI==touchResp(1)):find(BI==touchResp(2)));
        
        
        meanbase =  mean(mean(bl,2));
        stdbase = std(mean(bl,2));
        meanresp = mean(response,2);
        zscore = (meanresp - meanbase) ./ stdbase;
        
        [sorted,~,~] = binslin(angleatT,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
        
        samps = cellfun(@numel,sorted);
        %         selBins = find(samps>sum(samps)./100);
        selBins = find(samps>10);
        
        
        zraw = nan(length(sorted),1000);
        for b = 1:length(sorted)
            if sum(b == selBins)>0
                currvals = sorted{b};
                if ~isempty(sorted{b})
                    zraw(b,1:length(currvals)) = currvals';
                end
            end
        end
        
        [p,~,stats]=anova1(zraw',[],'off');
        vals =multcompare(stats,[],'off');
        otune(i,d)=p;
        
        if  otune(i,d)<.01
            [~,maxResponseIdx] = max(nanmean(zraw,2));
            for j = 1:size(zraw,1)
                [~,pz(j)] = ttest2(zraw(maxResponseIdx,:),zraw(j,:));
            end
            tuningSharpness{d}(:,i) = pz;
        end
        
        x=cellfun(@str2num,stats.gnames);
        
        %CI BINS
        cibins = nan(size(x,1),1);
        
        SEMg = nanstd(zraw,[],2) ./ sqrt(sum(~isnan(zraw),2));
        for p = 1:length(x)
            rawx = zraw(x(p),:);
            SEM = SEMg(x(p));
            ts = tinv(.95,sum(~isnan(rawx),2)-1);      % T-Score
            cibins(p,:) = ts.*SEM;   %confidence intervals
        end
        
        SEM = SEMg(x);
        y = stats.means;
        x=thetabounds(x)';
        
        xy = [x y'];
        
        hold on;
        %         shadedErrorBar(xy(:,1),normalize_var(xy(:,2),0,1),cibins,['-' colors{d}]);
        shadedErrorBar(xy(:,1),xy(:,2),cibins,['-' colors{d}]);
        hold on;
        set(gca,'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
        
        zresp{d} = xy;
    end
    
    tmp = cell2mat(zresp');
    set(gca,'ylim',[min(tmp(:,2))*.5 max(tmp(:,2))*1.5])
    
    [raw,~ ] = intersect(zresp{2}(:,1),selThetaBins);
    [modeled,~ ] = intersect(zresp{3}(:,1),selThetaBins);
    
    interIdx = ismember(zresp{2}(:,1),modeled);
    interIdx2 = ismember(zresp{3}(:,1),raw);
    
    [zresp{2}(interIdx,1) zresp{3}(interIdx2,1)]
    
    expVar(i) = corr(zresp{2}(interIdx,2),zresp{3}(interIdx2,2)).^2;
    
    title([ 'expVar=' num2str(expVar(i) )])
end


tuningCompare.otune.columnNames = namesGroup;
tuningCompare.otune.values = otune;
tuningCompare.sharpness.cellNames = namesGroup;
tuningCompare.sharpness.value = tuningSharpness;
suptitle('object location tuning of raw(black) vs modeled (red)')


% confusionMat of tuning
tunedCells = tuningCompare.otune.values<.01;
peakTuning = nan(length(glmModel),length(tuningCompare.otune.columnNames));
for k = 1:3
    [~,peakTuning(:,k) ] = max(tuningCompare.sharpness.value{k});
end

peakTuning = thetabounds(peakTuning);
peakTuning(~tunedCells)=nan;

devExplained = cellfun(@(x) nanmean(x.gof.devExplained),glmModel);
deThresh = .08;
selCells = find(devExplained>deThresh);

figure(130);clf;subplot(1,2,1);
scatter(peakTuning(selCells,2),peakTuning(selCells,3))
set(gca,'xlim',[min(peakTuning(:)) max(peakTuning(:))],'ylim',[min(peakTuning(:)) max(peakTuning(:))])
hold on;plot([min(peakTuning(:)) max(peakTuning(:))],[min(peakTuning(:)) max(peakTuning(:))],'-.k')
axis square
xlabel('raw tuning')
ylabel('predicted tuning')

subplot(1,2,2);
histogram(peakTuning(selCells,3)-peakTuning(selCells,2),-5:2.5:20)
set(gca,'xlim',[-5 20])
axis square
xlabel('degrees from raw tuning');ylabel('number of cells')

%goodness of fit of ol tuning and variance explained
figure(131);clf
scatter(devExplained,expVar)
axis square
xlabel('devExplained');ylabel('rsqd of OL tuning')


%% goodness of fit!
devExplained = cellfun(@(x) nanmean(x.gof.devExplained),glmModel);
for i = 9
    array = U{glmModel{i}.meta};
    if ~isfield(array.meta,'responseWindow')
        array.meta.responseWindow = [8 30];
    end
    
    BI = glmModel{i}.modelParams.buildIndices;
    touchResp = array.meta.responseWindow;
    spikes = squeeze(array.R_ntk(:,:,:));
    cellFR(i) = nanmean(spikes(:))*1000;
    
    figure(127);
    subplot(4,5,i)
    hold on;plot(BI,smooth(mean(glmModel{i}.predicted.spikeTestRaw *1000,1),glmModel{i}.modelParams.bf.bfwidth),'k','linewidth',2)
    hold on;plot(BI,mean(glmModel{i}.predicted.spikeProb *1000,1),'r','linewidth',2)
    ylabel('firing rate (hz)')
    
    smoothingWindow = glmModel{i}.modelParams.bf.bfwidth*2+1;
    rawResponse = smooth(mean(glmModel{i}.predicted.spikeTestRaw *1000,1)',smoothingWindow);
    modeledResponse = mean(glmModel{i}.predicted.spikeProb *1000,1)';
    tr1 = find(BI == array.meta.responseWindow(1));
    tr2 = find(BI == array.meta.responseWindow(2));
    
    correlationFit(i) = corr(rawResponse(tr1:tr2),modeledResponse(tr1:tr2));
    
    title(['gof = ' num2str(devExplained(i)) ' or ' num2str(correlationFit(i))])
    set(gca,'xtick',-25:25:50)
end
suptitle('touch response of raw(black) vs modeled (red')

%% kernels
figure(128);clf
BI = glmModel{i}.modelParams.buildIndices;
for i = 1:length(glmModel)
    subplot(4,5,i)
    coeffsToPlot = fields(glmModel{i}.coeffs);
    for u = 2:length(coeffsToPlot)
        coeffsToPlotName = coeffsToPlot{u};
        if strcmp(coeffsToPlotName,'touch') || strcmp(coeffsToPlotName,'touchDur')
            hold on; plot(BI,sum(glmModel{i}.coeffs.(coeffsToPlotName)'.*glmModel{i}.basisFunctions.touch,2))
        else
            hold on; plot(BI,sum(glmModel{i}.coeffs.(coeffsToPlotName)'.*glmModel{i}.basisFunctions.features,2))
        end
    end
    set(gca,'xtick',[-25:25:50])
end
%     ylabel('kernel weight')
%     xlabel('time from touch onset (ms)')
legend(coeffsToPlot(2:end))
suptitle('kernel weights in spike prediction')

%% feature importance in spike prediction for best cells
devExplained = cellfun(@(x) nanmean(x.gof.devExplained),glmModel);
deThresh = .08;
clear kernels
for i = 1:length(glmModel)
    array = U{glmModel{i}.meta};
    BI = glmModel{i}.modelParams.buildIndices;
    touchResp = array.meta.responseWindow;
    
    tr1 = find(BI==touchResp(1));
    tr2 = find(BI==touchResp(2));
    
    
    coeffsToPlot = fields(glmModel{i}.coeffs);
    for u = 2:length(coeffsToPlot)
        coeffsToPlotName = coeffsToPlot{u};
        if strcmp(coeffsToPlotName,'touch') || strcmp(coeffsToPlotName,'touchDur')
            kernels{i}(u-1,:) = sum(glmModel{i}.coeffs.(coeffsToPlotName)'.*glmModel{i}.basisFunctions.touch,2);
        else
            kernels{i}(u-1,:)  = sum(glmModel{i}.coeffs.(coeffsToPlotName)'.*glmModel{i}.basisFunctions.features,2);
        end
    end
    
    trKernels = nansum(abs(kernels{i}(:,tr1:tr2)),2);
    
    wtProportion(i,:) = trKernels ./ nansum(trKernels);
    
end

selCells = wtProportion(find(devExplained>deThresh),:);

figure(230);clf
shadedErrorBar(1:7,nanmean(selCells),nanstd(selCells)./sqrt(numel(find(devExplained>deThresh))))
set(gca,'xticklabel',coeffsToPlot(2:end))
ylabel('proportion of weight distribution')