
selectedCells = find(is_tuned==1);


%function parameters
viewWindow = [-25:50]; %viewing window around touch 
numTouchesPerBin = 75; %number of touches to assign in each bin for quantification. 
min_bins = 5; %minimum number of angle bins to consider quantifying 


%dependent function to id all touches and pre/post decision touches
preDecisionTouches = preDecisionTouchMat(uberarray);

timePoints =2:2:30;
modulationDepth = nan(length(selectedCells),length(timePoints)); 
touchDurations = cell(length(selectedCells),1);

for rec = 1:length(selectedCells)
    %stimulus and response variables definitions
    array = uberarray{selectedCells(rec)};
    [tVar] = atTouch_sorter(array,viewWindow,preDecisionTouches{selectedCells(rec)});
    
    if ~isempty(hilbert_feature)
        if strcmp(hilbert_feature,'angle')
            selected_feature = tVar.allTouches.S_ctk(:,1);
        elseif strcmp(hilbert_feature,'amplitude')
            selected_feature = tVar.allTouches.S_ctk(:,3);
        elseif strcmp(hilbert_feature,'midpoint')
            selected_feature = tVar.allTouches.S_ctk(:,4);
        elseif strcmp(hilbert_feature,'phase')
            selected_feature = tVar.allTouches.S_ctk(:,5);
        elseif strcmp(hilbert_feature,'curvature')
            selected_feature = tVar.allTouches.S_ctk(:,6);
        else
            error('select features of "angle", "amplitude", "midpoint", or "phase"')
        end
    else
        error('select features of "angle", "amplitude", "midpoint", or "phase"')
    end
    
    touchDurations{rec} = tVar.allTouches.dtS_ctk(:,1);
    
    for b = 2:length(timePoints)-1
        rw = find(viewWindow == 0) : find(viewWindow == timePoints(b));
        response = mean(tVar.allTouches.R_ntk(:,rw),2) * 1000;
        numBins = round(numel(selected_feature)./numTouchesPerBin);
        
        % Tuning in touch response window
        if strcmp(hilbert_feature,'phase')
            [sorted, sortedBy] = binslin(selected_feature,response,'equalE',13,-pi,pi);
        else
            [sorted, sortedBy] = binslin(selected_feature,response,'equalN',numBins);
        end
        
            
        quant_ol_p = anova1(cell2nanmat(sorted),[],'off');
        
        SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
        tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
        CI = SEM.*tscore;
        
        if numel(sortedBy)>min_bins && quant_ol_p<.01
            modulationDepth(rec,b) = (max(cellfun(@mean,sorted)) - min(cellfun(@mean,sorted))) ./ mean(cellfun(@mean,sorted));
        else 
            if b>=2 %this is to populate values in between tuning w/ 0s so we can interpolate
                if ~isnan(modulationDepth(rec,b-1))
                    modulationDepth(rec,b)=0;
                else
                    modulationDepth(rec,b) = nan;
                end
            end
        end
        

    end
    
end

%Interpolating modulation depth for 0'd values in between tuning
for i = 1:size(modulationDepth,1)
    if sum(modulationDepth(i,:)==0)>0
        tmp1 = modulationDepth(i,:);
        tmp1(isnan(tmp1))=[];
        nonzero_idx = find(tmp1>0);
        linear_interp = interp1(nonzero_idx,tmp1(nonzero_idx),1:length(tmp1));
        modulationDepth(i,~isnan(modulationDepth(i,:))) = linear_interp;
    end
end

%filling all nan values w/ 0s
unmod_neurons = sum(isnan(modulationDepth),2)==length(timePoints);
modulationDepth(isnan(modulationDepth))=0;
modulationDepth(unmod_neurons,:) = nan;

%% plotting features

%plotting paramaters
ms_bins_before = 3; %this is the ms before peak modulation to plot

%plot individual traces normalized and shifted to peak 
figure(8);clf
subplot(2,1,1)
shiftedModulation = cell(length(selectedCells),1);
peakModIdx = nan(length(selectedCells),1);
for g = 1:length(selectedCells)
    plotVals = normalize_var(modulationDepth(g,:),0,1);
    [~,peakModIdx(g)] = max(plotVals);
    if peakModIdx(g)>ms_bins_before
    hold on; plot(-ms_bins_before:-ms_bins_before+length(plotVals(peakModIdx(g)-ms_bins_before:end))-1,plotVals(peakModIdx(g)-ms_bins_before:end),'color',[.8 .8 .8])
    shiftedModulation{g} = plotVals(peakModIdx(g)-ms_bins_before:end);
    end
     
end
%plot population mean w/ SEM traces
peak_OL_integ_wdw = timePoints(peakModIdx(~unmod_neurons));
modMat = cell2nanmat(shiftedModulation'); 
meanMod = nanmean(modMat,2); 
semMod = nanstd(modMat,[],2)./sqrt(sum(~isnan(modMat),2)); 
hold on; shadedErrorBar(-ms_bins_before:-ms_bins_before+length(meanMod)-1,meanMod,semMod,'k'); 
set(gca,'xtick',-ms_bins_before:ms_bins_before:10,'xticklabel',ms_bins_before*-2:ms_bins_before*2:20,'ytick',0:.25:1,'xlim',[-ms_bins_before 9])
xlabel('time from peak modulation (ms)')
ylabel('normalized modulation index')
title(['mean time for peak modulation = ' num2str(mean(peak_OL_integ_wdw)) 'ms'])


subplot(2,1,2);
td_hist = cellfun(@(x) histcounts(x,-2.5:5:52.5) ./ numel(x),touchDurations,'uniformoutput',0);
bar(0:5:50,mean(cell2mat(td_hist)),'k')
set(gca,'xtick',0:10:50)
xlabel('duration of touch (ms)')
ylabel('proportion of touches')
title(['pop. median touch duration = ' num2str(mean(cellfun(@median,touchDurations))) 'ms'])

%plot individual raw traces
figure(9);clf
plot(timePoints,modulationDepth',0,1)
xlabel('integration time from touch onset (ms)')
ylabel('modulation depth')
title(['mean time for peak modulation = ' num2str(mean(peak_OL_integ_wdw)) 'ms'])

%slope of rise from previous ms bin to peak mod idx. 
slope = modMat(4,:)-modMat(3,:) ;