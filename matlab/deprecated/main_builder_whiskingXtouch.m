%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%% Top level parameters and definitions
touchWindow = [-25:50]; %window for analyses around touch

touchCells = touchCell(U,'off');
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99,'off');

%% Quantifying object location tuning
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
hilbertTouch = tuningQuantification(U,popV,selectedCells,fieldsList([ 1 3 4 5]),whichTouches,touchWindow,'off');

%% Quantifying whisking tuning
whisking = whisking_general(U,'off');
hilbertWhisking = whisking_hilbert(U,popV,'off');

%% modulation index 
% not sure if this is useful. Followed the calculation that Curtis and
% Kleinfeld 2009 used. (max(response) - min(response) ./ mean(response)
selectedCells = find(touchCells==1); %look at only touch cells because those are the only ones w/ OL tuning
wt_modulationIndex(hilbertTouch,hilbertWhisking,selectedCells)

%% comparison between whisking and touch for touchCells 
fieldsToCompare = fields(hilbertTouch.R_ntk.allTouches);
selectedCells = find(touchCells==1); %look at only touch cells because those are the only ones w/ OL tuning
pThresh = 0.05; 
THilbertTuningP = cell(1,length(fieldsToCompare));
Ttune = cell(1,length(fieldsToCompare)); 
WHilbertTuningP = cell(1,length(fieldsToCompare));
WHtune = cell(1,length(fieldsToCompare)); 
Tresolution = cell(1,length(fieldsToCompare)); 
Wresolution = cell(1,length(fieldsToCompare)); 

for g = 1:4
    figure(23);clf;
    currTArray = hilbertTouch.R_ntk.allTouches.(fieldsToCompare{g});
    currWArray = hilbertWhisking.R_ntk.(fieldsToCompare{g});
    stimulus = hilbertWhisking.S_ctk.(fieldsToCompare{g});
    spacing = abs(stimulus(1)-stimulus(2));
    
    chosenTArrays = currTArray(selectedCells);
%     TouchR = cellfun(@(x) nanmean(x,2),chosenTArrays,'uniformoutput',0);
    TouchR = cellfun(@(x) normalize_var(nanmean(x,2),0,1),chosenTArrays,'uniformoutput',0);
    TouchSEM = cellfun(@(x) nanstd(x,[],2) ./ sqrt(sum(~isnan(x),2)),chosenTArrays,'uniformoutput',0);
    chosenWArrays = currWArray(selectedCells);
%     WhiskR = cellfun(@(x) nanmean(x,1),chosenWArrays,'uniformoutput',0);
    WhiskR = cellfun(@(x) normalize_var(nanmean(x,1),0,1),chosenWArrays,'uniformoutput',0);
    WhiskSEM = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),chosenWArrays,'uniformoutput',0);
   
    for i = 1:length(WhiskR)
        [THilbertTuningP{g}(i),~,~] = anova1(chosenTArrays{i}',[],'off');
        [WHilbertTuningP{g}(i),~,stats] = anova1(chosenWArrays{i},[],'off');
        
        subplot(4,8,i)
        %analysis for touch tuning that is significant
        if THilbertTuningP{g}(i)<pThresh
             shadedErrorBar(stimulus,TouchR{i}*1000,TouchSEM{i}*1000,'b')
             [~,tuneIdx] = max(TouchR{i});
             Ttune{g}(i)  = stimulus(tuneIdx);
             
             %resolution for touch tuning
             tPvals = nan(1,size(chosenTArrays{i},1)); 
             for k = 1:size(chosenTArrays{i},1)
             [~,tPvals(k)] = ttest2(chosenTArrays{i}(tuneIdx,:),chosenTArrays{i}(k,:));
             end
             if sum(tPvals<pThresh)>0
                Tresolution{g}(i) = min(abs(find(tPvals<pThresh) -find(tPvals==1))) * spacing;
             else
                Tresolution{g}(i) = nan;
             end
             
             %resolution for whisk tuning
             [~,tuneIdx] = max(WhiskR{i});
             WHtune{g}(i)  = stimulus(tuneIdx);
             
             wPvals = nan(1,size(chosenWArrays{i},2)); 
             for k = 1:size(chosenWArrays{i},2)
             [~,wPvals(k)] = ttest2(chosenWArrays{i}(:,tuneIdx),chosenWArrays{i}(:,k));
             end
             if sum(wPvals<pThresh)>0
                Wresolution{g}(i) = min(abs(find(wPvals<pThresh) -find(wPvals==1))) * spacing;
             else 
                 Wresolution{g}(i) = nan;
             end

        else
             shadedErrorBar(stimulus,TouchR{i}*1000,TouchSEM{i}*1000,'k')
             Ttune{g}(i)  = nan;
             Tresolution{g}(i) = nan;
             Wresolution{g}(i) = nan;
             WHtune{g}(i) = nan;
        end
        
        
        %analysis for whisking tuning that is significant
        if WHilbertTuningP{g}(i)<pThresh
             hold on; shadedErrorBar(stimulus,WhiskR{i}*1000,WhiskSEM{i}*1000,'r')

        else
             hold on; shadedErrorBar(stimulus,WhiskR{i}*1000,WhiskSEM{i}*1000,'k')
             WHtune{g}(i)  = nan;
        end 
        xBounds = stimulus(~isnan(TouchR{i}*1000));
        set(gca,'xlim',[min(xBounds) max(xBounds)])
        set(gca,'ytick',[])
        
        if strcmp(fieldsToCompare{g},'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
        end

    end
    suptitle(fieldsToCompare{g})
    
    print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{g} '_normalized'],'-dpng')

%     figure(30);
%     subplot(2,2,g)
%     allMat = [(cell2mat(WhiskR)*1000)' cell2mat(TouchR')*1000];
%     allMat = allMat(~sum(isnan(allMat),2),:); 
%     scatter(allMat(:,1),allMat(:,2),'k')
%     minMax(1) = min(allMat(:));
%     minMax(2) = prctile(allMat(:),99);
%     hold on; plot(minMax,minMax,'-.k')
%     set(gca,'xlim',minMax,'ylim',minMax)
%     axis square
%     title(fieldsToCompare{g})
%     xlabel('whisking FR');ylabel('touch FR')
    
end

figure(9);clf
for i = 1:4
    HilbVals = {WHtune{i},Ttune{i},Wresolution{i},Tresolution{i}};
    tmpVals = cell2mat(cellfun(@(x) ~isnan(x) .* (1:size(x,2)),HilbVals,'uniformoutput',0)');
    keepIdx = find((1:length(tmpVals)) - mean(tmpVals) == 0);
    
    subplot(2,2,i)
    for d = 1:length(keepIdx)
        hold on; errorbar(HilbVals{1}(keepIdx(d)),HilbVals{2}(keepIdx(d)),HilbVals{3}(keepIdx(d)),'horizontal','ko')
        hold on; errorbar(HilbVals{1}(keepIdx(d)),HilbVals{2}(keepIdx(d)),HilbVals{4}(keepIdx(d)),'vertical','ko')
    end
    
    if i ~= 4
        [minmax(1), minmax(2)] = bounds([WHtune{i} Ttune{i}]);
        set(gca,'xlim',minmax,'ylim',minmax)
        hold on;plot(minmax,minmax,'-.k')
        
    else
        set(gca,'xtick',-pi:pi/2:pi,'ytick',-pi:pi/2:pi,'ylim',[-pi pi],'xlim',[-pi pi]...
            ,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
        hold on; plot([-pi pi],[-pi pi],'-.k')
    end
    axis square
    title([fieldsToCompare{i} ' n=' num2str(numel(keepIdx)) ' units'])
    xlabel('whisking pref');ylabel('touch pref')
end

% 
% %% comparison
% naiveVSexpert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);
% 
% OLexpert = numel(intersect(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1),find(naiveVSexpert==1))) ./ nansum(hilbertTouch.populationQuant.theta.matrix(1,:) == 1);
% OLnaive = numel(intersect(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1),find(naiveVSexpert==0))) ./ nansum(hilbertTouch.populationQuant.theta.matrix(1,:) == 1);
% figure(100);clf
% subplot(2,2,1);pie([OLexpert OLnaive],{'expert','naive'})
% title(['OL cells = ' num2str(nansum(hilbertTouch.populationQuant.theta.matrix(1,:) == 1))])
% 
% wEXCOL = numel(intersect(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1),find(whisking.matrix(1,:) == 1))) ./ numel(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1)) ;
% wINHOL = numel(intersect(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1),find(whisking.matrix(1,:) == -1))) ./ numel(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1));
% nsOL = numel(intersect(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1),find(whisking.matrix(1,:) == 0))) ./ numel(find(hilbertTouch.populationQuant.theta.matrix(1,:) == 1));
% 
% subplot(2,2,1);pie([wEXCOL wINHOL nsOL],{'whisking EXC','whisking INH','whisking ns'})
% title(['OL cells = ' num2str(nansum(hilbertTouch.populationQuant.theta.matrix(1,:) == 1))])
% 
% wEXCnanOL = numel(intersect(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:))),find(whisking.matrix(1,:) == 1))) ./ numel(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:)))) ;
% wINHnanOL = numel(intersect(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:))),find(whisking.matrix(1,:) == -1))) ./ numel(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:))));
% nsnanOL = numel(intersect(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:))),find(whisking.matrix(1,:) == 0))) ./ numel(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:))));
% 
% subplot(2,2,2);pie([wEXCnanOL wINHnanOL nsnanOL],{'whisking EXC','whisking INH','whisking ns'})
% title(['nonOL cells = ' num2str(numel(find(isnan(hilbertTouch.populationQuant.theta.matrix(1,:)))) )])
% 
% 
% wEXCexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == 1))) ./ numel(find(whisking.matrix(1,:) == 1));
% wEXCnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == 1))) ./ numel(find(whisking.matrix(1,:) == 1));
% subplot(2,2,3);pie([wEXCexpert wEXCnaive],{'expert','naive'})
% title(['whisking EXC cells = ' num2str(numel(find(whisking.matrix(1,:) == 1)))])
% 
% 
% wINHexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == -1))) ./ numel(find(whisking.matrix(1,:) == -1));
% wINHnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == -1))) ./ numel(find(whisking.matrix(1,:) == -1));
% subplot(2,2,4);pie([wINHexpert wINHnaive],{'expert','naive'})
% title(['whisking INH cells = ' num2str(numel(find(whisking.matrix(1,:) == -1)))])
% 
% 
% nsexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == 0))) ./ sum(naiveVSexpert==1);
% nsnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == 0))) ./ sum(naiveVSexpert==0);
% 
% 
% %% all neuron comparison 
% fieldsToCompare = fields(hilbertTouch.R_ntk.allTouches);
% for g = 1:4
%     currArrays = hilbertWhisking.R_ntk.(fieldsToCompare{g});
%     for i = 1:length(U)
%         [FullWhiskTuning{g}(i),~,~] = anova1(currArrays{i},[],'off');
%     end
% end
% 
% hilbertWhiskTuning = cell2mat(FullWhiskTuning')<pThresh;
% 
% OLtuning = [nan(4,sum(~(touchCells==1))) (cell2mat(THilbertTuning')<pThresh)];
% [~,idx] = sort(touchCells);
% comparisonMat = [touchCells(idx) ;whisking.matrix(:,idx) ; hilbertWhiskTuning(:,idx) ; OLtuning];
% 
% %add choice decoding too?

