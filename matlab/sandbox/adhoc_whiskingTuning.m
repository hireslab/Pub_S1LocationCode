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
whisking = whisking_general(U,'off');
hilbertWhisking = whisking_hilbert(U,popV,'off');

%%
[rc]= numSubplots(length(U))

fieldsToCompare = hilbertWhisking.cellVarNames;

WHilbertP = cell(1,numel(fieldsToCompare)); 
peakResponse = cell(1,numel(fieldsToCompare)); 

pctSamp = .01; 
smoothParam = .7;
alpha = .01; 
for i = 4
    figure(230);clf
    featureMean = cellfun(@(x) nanmean(x),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    featureSEM = cellfun(@(x) nanstd(x) ./ sqrt(sum(~isnan(x))),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    sampleSize = cellfun(@(x) sum(~isnan(x)),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    tscore = cellfun(@(x) tinv(.0975,x-1), sampleSize,'uniformoutput',0);
    CI = cellfun(@(x,y) x.*y,featureSEM,tscore,'uniformoutput',0);
    

    selBins= cellfun(@(x) sum(~isnan(x))> round(sum(~isnan(x(:)))*pctSamp),hilbertWhisking.R_ntk.(fieldsToCompare{i}),'uniformoutput',0);
    x = hilbertWhisking.S_ctk.(fieldsToCompare{i});
    for d = 1:length(U)
        [WHilbertP{i}(d),~,~] = anova1(hilbertWhisking.R_ntk.(fieldsToCompare{i}){d},[],'off');

        subplot(rc(1),rc(2),d);
        if strcmp(fieldsToCompare{i},'phase')
            f = fit((1:numel(x))',(featureMean{d}.*1000)','smoothingspline','SmoothingParam',smoothParam);
            fCI = fit((1:numel(x))',(CI{d}.*1000)','smoothingspline','SmoothingParam',smoothParam);
            fitline = f(1:numel(x));
            fCIline = fCI(1:numel(x));
            
            %             bar(x,featureMean{d}.*1000,'facecolor',[.8 .8 .8])
            if WHilbertP{i}(d) < alpha
                hold on; shadedErrorBar(x,fitline,fCIline,'g')
            else
                hold on; shadedErrorBar(x,fitline,fCIline,'k')
            end
            set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
            
            [~,idx] = max(fitline);
            peakResponse{i}(d) = x(idx);
        else
            xtv = x(selBins{d});
            
            %             bar(xtv,featureMean{d}(selBins{d}).*1000,'facecolor',[.8 .8 .8])
            f = fit(xtv',(featureMean{d}(selBins{d}).*1000)','smoothingspline','SmoothingParam',smoothParam);
            fCI = fit(xtv',(CI{d}(selBins{d}).*1000)','smoothingspline','SmoothingParam',smoothParam);
            fitline = f(xtv(1):xtv(end));
            fCIline = fCI(xtv(1):xtv(end));
            if WHilbertP{i}(d) < alpha
                hold on; shadedErrorBar(xtv(1):xtv(end),fitline,fCIline,'g')
            else
                hold on; shadedErrorBar(xtv(1):xtv(end),fitline,fCIline,'k')
            end
            
            set(gca,'xlim',[min(xtv) max(xtv)],'ylim',[0 (max(featureMean{d}(selBins{d}).*1000)+.01).*1.2])
            
            [~,idx] = max(fitline);
            ranges = xtv(1):xtv(end);
            peakResponse{i}(d) = ranges(idx);
        end
    end
    
    suptitle(fieldsToCompare{i})
    
   print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} '_whisking'],'-dpng')
   print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} '_whisking'],'-depsc')
end

%%
wOn = [find(whisking.matrix(1,:)==1) find(whisking.matrix(1,:)==-1)];
pTune = cell2mat(WHilbertP');
chosenPTune = pTune(:,wOn)<.01; %looking only at cells that are whisking excited and inhibited. 

angleOnly = intersect(find(chosenPTune(1,:)),find(sum(chosenPTune)==1));

ampOnly = intersect(find(chosenPTune(2,:)),find(sum(chosenPTune)==1));

mpOnly = intersect(find(chosenPTune(3,:)),find(sum(chosenPTune)==1));

pOnly = intersect(find(chosenPTune(4,:)),find(sum(chosenPTune)==1));

%to find bar/spline of above do wOn(angleOnly); 
twofeats = find(sum(chosenPTune)==2);
threefeats = find(sum(chosenPTune)==3);
allfeats = find(sum(chosenPTune)==4);
unmod = find(sum(chosenPTune)==0);

plotSort = [unmod angleOnly ampOnly mpOnly pOnly angleOnly twofeats threefeats allfeats];

sortedMap = chosenPTune(:,plotSort);
figure(9);clf
imagesc(sortedMap')
set(gca,'xtick',1:4,'xticklabel',fieldsToCompare)



%% what degree are these tuned cells max firing at? 
peaks = cell2mat(peakResponse');
whisking_peaks = peaks(:,wOn);

tunedAmp = whisking_peaks(1,chosenPTune(1,:));
tunedMP = whisking_peaks(2,chosenPTune(2,:));
tunedPhase = whisking_peaks(3,chosenPTune(3,:));
tunedAngle = whisking_peaks(4,chosenPTune(4,:));

figure(53);clf
for i = 1:4
    tuning = whisking_peaks(i,chosenPTune(i,:));
    subplot(2,2,i);
    if i ==3
        histogram(tuning,-pi:pi/4:pi,'facecolor','k','facealpha',1)
    else
    histogram(tuning,min(tuning):4:max(tuning),'facecolor','k','facealpha',1)
    end
    
end
    
    
    

%% slow x fast tuning 
for g = 1:length(U)
    phase = hilbertWhisking.S_ctk.raw{g}(:,3);
    amp = hilbertWhisking.S_ctk.raw{g}(:,1);
    midpoint = hilbertWhisking.S_ctk.raw{g}(:,2); 
    angle = hilbertWhisking.S_ctk.raw{g}(:,4); 
    spikes = hilbertWhisking.R_ntk.raw{g}(:);
    spksT = find(spikes == 1);
    noSpksT = find(spikes == 0); 
    
    featX = phase; 
    featY = angle; 
    featZ = phase; 
    
    if ~isempty(spksT)
    figure(230);clf
    hold on; scatter(featX(noSpksT),featY(noSpksT),'ko','markeredgecolor',[.9 .9 .9])
    scatter(featX(spksT),featY(spksT),'ko','filled')
    title(num2str(corr(featX(spksT),featY(spksT),'rows','complete')));
    end

    pause
    
    
end
    








