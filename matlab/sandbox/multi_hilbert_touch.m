%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% 1 dimensional hilbert features
touchWindow = [-25:50]; %window for analyses around touch

touchCells = touchCell(U,'off');
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99,'off');

% Quantifying hilbert touch tuning
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
hilbertTouch = tuningQuantification(U,popV,selectedCells,fieldsList([ 1 3 4 5]),whichTouches,touchWindow,'off');

%%
[rc]= numSubplots(length(U))

fieldsToCompare = fields(hilbertTouch.R_ntk.allTouches);

THilbertP = cell(1,numel(fieldsToCompare)); 
peakResponse = cell(1,numel(fieldsToCompare)); 

pctSamp = .01; 
smoothParam = .7;
alpha = .01; 
for i = 1:4
    figure(230);clf
    featureMean = cellfun(@(x) nanmean(x,2),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    featureSEM = cellfun(@(x) nanstd(x,[],2) ./ sqrt(sum(~isnan(x),2)),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    sampleSize = cellfun(@(x) sum(~isnan(x),2),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    tscore = cellfun(@(x) tinv(.0975,x-1), sampleSize,'uniformoutput',0);
    CI = cellfun(@(x,y) x.*y,featureSEM,tscore,'uniformoutput',0);
    

    selBins= cellfun(@(x) sum(~isnan(x),2)> round(sum(~isnan(x(:)))*pctSamp),hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}),'uniformoutput',0);
    x = hilbertTouch.S_ctk.(fieldsToCompare{i});
    for d = 1:length(U)
        if ~isempty(featureMean{d})
            
            [THilbertP{i}(d),~,~] = anova1(hilbertTouch.R_ntk.allTouches.(fieldsToCompare{i}){d}',[],'off');
            
            subplot(rc(1),rc(2),d);
            if strcmp(fieldsToCompare{i},'phase')
                x = linspace(-pi,pi,12);
                xtv = x(selBins{d})';
                
                %bar(x,featureMean{d}.*1000,'facecolor',[.8 .8 .8])
                f = fit(xtv,(featureMean{d}(selBins{d}).*1000),'smoothingspline','SmoothingParam',smoothParam);
                fCI = fit(xtv,(CI{d}(selBins{d}).*1000),'smoothingspline','SmoothingParam',smoothParam);
                fitline = f(x);
                fCIline = fCI(x);
                
                if THilbertP{i}(d)< alpha
                    hold on; shadedErrorBar(x,fitline,fCIline,'g')
                else
                    hold on; shadedErrorBar(x,fitline,fCIline,'k')
                end
                
                set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
                [~,idx] = max(fitline);
                peakResponse{i}(d) = x(idx);
            else
                xtv = x(selBins{d})';
                
%               bar(xtv,featureMean{d}(selBins{d}).*1000,'facecolor',[.8 .8 .8])
                f = fit(xtv,(featureMean{d}(selBins{d}).*1000),'smoothingspline','SmoothingParam',smoothParam);
                fCI = fit(xtv,(CI{d}(selBins{d}).*1000),'smoothingspline','SmoothingParam',smoothParam);
                fitline = f(xtv(1):xtv(end));
                fCIline = fCI(xtv(1):xtv(end));
                
                if THilbertP{i}(d)< alpha
                    hold on; shadedErrorBar(xtv(1):xtv(end),fitline,fCIline,'g')
                else
                    hold on; shadedErrorBar(xtv(1):xtv(end),fitline,fCIline,'k')
                end
                
                set(gca,'xlim',[min(xtv) max(xtv)],'ylim',[0 (max(featureMean{d}(selBins{d}).*1000)+.01).*1.2])
                [~,idx] = max(fitline);
                ranges = xtv(1):xtv(end);
                peakResponse{i}(d) = ranges(idx);
            end
        else 
            THilbertP{i}(d) = nan; 
        end
    end
    
    suptitle(fieldsToCompare{i})
    
%    print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} 'touch'],'-dpng')
%    print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{i} 'touch'],'-depsc')
end

OLcells = THilbertP{1}<alpha;


%% %2dimensional or 3dimensional scatter of features at touch
% Defining touch response
U = defTouchResponse(U,.99,'off');

% viewing window around touch 
touchWindow = [-25:50];

clear hilbertTouch 

preDecisionTouches = preDecisionTouchMat(U);
for rec = 1:length(U)
 [tVar] = atTouch_sorter(U{rec},touchWindow,preDecisionTouches{rec});

 if isfield(U{rec}.meta,'responseWindow')
     tResponse = U{rec}.meta.responseWindow+find(touchWindow==0) ;
%      tResponse = [-25 -5] + find(touchWindow==0); %control to find off touch responses
 else
     tResponse = [5 35] + find(touchWindow==0);
 end
  response = mean(tVar.allTouches.R_ntk(:,tResponse(1):tResponse(2)),2);
 
  hilbertTouch.S_ctk{rec} = tVar.allTouches.S_ctk(:,[3, 4, 5, 1]);
  hilbertTouch.R_ntk{rec} = response; 
  
end

figure(2332);clf
figure(123);clf
figure(124);clf
for k = 1:length(U)
    phase = hilbertTouch.S_ctk{k}(:,3);
    amp = hilbertTouch.S_ctk{k}(:,1);
    midpoint = hilbertTouch.S_ctk{k}(:,2);
    angle = hilbertTouch.S_ctk{k}(:,4);
    
    responses = normalize_var(hilbertTouch.R_ntk{k},.05,.95);
    rawResponses = hilbertTouch.R_ntk{k};
    [~,idx] = sort(responses);
    
    featX = phase;
    featY = angle;
    featZ = angle;
    
    colors = repmat(1 - responses,1,3);
    
    if isfield(U{k}.meta,'responseWindow')
        figure(2332);subplot(6,10,k)
        scatter(featX(idx),featY(idx),50,colors(idx,:),'filled');
        set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
        
        protraction = phase<=0;
        retraction = phase>0;
        
        bounds = popV{1}.allTouches.theta.bounds;
        [pSorted ] = binslin(angle(protraction),rawResponses(protraction)*1000,'equalE',numel(bounds),bounds(1),bounds(end));
        [rSorted ] = binslin(angle(retraction),rawResponses(retraction)*1000,'equalE',numel(bounds),bounds(1),bounds(end));
        
        pResponse = cellfun(@nanmean,pSorted);
        pSEM = cellfun(@(x) nanstd(x)./numel(x),pSorted);
        rResponse = cellfun(@nanmean,rSorted);
        rSEM = cellfun(@(x) nanstd(x)./numel(x),rSorted);
        x = linspace(-99,99,100); 
        sampledBounds = x(~isnan(pSEM) + ~isnan(rSEM) > 0 );
        
        figure(123);subplot(6,10,k)
%         shadedErrorBar(linspace(-99,99,100),smooth(pResponse),smooth(pSEM),'b');
%         hold on; shadedErrorBar(linspace(-99,99,100),smooth(rResponse),smooth(rSEM),'r');
        shadedErrorBar(linspace(-99,99,100),(pResponse),(pSEM),'b');
        hold on; shadedErrorBar(linspace(-99,99,100),(rResponse),(rSEM),'r');
        set(gca,'xlim',[min(sampledBounds) max(sampledBounds)])
        
         figure(124);subplot(6,10,k)
%         shadedErrorBar(linspace(-99,99,100),smooth(pResponse),smooth(pSEM),'b');
%         hold on; shadedErrorBar(linspace(-99,99,100),smooth(rResponse),smooth(rSEM),'r');
        plot(linspace(-99,99,100),normalize_var((pResponse),0,1),'b');
        hold on; plot(linspace(-99,99,100),normalize_var(rResponse,0,1),'r');
        set(gca,'xlim',[min(sampledBounds) max(sampledBounds)])
        
    end
    
end
