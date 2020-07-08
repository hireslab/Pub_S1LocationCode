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

%% phase at touch tuning
preDecisionTouches = preDecisionTouchMat(U);
phaseBins = 12; %would recommend even numbers;
for rec = 1:length(U)
 [tVar] = atTouch_sorter(U{rec},touchWindow,preDecisionTouches{rec});
 
 phase = tVar.allTouches.S_ctk(:,5);
 if isfield(U{rec}.meta,'responseWindow')
%      tResponse = U{rec}.meta.responseWindow+find(touchWindow==0) ;
     tResponse = [-25 -5] + find(touchWindow==0); %control to find off touch responses
 else
     tResponse = [5 35] + find(touchWindow==0);
 end
  response = mean(tVar.allTouches.R_ntk(:,tResponse(1):tResponse(2)),2);
 
 [sorted,sby] = binslin(phase,response,'equalE',phaseBins+1,-pi,pi);
 hilbertTouch.R_ntk.phase{rec} = cell2nanmat(sorted);
 hilbertTouch.S_ctk.phase = linspace(-pi,pi,phaseBins);
end

%% decoding touch position by phase
% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 10;

%using the median number of samples in each bin to dictate resampling size
resampNumTMP = cell2mat(cellfun(@(x) sum(~isnan(x)),hilbertTouch.R_ntk.phase(selectedCells),'uniformoutput',0)');
%visualing sampling distribution across all cells
figure(9);clf;shadedErrorBar(linspace(-pi,pi,12),mean(resampNumTMP),std(resampNumTMP))
set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'},'xlim',[-pi pi])
ylabel('mean num samples')

% resampNum = round(mean(resampNumTMP(:)));
resampeNum = 50;

boostedDmatX = cell(1,length(selectedCells)); 

%resampling all touch bins to have the same number of "touches"
for i = 1:length(selectedCells)
    cellfeature = hilbertTouch.R_ntk.phase{selectedCells(i)};
    
    if sum(sum(~isnan(cellfeature))==0)>0
        disp(['skipping cell ' num2str(selectedCells(i)) ' because unsampled bin'])
        resampedSpikes = nan(resampNum,size(cellfeature,2));
    else
        resampedSpikes = nan(resampNum,size(cellfeature,2));
        for b = 1:size(cellfeature,2)
            [~,bootsam] = bootstrp(round(resampNum*1.5),@mean,cellfeature(:,b));
            binresampedSpikes = cellfeature(bootsam,b);
            binresampedSpikes = binresampedSpikes(~isnan(binresampedSpikes));
            
            resampedSpikes(:,b) = binresampedSpikes(1:resampNum);
        end
    end
    boostedDmatX{i} = resampedSpikes;
    boostedDmatY = repmat(1:size(cellfeature,2),resampNum,1);
end

%circularly permuting data for varied tunings and shuffling data for
%variability

numCellsToSample = [5 10 15 20 30 40 50 75 100 250 500 750 1000];

selCells = datasample(1:length(selectedCells),max(numCellsToSample));
resampCells = boostedDmatX(selCells);
shiftTuningDmatX = cellfun(@(x) x(:,circshift(1:size(x,2),randperm(size(x,2)))),resampCells,'uniformoutput',0);
shuffledResponseDmatX = cellfun(@(x) x(randperm(size(x,1)),:),shiftTuningDmatX,'uniformoutput',0);

%multinomial model 
decodingResolutionMean = zeros(1,length(numCellsToSample)); 
decodingResolutionSEM = zeros(1,length(numCellsToSample)); 

for g = 1:length(numCellsToSample)
    
    selCells = datasample(1:length(shuffledResponseDmatX),numCellsToSample(g));
    selCellsDmatX = shuffledResponseDmatX(selCells);

    selCellsDmatY = boostedDmatY;
    
    DmatX = cell2mat(cellfun(@(x) x(:),selCellsDmatX,'uniformoutput',0));
    DmatX = (DmatX - mean(DmatX)) ./ std(DmatX); %standardization
    DmatY = selCellsDmatY(:);
    
    mdl.io.X = DmatX;
    mdl.io.Y.normal = DmatY;
    
     mdl = multinomialModel(mdl,DmatX,DmatY,glmnetOpt);
%     decoderPerformance(mdl)
%     
    binResolution = 360/numel(unique(DmatY));
    decodingResolutionMean(g) = mean(abs(mdl.io.trueXpredicted(:,1) - mdl.io.trueXpredicted(:,2))) * binResolution;
    decodingResolutionSEM(g) = (std(abs(mdl.io.trueXpredicted(:,1) - mdl.io.trueXpredicted(:,2))) ./ glmnetOpt.numIterations) * binResolution;
    
end

figure(23);clf
shadedErrorBar(numCellsToSample,decodingResolutionMean,decodingResolutionSEM,'k')
set(gca,'ylim',[0 160],'xlim',[0 max(numCellsToSample)])
xlabel('number of cells')
ylabel('decoding resolution (degrees of phase)')








    

