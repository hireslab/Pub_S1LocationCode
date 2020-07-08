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

%% Quantifying whisking tuning
whisking = whisking_general(U,'off');
hilbertWhisking = whisking_hilbert(U,popV,'off');

%% decoding whisker position
featureName = {'phase','theta'};
numBins =12 ; %only for theta;

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 3;


%theta set up
theta99 = cellfun(@(x) prctile(x,99),hilbertWhisking.S_ctk.raw,'uniformoutput',0);
thetaRevamped = cellfun(@(x,y) x(:,1)>y(1),hilbertWhisking.S_ctk.raw,theta99,'uniformoutput',0);
theta = hilbertWhisking.S_ctk.raw;
for b = 1:length(thetaRevamped)
    theta{b}(thetaRevamped{b}) = theta99{b}(1);
end
lims = cellfun(@(x) range(x(:,1)), theta);

theta = cellfun(@(x) normalize_var(x(:,1),-1,1),theta,'uniformoutput',0);
thetaResponse = cellfun(@(x,y) binslin(x,y(:),'equalE',numBins+1,-1,1),theta,hilbertWhisking.R_ntk.raw,'uniformoutput',0);

% figure;histogram(cell2mat(cellfun(@(x)
% sum(~isnan(x)),hilbertWhisking.R_ntk.phase,'uniformoutput',0)))
% %distribution of num samples in each bin
%using the median number of samples in each bin to dictate resampling size
% resampNum = round(nanmedian(cell2mat(cellfun(@(x) sum(~isnan(x)),hilbertWhisking.R_ntk.phase,'uniformoutput',0))));

 resampNum = 4000;
    
for k = 1
    
    boostedDmatX = cell(1,length(U));
    boostedDmatY = cell(1,length(U));
    %resampling spikes so all bins have equal sample sizes
    for i = 1:length(U)
        if strcmp(featureName{k},'phase')
            cellfeature = hilbertWhisking.R_ntk.phase{i};
        elseif strcmp(featureName{k},'theta')
            cellfeature = cell2nanmat(thetaResponse{i});
        end
        
        resampedSpikes = nan(resampNum,size(cellfeature,2));
        for b = 1:size(cellfeature,2)
            [~,bootsam] = bootstrp(50,@mean,cellfeature(:,b));
            binresampedSpikes = cellfeature(bootsam,b);
            binresampedSpikes = binresampedSpikes(~isnan(binresampedSpikes));
            
            resampedSpikes(:,b) = binresampedSpikes(1:resampNum);
        end
        
        boostedDmatX{i} = resampedSpikes;
        boostedDmatY = repmat(1:size(cellfeature,2),resampNum,1);
    end
    
    numCellsToSample = [10 25 50 100];
    
    %circularly permuting tuning so tuning curves are equally distributed
    %across all sampled positions
    selCells = datasample(1:length(U),max(numCellsToSample));
    resampCells = boostedDmatX(selCells);
    shiftTuningDmatX = cellfun(@(x) x(:,circshift(1:size(x,2),randperm(size(x,2)))),resampCells,'uniformoutput',0);
    variationDmatX = cellfun(@(x) x(randperm(size(x,1)),:),shiftTuningDmatX,'uniformoutput',0);
    
 
    
    sampleIterations = 10;  %number of iterations of used cells 
    true = cell(length(numCellsToSample),1); 
    predicted = cell(length(numCellsToSample),1);     
    for g = 1:length(numCellsToSample)
        
        for u=1:sampleIterations
        selCells = datasample(1:length(variationDmatX),numCellsToSample(g));
        selCellsDmatX = variationDmatX(selCells);
        
        
        DmatX = cell2mat(cellfun(@(x) x(:),selCellsDmatX,'uniformoutput',0));
        DmatX = (DmatX - nanmean(DmatX)) ./ nanstd(DmatX); %standardization
        
        DmatX = DmatX(:,sum(isnan(DmatX))==0); 
%         DmatX = DmatX(reshape(randperm(numel(DmatX)),size(DmatX))); %shuffling Dmat X 
        
        DmatY = boostedDmatY(:);
        
        mdl.io.X = DmatX;
        mdl.io.Y.normal = DmatY;
        
        mdl = multinomialModel(mdl,DmatX,DmatY,glmnetOpt);
        
        if strcmp(featureName{k},'phase')
            binResolution = 360/numel(unique(DmatY));
        elseif strcmp(featureName{k},'theta')
            binResolution = mean(lims ./ numel(unique(DmatY)));
        end

        true{g}(:,u) = mdl.io.trueXpredicted(:,1);
        predicted{g}(:,u) = mdl.io.trueXpredicted(:,2); 

        end

    end
    
    decodingResolutionMean = cellfun(@(x,y) mean(abs(x(:)-y(:))) * binResolution,true,predicted);
    decodingResolutionSEM  = cellfun(@(x,y) std(abs(x(:)-y(:))) ./ sqrt(sampleIterations), true,predicted);
  
    figure(23);clf
    shadedErrorBar(numCellsToSample,smooth(decodingResolutionMean),decodingResolutionSEM,'k')
    set(gca,'ylim',[0 160],'xlim',[0 100])
    xlabel('number of cells')
    ylabel('decoding resolution (degrees of phase)')
    
    
end
