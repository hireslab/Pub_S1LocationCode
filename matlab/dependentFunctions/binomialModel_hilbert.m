function mdl = binomialModel_hilbert(DmatX,DmatY,selectedArray,glmnetOpt,mdl)


%GLMNET
sangles = cell(1,glmnetOpt.numIterations);
sHeat = cell(1,glmnetOpt.numIterations);
sFitCoeffs = cell(1,glmnetOpt.numIterations);
fitDevExplained = nan(1,glmnetOpt.numIterations);
fitDevExplainedFull = nan(1,glmnetOpt.numIterations);
mccDiff = nan(glmnetOpt.numIterations,2);
numTrials = size(DmatY,1)/length(glmnetOpt.buildIndices);

for p = 1:glmnetOpt.numIterations
    disp(['running glmnet with elasticNet ' num2str(glmnetOpt.alpha) ' for iteration ' num2str(p) '/' num2str(glmnetOpt.numIterations)])
    % test and train sets
    % might want to include stratification to make sure all ranges
    % are used
    exampleIdx = 1:numTrials;
    shuffIdx = exampleIdx(randperm(length(exampleIdx)));
    trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
    testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);
    
    trainIdxStart = ((trainIdxStartRaw-1).*length(glmnetOpt.buildIndices))+1;
    testIdxStart = ((testIdxStartRaw-1).*length(glmnetOpt.buildIndices))+1;
    
    trainIdx = trainIdxStart'+(0:length(glmnetOpt.buildIndices)-1);
    testIdx = testIdxStart'+(0:length(glmnetOpt.buildIndices)-1);
    
    trainDmatX = DmatX(trainIdx',:);
    trainDmatY = DmatY(trainIdx',:);
    testDmatX = DmatX(testIdx',:);
    testDmatY = DmatY(testIdx',:);
    
    %Check that design matrix is properly indexed when we split to test
    %and train set
    %                 figure(19);clf
    %                 numTrialsToPlot = 3;
    %                 startTouch = datasample(1:length(trainIdxStart)-numTrialsToPlot,1);
    %                 imagesc(trainDmatX(length(glmnetOpt.buildIndices)*startTouch+1:(length(glmnetOpt.buildIndices)*startTouch+1)+length(glmnetOpt.buildIndices).*numTrialsToPlot,:))
    %                 set(gca,'xtick',[],'ytick',[])
    %
    cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
%     cvglmnetPlot(cv)
    
    fitLambda = cv.lambda_1se;
    iLambda = find(cv.lambda == fitLambda);
    fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
    predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
    predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
    
    %This is to check that our deviance calculation is correct. Trained
    %model outputs deviance explained but we need for test model which
    %we show below.
    %         trainnullLL = sum(log(binopdf(trainDmatY,1,(ones(length(trainDmatY),1)*mean(trainDmatY)))));
    %         trainpredProb = sigmoid([ones(length(trainDmatX),1) trainDmatX] *fitCoeffs);
    %         trainfullLL = sum(log(binopdf(trainDmatY,1, trainpredProb)));
    %         trainDevianceExplained = 1-(-trainfullLL)/(-trainnullLL);
    %         devExplained = cv.glmnet_fit.dev(iLambda);
    
    % goodness of fit using deviance explained in selected touch
    % response windows
    mu = mean(testDmatY);
    reshapedDmatY = reshape(testDmatY,length(glmnetOpt.buildIndices),length(testDmatY)./length(glmnetOpt.buildIndices));
    reshapedPredProb = reshape(predProb,length(glmnetOpt.buildIndices),length(predProb)./length(glmnetOpt.buildIndices));
    selIdx = [find(glmnetOpt.buildIndices==selectedArray.meta.responseWindow(1)):find(glmnetOpt.buildIndices==selectedArray.meta.responseWindow(2))];
    selectedDmatY = reshapedDmatY(selIdx,:);
    selectedPredProb = reshapedPredProb(selIdx,:);
    
    fullLL = sum(log(binopdf(selectedDmatY(:),1,selectedPredProb(:))));
    nullLL = sum(log(binopdf(selectedDmatY(:),1,(ones(length(selectedPredProb(:)),1)*mu))));
    saturatedLL = sum(log(binopdf(selectedDmatY(:) ,1,selectedDmatY(:) )));
    
    fitDevExplained(p) = 1 - (saturatedLL-fullLL)/(saturatedLL-nullLL);
    
    %goodness of fit using all indices
    fullLL = sum(log(binopdf(testDmatY,1,predProb)));
    nullLL = sum(log(binopdf(testDmatY,1,(ones(length(predProb),1)*mu))));
    saturatedLL = sum(log(binopdf(testDmatY,1,testDmatY)));
    fitDevExplainedFull(p) = 1 - (saturatedLL-fullLL)/(saturatedLL-nullLL)';
    
    % goodness of fit using MCC
    mu = mean(testDmatY);
    pcBins = 0:.05:1;
    spksProbNull = (ones(length(predProb),1)*mu) > pcBins;
    spksProbFull = predProb>pcBins;
    for j = 1:length(0:.05:1)
        mccFull(j) = mccCalculator(testDmatY,double(spksProbFull(:,j)));
        mccNull(j) = mccCalculator(testDmatY,double(spksProbNull(:,j)));
    end
    [mccDiff(p,1), maxidx] = max(mccFull-mccNull);
    mccDiff(p,2) = pcBins(maxidx);
    
    %variables for recreating heat map
    rawFR = reshape(testDmatY,length(glmnetOpt.buildIndices),length(testIdxStart));
    predFR = reshape(predProb,length(glmnetOpt.buildIndices),length(testIdxStart));
    bounds = [-100:2:100];
    rawAngle = mdl.raw.angle;
    [sangles{p},sidx] = sort(rawAngle(testIdxStartRaw));
    sHeat{p} = predFR(:,sidx)';
    sRaw{p} = rawFR(:,sidx)';
    sFitCoeffs{p} = fitCoeffs;
    
end


%constructing output structure of model
mdl.modelParams = glmnetOpt;

touchShiftIdx = mdl.basisFunctions.touch;

coeffsIter = cell2mat(sFitCoeffs);
mdl.coeffs.raw = cell2mat(sFitCoeffs);


featureLagsEnds = [cumsum(mdl.io.selectedFeatures.dims)+1];
featureLagsStarts = [2 cumsum(mdl.io.selectedFeatures.dims)+2];

for j = 1:length(featureLagsEnds)
    mdl.coeffs.(mdl.io.selectedFeatures.name{j}) = mean(coeffsIter(featureLagsStarts(j):featureLagsEnds(j),:),2);
end

mdl.raw.spikes = reshape(DmatY,length(glmnetOpt.buildIndices),length(DmatY)/length(glmnetOpt.buildIndices))';

mdl.predicted.spikeTestRaw = cell2mat(sRaw');
mdl.predicted.spikeProb = cell2mat(sHeat');
mdl.predicted.angles = cell2mat(sangles');

mdl.gof.devExplained = fitDevExplained;
mdl.gof.mcc = mccDiff;


%predictedHeat
[sorted, sortedBy ,~]=binslin(cell2mat(sangles'),cell2mat(sHeat'),'equalE',numel(bounds)+1,bounds(1),bounds(end));
[sortedRaw, sortedByRaw ,~]=binslin(cell2mat(sangles'),cell2mat(sRaw'),'equalE',numel(bounds)+1,bounds(1),bounds(end));
selBins = ~cellfun(@isempty ,sortedBy);

spksBinned = cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
spksBinnedRaw = cell2mat(cellfun(@(x) mean(x,1),sortedRaw,'uniformoutput',0));

heatPlot = spksBinned(selBins,:);
heatPlotRaw = spksBinnedRaw(selBins,:);

modelBounds = bounds(selBins);

normHeatPlot =(heatPlot - min(heatPlot(:))) ./ (max(heatPlot(:))-min(heatPlot(:)));
normRawPlot = (heatPlotRaw - min(heatPlotRaw(:))) ./ (max(heatPlotRaw(:))-min(heatPlotRaw(:)));

mdl.heat.angle = modelBounds;
mdl.heat.matrix = normHeatPlot;
mdl.heat.matrixRaw = normRawPlot;


