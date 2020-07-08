function mdl = gaussianModel(DmatX,DmatY,selectedArray,glmnetOpt,mdl)

numTrials = size(DmatY,1);
fitDevExplained = nan(1,glmnetOpt.numIterations);

for p = 1:glmnetOpt.numIterations
    disp(['running glmnet with elasticNet ' num2str(glmnetOpt.alpha) ' for iteration ' num2str(p) '/' num2str(glmnetOpt.numIterations)])
    % test and train sets
    % might want to include stratification to make sure all ranges
    % are used
    exampleIdx = 1:numTrials;
    shuffIdx = exampleIdx(randperm(length(exampleIdx)));
    trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
    testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);
    
    trainDmatX = DmatX(trainIdxStartRaw,:);
    trainDmatY = DmatY(trainIdxStartRaw,:);
    testDmatX = DmatX(testIdxStartRaw,:);
    testDmatY = DmatY(testIdxStartRaw,:);
    
    %Check that design matrix is properly indexed when we split to test
    %and train set
    %                 figure(19);clf
    %                 numTrialsToPlot = 3;
    %                 startTouch = datasample(1:length(trainIdxStart)-numTrialsToPlot,1);
    %                 imagesc(trainDmatX(length(glmnetOpt.buildIndices)*startTouch+1:(length(glmnetOpt.buildIndices)*startTouch+1)+length(glmnetOpt.buildIndices).*numTrialsToPlot,:))
    %                 set(gca,'xtick',[],'ytick',[])
    %
    cv = cvglmnet(trainDmatX,trainDmatY,'gaussian',glmnetOpt,[],glmnetOpt.xfoldCV);
    %     cvglmnetPlot(cv)
    
    fitLambda = cv.lambda_1se;
    iLambda = find(cv.lambda == fitLambda);
    fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
    
    %This is to check that our deviance calculation is correct. Trained
    %model outputs deviance explained but we need for test model which
    %we show below.
%     model = exp([ones(length(trainDmatX),1),trainDmatX]*fitCoeffs); %exponential link function 
%     mu = mean(trainDmatY); % null poisson parameter
%     trainnullLL = sum(log(poisspdf(trainDmatY,mu)));
%     saturatedLogLikelihood = sum(log(poisspdf(trainDmatY,trainDmatY)));
%     trainfullLL = sum(log(poisspdf(trainDmatY,model)));
%     trainDevExplained(p) = (trainfullLL - trainnullLL)/(saturatedLogLikelihood - trainnullLL);
%     devExplained(p) = cv.glmnet_fit.dev(iLambda);
    
    
    model = [ones(length(testDmatX),1),testDmatX]*fitCoeffs;
    mu = mean(testDmatY); % null poisson parameter
    nullLogLikelihood = sum(log(normpdf(testDmatY,mu)));
    saturatedLogLikelihood = sum(log(normpdf(testDmatY,testDmatY)));
    fullLogLikelihood = sum(log(normpdf(testDmatY,model)));
    fitDevExplained(p) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    aic(p) = -2.*(fullLogLikelihood) + (2*size(DmatX,2)); 
    devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
    
    
    %   %variables for recreating heat map
    sInput{p} = testDmatX;
    sHeat{p} = model;
    sRaw{p} = testDmatY;
    sPole{p} = mdl.raw.trimmedPole(testIdxStartRaw);
    sFitCoeffs{p} = fitCoeffs;
end


%constructing output structure of model
mdl.modelParams = glmnetOpt;

mdl.coeffs.raw = cell2mat(sFitCoeffs);

mdl.predicted.inputX = sInput;
mdl.predicted.spikeTestRaw = sRaw;
mdl.predicted.spikeProb = sHeat;
mdl.predicted.pole = sPole;

mdl.gof.devExplained = fitDevExplained;
mdl.gof.aic = aic; 
