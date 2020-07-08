%% Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% neurons and their ability to decode choice 
%using average fr from pOnset to first lick or using the spike train from pole onset to median lick time. 

glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.05;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 5;
mcc = cell(1,2); 

for i = 1:length(U)
    
    if strcmp(U{i}.meta.layer,'BVL5b')        
        disp(['iterating for cell ' num2str(i) '/' num2str(length(U))])
        
        masks = maskBuilder(U{i});
        spks = squeeze(U{i}.R_ntk);
        
        meanPreDFR = (nanmean(spks.* masks.onsettolick) * 1000)';
        spikeTrainPreDecision = (spks.* masks.onsettoMedianlick)';
        selectedTPs = sum((~isnan(spikeTrainPreDecision)),1) == U{i}.k; 
        spikeTrainPreDecision = spikeTrainPreDecision(:,selectedTPs);
        
        DmatX = {meanPreDFR,spikeTrainPreDecision}; 
        DmatY = ~(masks.nonlickTrials)';
        
        if sum(DmatY==0) < 5
            disp(['skipping cell ' num2str(i) ' because not enough non-lick trials'])
        else
            for g = 1:length(DmatX)
                for d = 1:glmnetOpt.numIterations
                    [trainIdx,testIdx] = stratifiedIndices(DmatY,.7);
                    
                    trainDmatX = DmatX{g}(cell2mat(trainIdx'),:);
                    testDmatX = DmatX{g}(cell2mat(testIdx'),:);
                    trainDmatY = DmatY(cell2mat(trainIdx'),:);
                    testDmatY = DmatY(cell2mat(testIdx'),:);
                    
                    cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
                    
                    fitLambda = cv.lambda_1se;
                    iLambda = find(cv.lambda == fitLambda);
                    fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
                    predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
                    predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
                    
                    mcc{g}{i}(d) = mccCalculator(testDmatY,predProb>.5);
                end
            end
        end
    end
    
end

meanMCC = cellfun(@(x) cellfun(@mean,x) , mcc,'uniformoutput',0);

choicePred.rowNames = {'mean FR pre decision','spike train pre decision','mouse accuracy'}; 
choicePred.mat = [cell2mat(meanMCC') ; cellfun(@(x) mean(x.meta.trialCorrect),U)];

%% Across what timescales is object locatoin decoding best decoded from? 
%using the mean fr or spike train from pole onset to average pole down to decode pole location

glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.05;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 5;
RsqdTrialResponses = cell(1,2); 
RsqdTouchResponses = cell(1,2); 

for i = 1:length(U)
    disp(['iterating for cell ' num2str(i) '/' num2str(length(U))])
    
    masks = maskBuilder(U{i});
    spks = squeeze(U{i}.R_ntk);
    
    meanFR = (nanmean(spks.* masks.availtoMeanOffset) * 1000)';
    spikeTrain = (spks.* masks.availtoMeanOffset)';
    selectedTPs = sum((~isnan(spikeTrain)),1) == U{i}.k;
    spikeTrain = spikeTrain(:,selectedTPs);
    
    DmatX = {spikeTrain,meanFR};
    DmatY = U{i}.meta.motorPosition';
%     
%     figure(23);subplot(5,11,i)
%     scatter(DmatX{1},DmatY,'.k'); 
    
    for g = 1:length(DmatX)
        for d = 1:glmnetOpt.numIterations
            [trainIdx,testIdx] = randomizedIndices(DmatY,.7);
            
            trainDmatX = DmatX{g}(trainIdx',:);
            testDmatX = DmatX{g}(testIdx',:);
            trainDmatY = DmatY(trainIdx',:);
            testDmatY = DmatY(testIdx',:);
            
            cv = cvglmnet(trainDmatX,trainDmatY,'gaussian',glmnetOpt,[],glmnetOpt.xfoldCV);
            
            fitLambda = cv.lambda_1se;
            iLambda = find(cv.lambda == fitLambda);
            fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
            predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
            
            SSR = sum((predicts-mean(DmatY)).^2);
            SSE = sum((testDmatY-predicts).^2);
            SSTO = SSR + SSE;
            RsqdTrialResponses{g}{i}(d) = SSR./SSTO;
        end
    end
    
end
meanTrialRsqd = cellfun(@(x) cellfun(@mean,x) , RsqdTrialResponses,'uniformoutput',0);


%using the mean fr of all responses per trial OR each touch individually 
viewWindow = -25:50; 
U = defTouchResponse(U,.99,'off');
popPredicts = cell(1,length(U));
preDecisionTouches = preDecisionTouchMat(U);

for i = 1:length(U) 
    if isfield(U{i}.meta,'responseWindow')
        disp(['iterating for cell ' num2str(i) '/' num2str(length(U))])
        
        tV = atTouch_sorter(U{i},viewWindow,preDecisionTouches{i});
        
        spikeTrainResponse = tV.allTouches.spikeMat(:,find(viewWindow==0)+U{i}.meta.responseWindow(1) : find(viewWindow==0)+U{i}.meta.responseWindow(2)); %capture responses within tuch response window
%         spikeTrainResponse = tV.allTouches.spikeMat(:,find(viewWindow==0)+5: find(viewWindow==0)+35); %use all responses from 5:35ms post touch
        meanResponse = nanmean(spikeTrainResponse,2);
        motors = tV.allTouches.variables(:,7);
        Umotors = unique(motors);

        allTouchMean = nan(length(Umotors),1); 
        for f = 1:length(Umotors)
            allTouchMean(f) = mean(meanResponse(find(motors == Umotors(f))));
        end
 
        DmatX = {spikeTrainResponse,meanResponse,allTouchMean};
        DmatY = {motors,motors,Umotors}; %angle at touch
        
        for g = 1:length(DmatX)
            for d = 1:glmnetOpt.numIterations
                [trainIdx,testIdx] = randomizedIndices(DmatY{g},.7);
                
                trainDmatX = DmatX{g}(trainIdx',:);
                testDmatX = DmatX{g}(testIdx',:);
                trainDmatY = DmatY{g}(trainIdx',:);
                testDmatY = DmatY{g}(testIdx',:);
                
                cv = cvglmnet(trainDmatX,trainDmatY,'gaussian',glmnetOpt,[],glmnetOpt.xfoldCV);
                
                fitLambda = cv.lambda_1se;
                iLambda = find(cv.lambda == fitLambda);
                fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
                predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
                
               tXp = [ testDmatY predicts];
                
                SSR = sum((predicts-mean(DmatY{g})).^2);
                SSE = sum((testDmatY-predicts).^2);
                SSTO = SSR + SSE;
                RsqdTouchResponses{g}{i}(d) = SSR./SSTO;
            end
             popPredicts{i}{g} = tXp;
        end
    end
end


meanTouchRsqd = cellfun(@(x) cellfun(@mean,x) , RsqdTouchResponses,'uniformoutput',0);

positionDecoding.rowNames = {'spikeT avail','FRm avail','spikeT touch','FRm touch','FRm touches/trial'};
positionDecoding.matrix = cell2mat([meanTrialRsqd meanTouchRsqd]');

figure(20);clf;scatter(repmat(1:5,1,numel(U)),positionDecoding.matrix(:),'filled','markerfacecolor',[.7 .7 .7])
hold on; errorbar(1:5,nanmean(positionDecoding.matrix,2),nanstd(positionDecoding.matrix,[],2)./nansum(~isnan(positionDecoding.matrix),2),'ko')
set(gca,'xtick',1:5,'xticklabel',positionDecoding.rowNames)
    
    

