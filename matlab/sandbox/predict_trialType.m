%Load whisking and neural time series struct
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
% U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
tStruct = object_location_quantification(U,selectedCells,'pole');
touchCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U));
nontouchCells = find(~cellfun(@(x) isfield(x.meta,'responseWindow'),U));

%%
selectedUnits = 1:length(U);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 15;

glmModel = cell(1,numel(selectedUnits));
for i = 1:length(selectedUnits)
    cellNum = selectedUnits(i);
    mask = maskBuilder(U{cellNum});
    spikes = squeeze(U{cellNum}.R_ntk);
    DmatY = U{cellNum}.meta.trialType';
    
    nontouchmask = double(~(mask.touch==1));
    nontouchmask(nontouchmask==0) = nan;
    
    touch = nansum(mask.avail .* mask.touch .* spikes)'; %touches in pole avail period
    nontouch = nansum(mask.avail .* nontouchmask .* spikes)'; %nontouches in pole avail period
    all = nansum(mask.avail .* spikes)'; %pole avail
    
    TypeNames = {'all','touch','nontouch'};
    Types = {all,touch, nontouch};
    
    for g = 1:numel(Types)
        DmatX = Types{g};
        [unique_values,num] = unique(DmatX);
        
        if (sum(DmatX~=0)/numel(DmatX))>.2
            predProb = cell(1,glmnetOpt.numIterations);
            true = cell(1,glmnetOpt.numIterations);
            for p = 1:glmnetOpt.numIterations
                numTrials = size(DmatY,1);
                exampleIdx = 1:numTrials;
                shuffIdx = exampleIdx(randperm(length(exampleIdx)));
                trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
                testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);
                
                trainDmatX = DmatX(trainIdxStartRaw,:);
                trainDmatY = DmatY(trainIdxStartRaw,:);
                testDmatX = DmatX(testIdxStartRaw,:);
                testDmatY = DmatY(testIdxStartRaw,:);
                
                [~,num] = unique(trainDmatX);
                
                cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
                
                fitLambda = cv.lambda_1se;
                iLambda = find(cv.lambda == fitLambda);
                fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
                predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
                true{p} = testDmatY;
                predProb{p} = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
                
            end
            
            glmModel{i}.cellNumber = cellNum;
            glmModel{i}.test.true.(TypeNames{g}) = true;
            glmModel{i}.test.predProb.(TypeNames{g}) = predProb;
        else
            disp(['skipping cell ' num2str(i) ' b/c too few unique responses'])
            glmModel{i}.cellNumber = cellNum;
            glmModel{i}.test.true.(TypeNames{g}) = {DmatY};
            glmModel{i}.test.predProb.(TypeNames{g}) = {ones(length(DmatY),1).*.5};
        end
        glmModel{i}.test.typeNames = TypeNames;
    end
    
end

%% AUC
[tpr{1},far{1}] = cellfun(@(x) roc(cell2mat(x.test.true.all')',cell2mat(x.test.predProb.all')'),glmModel,'uniformoutput',0);
[tpr{2},far{2}] = cellfun(@(x) roc(cell2mat(x.test.true.touch')',cell2mat(x.test.predProb.touch')'),glmModel,'uniformoutput',0);
[tpr{3},far{3}] = cellfun(@(x) roc(cell2mat(x.test.true.nontouch')',cell2mat(x.test.predProb.nontouch')'),glmModel,'uniformoutput',0);


auc = cellfun(@(x,y) 1-cellfun(@(a,b) trapz(a,b),x,y),tpr,far,'uniformoutput',0);

figure(8);clf
for b = 1:length(auc)
    subplot(1,3,b)
    for g = 1:length(tpr{b})
        if auc{b}(g)>.5
            hold on; plot(far{b}{g},tpr{b}{g},'k')
        else
            hold on; plot(far{b}{g},tpr{b}{g},'color',[.8 .8 .8])
        end
    end
    title(['predictive units = ' num2str(sum(auc{b}>.5)) '/' num2str(numel(selectedUnits)) '. mean auc pred = ' num2str(mean(auc{b}(auc{b}>.5)))])
    axis square
    xlabel('false positive rate');
    ylabel('true positive rate');
    set(gca,'xtick',0:.5:1,'ytick',0:.5:1)
end

%%
aucAll = auc

figure(48);clf
scatter(aucOutside(auc>.5),aucInside(auc>.5),'filled','k')
hold on; plot([0 1],[0 1],'-.k')
set(gca,'xlim',[.4 1],'ylim',[.4 1],'ytick',0:.25:1,'xtick',0:.25:1)
axis square

