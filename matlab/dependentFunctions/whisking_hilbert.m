function hilbertWhisking = whisking_hilbert(array,popV,displayOpt)

if (nargin < 2), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));
%%
% numBins = 20; 
pThreshold = 0.01; 

% whiskingCells = find(whisking.matrix(1,:)==-1);

for i = 1:length(array)
    currArray = array{i};
    
    
    [objmask]= maskBuilder(currArray);
    angle = squeeze(currArray.S_ctk(1,:,:));
    amplitude = squeeze(currArray.S_ctk(3,:,:));
    midpoint = squeeze(currArray.S_ctk(4,:,:));
    phase = squeeze(currArray.S_ctk(5,:,:));
    spikes = squeeze(currArray.R_ntk(:,:,:));
    
    
    twFilter = objmask.whisking .* objmask.touch;
    
    filt_amplitude = twFilter .* amplitude;
    filt_midpoint = twFilter .* midpoint;
    filt_phase = twFilter .* phase;
    filt_angle = twFilter .* angle;
    filt_spikes = twFilter .* spikes;
    
    inputX = [filt_amplitude(:) filt_midpoint(:) filt_phase(:) filt_angle(:)];

%     for k = 1:size(inputX,2)
%         [response{k},stimulus{k}] = binslin(inputX(:,k),filt_spikes(:),'equalE',numBins,min(inputX(:,k)), max(inputX(:,k)));
%     end
%     
    bigBounds = popV{1}.allTouches.theta.bounds;
    littleBounds = popV{1}.allTouches.phase.bounds;
    response = cell(1,size(inputX,2)); 
    stimulus = cell(1,size(inputX,2)); 
    for k = 1:size(inputX,2)
        if k ~= 3
        [response{k},stimulus{k}] = binslin(inputX(:,k),filt_spikes(:),'equalE',numel(bigBounds),bigBounds(1),bigBounds(end));
        else
        [response{k},stimulus{k}] = binslin(inputX(:,k),filt_spikes(:),'equalE',numel(littleBounds)+1,littleBounds(1),littleBounds(end));    
        end
    end
   
    average_response = cellfun(@(x) cellfun(@nanmean,x), response,'uniformoutput',0);
    average_stimulus = cellfun(@(x) cellfun(@nanmean,x), stimulus,'uniformoutput',0);
    SEM_response = cellfun(@(x) cellfun(@(y) nanstd(y)./sqrt(numel(y)),x), response,'uniformoutput',0);
    
    %find to see if bins are sig different in each whisking feature 
    pcompare = cell(length(response),1) ; 
    raw_response = cell(1,length(average_response)); 
    for j = 1:length(response) 
        compare_matrix = cell2nanmat(response{j});
        p = anova1(compare_matrix,[],'off');
        
        if p<pThreshold     
            [~,maxIdx] = max(average_response{j});
            peakResponse = compare_matrix(:,maxIdx);
            peakResponse = peakResponse(~isnan(peakResponse)); 
            
            for g = 1:size(compare_matrix,2)
                testResponse = compare_matrix(:,g);
                testResponse = testResponse(~isnan(testResponse)); 
                 
                [~,pcompare{j}(g)] = ttest2(peakResponse,testResponse);
            end
        end
        raw_response{j} = compare_matrix; 
    end

    if willdisplay
        figure(8);clf
        featureNames = {'amp','midpoint','phase','angle'};
        for j = 1:size(average_response,2)
            subplot(2,2,j)
            shadedErrorBar(average_stimulus{j},average_response{j}*1000,SEM_response{j}*1000)
            set(gca,'xlim',[min(average_stimulus{j}) max(average_stimulus{j})])
            title(featureNames{j})
            
            if ~isempty(pcompare{j})
                peakResponse = find(pcompare{j}==1); 
                sigDiffResponses = find(pcompare{j}<pThreshold==1); 
                [~,idx] = min(abs(peakResponse-sigDiffResponses));
                firstDiff = sigDiffResponses(idx);
                
                hold on; scatter(average_stimulus{j}(peakResponse),average_response{j}(peakResponse)*1000,'filled','g')
                scatter(average_stimulus{j}(firstDiff),average_response{j}(firstDiff)*1000,'filled','r')
            end
        end
    end
    
    hilbertWhisking.cellVarNames = {'theta','amp','midpoint','phase'};
    hilbertWhisking.S_ctk.raw{i} = inputX;
    hilbertWhisking.R_ntk.raw{i} = filt_spikes; 
    
    hilbertWhisking.S_ctk.theta = mean([bigBounds(1:end-1);bigBounds(2:end)]); 
    hilbertWhisking.S_ctk.amp = mean([bigBounds(1:end-1);bigBounds(2:end)]);
    hilbertWhisking.S_ctk.midpoint = mean([bigBounds(1:end-1);bigBounds(2:end)]); 
    hilbertWhisking.S_ctk.phase = littleBounds; %mean([littleBounds(1:end-1);littleBounds(2:end)]); 
   
    hilbertWhisking.R_ntk.theta{i} = raw_response{4}; 
    hilbertWhisking.R_ntk.amp{i}  = raw_response{1}; 
    hilbertWhisking.R_ntk.midpoint{i} = raw_response{2};
    hilbertWhisking.R_ntk.phase{i} = raw_response{3};
    
    
end

%% GLMNET ATTEMPT BUT FRs may be too low. If wanting to do consider
%%binning and use poisson family option
%     glmnetOpt = glmnetSet;
% glmnetOpt.standardize = 0; %set to 0 b/c already standardized
% glmnetOpt.alpha = 0.95;
% glmnetOpt.xfoldCV = 5;
% glmnetOpt.numIterations = 10;

%     dmatXY = [filt_amplitude(:) filt_midpoint(:) filt_phase(:) filt_angle(:) filt_spikes(:)];
%     dmatXY = dmatXY(~sum(isnan(dmatXY),2)>0 ,:) ;
%
%     numTP = size(dmatXY,1);
%     numFeatures = size(dmatXY,2)-1;
%
%     %hankel matrix construction
%     lag = -20:2:20;
%     lagMat = nan(length(lag),numTP+1) ;
%     for k = 1:length(lag)
%         lagMat(k,:) = lag(k):numTP+lag(k);
%     end
%
%     %toss 0, negative, and exceeding indices
%     lagMat(:,sum(lagMat<=0)>0)=[];
%     lagMat(:,sum(lagMat>numTP)>0)=[];
%
%
%     midIdx = lagMat(find(lag==0),:);
%
%     idxDmat = dmatXY(lagMat,1:numFeatures);
%     DmatX = reshape(idxDmat,size(lagMat,2),length(lag) * numFeatures);
%     DmatXNorm = (DmatX-mean(DmatX)) ./ std(DmatX);
%     DmatY = dmatXY(midIdx,end);
%
%
%
%     cv = cvglmnet(DmatXNorm,DmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
%     cvglmnetPlot(cv)
%
%     fitLambda = cv.lambda_1se;
%     iLambda = find(cv.lambda == fitLambda);
%     fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
%     predicts = cvglmnetPredict(cv,DmatXNorm,fitLambda); % this outputs testDmatX * fitCoeffs
%     predProb = sigmoid(predicts);

