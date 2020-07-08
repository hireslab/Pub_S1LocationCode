function [glmModel] = designMatrixBlocks_v2(selectedArray,glmnetOpt,glmModel)

preDecisionTouches = preDecisionTouchMat(selectedArray);
viewWindow = -25:50;

basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);

for i = 1:length(selectedArray)
    currentCell = selectedArray{i};
    
    [tVar] = atTouch_sorter(selectedArray{i},viewWindow,preDecisionTouches{i}); %use this to fill velocity pre-touch 
    
    %Defining features
    touchMat = zeros(currentCell.t,currentCell.k); 
    touchIdx = [find(currentCell.S_ctk(9,:,:)==1) ;find(currentCell.S_ctk(12,:,:)==1)];
    touchMat(touchIdx) = 1;
    spikes = squeeze(currentCell.R_ntk);
    phase = squeeze(currentCell.S_ctk(5,:,:)) .* touchMat; %these are at touch features we care about
    amplitude = squeeze(currentCell.S_ctk(3,:,:)) .* touchMat; %these are at touch features we care about
    midpoint = squeeze(currentCell.S_ctk(4,:,:)) .* touchMat; %these are at touch features we care about
    angle = squeeze(currentCell.S_ctk(1,:,:)) .* touchMat; %these are at touch features we care about
    kappa = squeeze(currentCell.S_ctk(17,:,:)) .* touchMat; %kappa at touch
    pt_velocity = touchMat; 
    pt_velocity(pt_velocity==1) =  tVar.allTouches.S_ctk(:,2); 
    
    DKappa = squeeze(currentCell.S_ctk(6,:,:));
    DTheta = squeeze(currentCell.S_ctk(18,:,:));
    
    %Get rid of touches with nan values
    nanidx = unique([find(isnan(midpoint(touchIdx))) ; find(isnan(amplitude(touchIdx))) ; find(isnan(phase(touchIdx)))] );
    touchIdx(nanidx)=[];
    
    %convoluting features with basis fucntions
    angle_conv = conv2(basisFunction,1,angle,'same');
    phase_conv = conv2(basisFunction,1,phase,'same');
    amplitude_conv = conv2(basisFunction,1,amplitude,'same');
    midpoint_conv = conv2(basisFunction,1,midpoint,'same');
    kappa_conv = conv2(basisFunction,1,kappa,'same');
    ptVel_conv = conv2(basisFunction,1,pt_velocity,'same'); 
    DKappa_conv = DKappa;
    DTheta_conv = DTheta; 
%     DKappa_conv = conv2(basisFunction,1,DKappa,'same'); %no convolution for Dkappa

    %INDEX BUILDER FOR LAGS
    startIdx = glmnetOpt.buildIndices'+touchIdx';
    
    %LAG AND LAGS APART
    shiftIdx = nan(length(glmnetOpt.bf.indicesToAdd ),numel(startIdx));
    touchShiftIdxRaw=nan(numel(glmnetOpt.buildIndices),length(glmnetOpt.bf.indicesToAdd ));
    for b = 1:length(glmnetOpt.bf.indicesToAdd )
        tmpIdx{b} = startIdx + glmnetOpt.bf.indicesToAdd (b);
        shiftIdx(b,:) = tmpIdx{b}(:);
        
        touchShiftIdxRaw(:,b) = (glmnetOpt.buildIndices + glmnetOpt.bf.indicesToAdd (b))==0;
        touchShiftIdxRaw(find(touchShiftIdxRaw(:,b)==1)+(-glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth),b) = basisFunction;
    end
    
    touchShiftIdx = touchShiftIdxRaw(:,(glmnetOpt.bf.indicesToAdd <=0));
    
    %Feature blocks for design matrix construction 
    glmModel{i}.io.components.touch = repmat(touchShiftIdx,length(touchIdx),1);
    glmModel{i}.io.components.angle = angle_conv(shiftIdx');
    glmModel{i}.io.components.phase = phase_conv(shiftIdx');
    glmModel{i}.io.components.amplitude = amplitude_conv(shiftIdx');
    glmModel{i}.io.components.midpoint = midpoint_conv(shiftIdx');
    glmModel{i}.io.components.kappa = kappa_conv(shiftIdx'); 
    glmModel{i}.io.components.pt_velocity = ptVel_conv(shiftIdx'); 
    glmModel{i}.io.components.DKappa = DKappa_conv(shiftIdx');
    glmModel{i}.io.components.DTheta = DTheta_conv(shiftIdx'); 
    
    glmModel{i}.basisFunctions.touch = touchShiftIdx;
    glmModel{i}.basisFunctions.features = touchShiftIdxRaw;
    
    DmatY = spikes(startIdx(:));
    rawAngle = angle(touchIdx);
    glmModel{i}.io.DmatY = DmatY;
    glmModel{i}.raw.angle = rawAngle;
    
end