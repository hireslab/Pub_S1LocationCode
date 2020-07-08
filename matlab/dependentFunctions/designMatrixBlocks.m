function [glmModel] = designMatrixBlocks(selectedArray,glmnetOpt,glmModel)

basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);

for i = 1:length(selectedArray)
    currentCell = selectedArray{i};
    
    %Defining features
    touchDurMat = zeros(currentCell.t,currentCell.k); 
    touchIdx = [find(currentCell.S_ctk(9,:,:)==1) ;find(currentCell.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(currentCell.S_ctk(10,:,:)==1) ;find(currentCell.S_ctk(13,:,:)==1)];
    spikes = squeeze(currentCell.R_ntk);
    phase = squeeze(currentCell.S_ctk(5,:,:));
    amplitude = squeeze(currentCell.S_ctk(3,:,:));
    midpoint = squeeze(currentCell.S_ctk(4,:,:));
    angle = squeeze(currentCell.S_ctk(1,:,:));
    velocity = squeeze(currentCell.S_ctk(2,:,:));
    curvature = squeeze(currentCell.S_ctk(17,:,:));
    DKappa = squeeze(currentCell.S_ctk(6,:,:));
    for d = 1:length(touchIdx)
        touchDurMat(touchIdx(d):touchOffIdx(d)) = 1;
    end
    
    %Get rid of touches with nan values
    nanidx = unique([find(isnan(midpoint(touchIdx))) ; find(isnan(amplitude(touchIdx))) ; find(isnan(phase(touchIdx))) ; find(isnan(curvature(touchIdx))) ] );
    touchIdx(nanidx)=[];
    
    %convoluting features with basis fucntions
    angle_conv = conv2(basisFunction,1,angle,'same');
    velocity_conv = conv2(basisFunction,1,velocity,'same'); 
    phase_conv = conv2(basisFunction,1,phase,'same');
    amplitude_conv = conv2(basisFunction,1,amplitude,'same');
    midpoint_conv = conv2(basisFunction,1,midpoint,'same');
    curvature_conv = conv2(basisFunction,1,curvature,'same');
    DKappa_conv = conv2(basisFunction,1,DKappa,'same');
    touchDur_conv = conv2(basisFunction,1,touchDurMat,'same'); 

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
    glmModel{i}.io.components.touchDur = touchDur_conv(shiftIdx(glmnetOpt.bf.indicesToAdd <=0,:)'); 
    glmModel{i}.io.components.angle = angle_conv(shiftIdx');
    glmModel{i}.io.components.phase = phase_conv(shiftIdx');
    glmModel{i}.io.components.amplitude = amplitude_conv(shiftIdx');
    glmModel{i}.io.components.midpoint = midpoint_conv(shiftIdx');
    glmModel{i}.io.components.curvature = curvature_conv(shiftIdx');
    glmModel{i}.io.components.DKappa = DKappa_conv(shiftIdx');
    glmModel{i}.io.components.velocity = velocity_conv(shiftIdx'); 
    
    glmModel{i}.basisFunctions.touch = touchShiftIdx;
    glmModel{i}.basisFunctions.features = touchShiftIdxRaw;
    
    DmatY = spikes(startIdx(:));
    rawAngle = angle(touchIdx);
    glmModel{i}.io.DmatY = DmatY;
    glmModel{i}.raw.angle = rawAngle;
    
end