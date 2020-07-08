% This function will find the indices of touch, the spikes around that
% touch given by the RANGE (ie -25:50ms) you provide. Lastly it'll find the
% theta, phase, amplitude, setpoint, max kappa, and pre touch velocity.
%
% Simple organization tool for visualizing parameters at touch onset
%
% OUTPUT: First 6 columns will be values for the variables 
% 1) THETA 2) AMP 3) SETPOINT 4) PHASE 5) MAX KAPPA 6) PRE TOUCH VELOCITY 
% Last columns will be the spikes around your given window 
% Will output values for ALL TOUCHES, PREDECISION TOUCHES, and POST
% DECISION TOUCHES

%INPUT: vector with time points you want to view (ie [-25:50])


function [tVar] = atTouch_sorter(array,range,preDecisionMat)

if (nargin < 3) 
    preDecisionMat = []; 
end


%find touch indices 
touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
spikes = squeeze(array.R_ntk);

touchOnIdx = sort(touchOnIdx(touchOnIdx<(numel(spikes)-range(end))));
touchOffIdx = sort(touchOffIdx(1:length(touchOnIdx)));

if ~isempty(preDecisionMat)
    pretOnIdx = find(preDecisionMat==1);
    postOnIdx = setdiff(touchOnIdx,pretOnIdx);
    
    pretOffIdx = touchOffIdx(ismember(touchOnIdx,pretOnIdx));
    posttOffIdx =  setdiff(touchOffIdx,pretOffIdx);
    
    tOnIndices = {touchOnIdx,pretOnIdx,postOnIdx};
    tOffIndices = {touchOffIdx,pretOffIdx,posttOffIdx};
else
    tOnIndices = {touchOnIdx};
    tOffIndices = {touchOffIdx};
end

%R_ntk spikes all around touches 
spikesIdx = cellfun(@(x) x+repmat(range,length(x),1),tOnIndices,'uniformoutput',0);
spikesAligned = cellfun(@(x) spikes(x),spikesIdx,'uniformoutput',0);

%S_ctk stimulus mat for features around touches 
varDesign = cell(1,3);
dtDesign = cell(1,3); 
dtDesign_R_ntk = cell(1,3); 

for g = 1:length(tOnIndices)
    touchTnums = ceil(tOnIndices{g}./array.t);
    varDesign{g}(:,1) = array.S_ctk(1,tOnIndices{g});%theta at touch
    varDesign{g}(:,3) = array.S_ctk(3,tOnIndices{g});%amp at at touch
    varDesign{g}(:,4) = array.S_ctk(4,tOnIndices{g});%setpoint at touch
    varDesign{g}(:,5) = array.S_ctk(5,tOnIndices{g});%phase at touch
    try
        varDesign{g}(:,6) = array.S_ctk(17,tOnIndices{g});%curvature at touch
    catch 
        varDesign{g}(:,6) = nan(numel(array.S_ctk(5,tOnIndices{g})),1);
    end
    for i = 1:length(tOnIndices{g})
        varDesign{g}(i,2)=mean(array.S_ctk(2,tOnIndices{g}(i)-5:tOnIndices{g}(i)-1)); %finds mean of velocity (-4:-1ms) before touch
        % forcing protraction and retraction calculation of dkappa 
        try %adding this since andrew's data doesnt have these variables
            kwin=array.S_ctk(19,tOnIndices{g}(i):tOffIndices{g}(i)); %get values in touch window
            if varDesign{g}(i,2)>0 && varDesign{g}(i,5)<0 %protraction
                [minValue ,~] = min(kwin);
                dtDesign{g}(i,2)=minValue;
            elseif varDesign{g}(i,2)<0 && varDesign{g}(i,5)>0 %retraction
                [maxValue ,~] = max(kwin);
                dtDesign{g}(i,2)=maxValue;
            else
                [~ ,maxidx] = max(abs(kwin)); %find idx of max kappa within each touch window, neg or pos
                dtDesign{g}(i,2)=kwin(maxidx); %use idx to pull max kappa
            end
            
            
            ktheta=array.S_ctk(18,tOnIndices{g}(i):tOffIndices{g}(i)); %get values in touch window
            [~ ,maxidx] = max(abs(ktheta)); %find idx of max theta within each touch window, neg or pos
            dtDesign{g}(i,3)=ktheta(maxidx); %use idx to pull max theta
            
            dtDesign_R_ntk{g}(i) = sum(spikes(tOnIndices{g}(i):tOffIndices{g}(i))); %sum spikes during touch duration
        catch 
            dtDesign{g}(i,1:3) = nan(1,3);
            dtDesign_R_ntk{g}(i) = sum(spikes(tOnIndices{g}(i):tOffIndices{g}(i)));

        end
    end
    try
        varDesign{g}(:,7) = array.meta.motorPosition(touchTnums);
    catch
        varDesign{g}(:,7)  = nan(1,numel(touchTnums));
    end
    dtDesign{g}(:,1) = tOffIndices{g} - tOnIndices{g};
end

%composing output structure
tVar.allTouches.itNames = {'theta','pre-touch velocity','amp','midpoint','phase','K','motorPosition'};
tVar.allTouches.dtNames = {'touchDuration','max dK','max dTheta'};
tVar.allTouches.S_ctk = varDesign{1};  
tVar.allTouches.R_ntk = spikesAligned{1};
tVar.allTouches.dtS_ctk = dtDesign{1};
tVar.allTouches.dtR_ntk = dtDesign_R_ntk{1}; %these are all the spikes that occur DURING the whole touch window
tVar.allTouches.touchOnIdx = touchOnIdx;

if ~isempty(preDecisionMat)
    tVar.preDecisionTouches.itNames = {'theta','pre-touch velocity','amp','midpoint','phase','K','motorPosition'};
    tVar.preDecisionTouches.dtNames = {'touchDuration','max dK','max dTheta'};
    tVar.preDecisionTouches.S_ctk = varDesign{2};
    tVar.preDecisionTouches.R_ntk = spikesAligned{2};
    tVar.preDecisionTouches.dtS_ctk = dtDesign{2};
    tVar.preDecisionTouches.dtR_ntk = dtDesign_R_ntk{2};
    
    tVar.postDecisionTouches.itNames = {'theta','pre-touch velocity','amp','midpoint','phase','K','motorPosition'};
    tVar.postDecisionTouches.dtNames = {'touchDuration','max dK','max dTheta'};
    tVar.postDecisionTouches.S_ctk = varDesign{3};
    tVar.postDecisionTouches.R_ntk = spikesAligned{3};
    tVar.postDecisionTouches.dtS_ctk = dtDesign{3};
    tVar.preDecisionTouches.dtR_ntk = dtDesign_R_ntk{3};
end


