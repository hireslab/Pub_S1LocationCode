%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells

%% Defining touch response
U = defTouchResponse(U,.99,'off');
touchWindow = [-25:50]; %window for analyses around touch

preDecisionTouches = preDecisionTouchMat(U);
 
for rec = 1:length(U)
 [tVar] = atTouch_sorter(U{rec},touchWindow,preDecisionTouches{rec});
 
 velocity = tVar.allTouches.S_ctk(:,2);
 if isfield(U{rec}.meta,'responseWindow')
%      tResponse = U{rec}.meta.responseWindow+find(touchWindow==0) ;
     tResponse = [-25 -5] + find(touchWindow==0); %control to find off touch responses
 else
     tResponse = [5 35] + find(touchWindow==0);
 end
  response = mean(tVar.allTouches.R_ntk(:,tResponse(1):tResponse(2)),2);
 
 [sorted,sby] = binslin(velocity,response,'equalN',10);
 
 vel.R_ntk{rec} = sorted;
 vel.S_ctk{rec} = sby;
 
end

%%
responses = cellfun(@(x) cellfun(@(y) nanmean(y),x),vel.R_ntk,'uniformoutput',0);
SEM = cellfun(@(x) cellfun(@(y) nanstd(y) ./ sqrt(numel(y)),x),vel.R_ntk,'uniformoutput',0);
stimulus = cellfun(@(x) cellfun(@(y) nanmedian(y),x),vel.S_ctk,'uniformoutput',0);

figure(8);clf
for g = 1:length(responses)
    subplot(6,10,g)
    if isfield(U{g}.meta,'responseWindow')
        shadedErrorBar(stimulus{g},responses{g}*1000,SEM{g}*1000,'b');
    else
        shadedErrorBar(stimulus{g},responses{g}*1000,SEM{g}*1000,'k');
    end
    hold on;
    plot([0 0],[0 (max(responses{g}*1000)+.01)*1.2],'-.k')
    set(gca,'ylim',[0 (max(responses{g}*1000)+.01)*1.2])
end
    