%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_clean.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 

% U = defTouchResponse(U,.95,'off');

%% Calculations table 
%VARIABLE DEFINITIONS: 
masks = cellfun(@(x) maskBuilder(x),U,'uniformoutput',0);
spks_in_touch = cellfun(@(x,y) nansum(nansum(squeeze(x.R_ntk).*(isnan(y.touch)))),U,masks);
spks_in_whisking = cellfun(@(x,y) nansum(nansum(squeeze(x.R_ntk).*y.whisking .*y.touch)),U,masks);
whisking_tp = cellfun(@(x,y) nansum(nansum(y.whisking .*y.touch)),U,masks);
quiet_tp = cellfun(@(x,y) nansum(nansum(y.quiet .*y.touch)),U,masks);
spks_in_all = cellfun(@(x) nansum(squeeze(x.R_ntk(:))),U);

touchUnits = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U)~=0);
pole_tuned = object_location_quantification(U,touchUnits,'pole','off'); %for old see object_location_v1.0
tuned_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
non_tuned_units = setdiff(touchUnits,tuned_units); 

%DEPENDENT FUNCTIONS: 
glmnetOpt.buildIndices = [-25:50]; %Indices around touch
glmnetOpt.touchDirection = 'protraction';
glmModel = [];
[glmModel] = designMatrixBlocks_simplified(U(touchUnits),glmnetOpt,glmModel);
excitThresh =  num2cell(cellfun(@(x) x.meta.touchProperties.baseline(1),U(touchUnits)));
excitThresh_ol_units =  num2cell(cellfun(@(x) x.meta.touchProperties.baseline(1),U(tuned_units)));
excitThresh_non_ol_units =  num2cell(cellfun(@(x) x.meta.touchProperties.baseline(1),U(non_tuned_units)));

%CALCULATIONS: 
%non touch specific calculations
mean_fr = cellfun(@(x) mean(squeeze(x.R_ntk(:))),U)*1000; %mean firing rate of all units 
whisking_fr = (spks_in_whisking./whisking_tp)*1000; %whisking firing rate
non_whisking_fr = ((spks_in_all-spks_in_touch-spks_in_whisking)./quiet_tp)*1000; %non whisking firing rate 
prop_touch = spks_in_touch./spks_in_all; %proportion of spikes attributed to touch 
prop_whisking_touch = (spks_in_touch+spks_in_whisking)./spks_in_all;%proportion of spikes attributed to whisking + touch 

%touch specific calculations 
onset_latency = cellfun(@(x) x.meta.touchProperties.responseWindow(1),U(touchUnits)); % time of touch onset 
resp_window_length = cellfun(@(x) range(x.meta.touchProperties.responseWindow),U(touchUnits)); % length of touch response
touch_response_spks = cellfun(@(x,y) mean(x.io.DmatY),glmModel); %average number of spikes in touch response window 
pResponse_touch = cellfun(@(x,y) mean((x.io.DmatY*1000)>y),glmModel,excitThresh); %probability of generating spiking response > baseline + 95%CI
pResponse_peak = cellfun(@(x,y) mean((x.calculations.responses_at_peak)>y),pole_tuned(touchUnits),excitThresh);%probability of generating spiking response > baseline + 95%CI @ peak 
pResponse_trough = cellfun(@(x,y) mean((x.calculations.responses_at_trough)>y),pole_tuned(touchUnits),excitThresh);%probability of generating spiking response > baseline + 95%CI @ trough
Response_peak = cellfun(@(x) nanmean(x.calculations.responses_at_peak),pole_tuned(touchUnits));%spiking response
Response_trough = cellfun(@(x) nanmean(x.calculations.responses_at_trough),pole_tuned(touchUnits));%spiking response 
% Response_peak = cellfun(@(x,y) nanmean(x.calculations.responses_at_peak(x.calculations.responses_at_peak>y)),pole_tuned(touchUnits),excitThresh);%spiking response > baseline + 95%CI @ peak 
% Response_trough = cellfun(@(x,y) nanmean(x.calculations.responses_at_trough(x.calculations.responses_at_trough>y)),pole_tuned(touchUnits),excitThresh);%spiking response > baseline + 95%CI @ trough

%ol unit sig
[~,pResponse_sigs] = cellfun(@(x,y) ttest2(x.calculations.responses_at_peak>y,x.calculations.responses_at_trough>y),pole_tuned(tuned_units),excitThresh_ol_units); %
% [~,response_sigs] = cellfun(@(x,y) ttest2(x.calculations.responses_at_peak(x.calculations.responses_at_peak>y),x.calculations.responses_at_trough(x.calculations.responses_at_trough>y)),pole_tuned(tuned_units),excitThresh_ol_units); %at alpha .05
[~,response_sigs] = cellfun(@(x) ttest2(x.calculations.responses_at_peak,x.calculations.responses_at_trough),pole_tuned(tuned_units)); %at alpha .05


%non ol unit sig 
[~,pNResponse_Nsigs] = cellfun(@(x,y) ttest2(x.calculations.responses_at_peak>y,x.calculations.responses_at_trough>y),pole_tuned(non_tuned_units),excitThresh_non_ol_units); %at alpha .05
[~,Nresponse_sigs] = cellfun(@(x,y) ttest2(x.calculations.responses_at_peak(x.calculations.responses_at_peak>y),x.calculations.responses_at_trough(x.calculations.responses_at_trough>y)),pole_tuned(non_tuned_units),excitThresh_non_ol_units); %at alpha .05


%% DATA TABLE BUILDING: 
%all active units 
% all_properties = {whisking_fr,non_whisking_fr,prop_touch,prop_whisking_touch};
% mean_all = [cellfun(@mean, all_properties(~cellfun(@isempty,all_properties))) nan(1,5)]';
% std_all = [cellfun(@std, all_properties(~cellfun(@isempty,all_properties))) nan(1,5)]';
% median_all= [cellfun(@median, all_properties(~cellfun(@isempty,all_properties))) nan(1,5)]';
% [range_all(:,1)] = [cellfun(@min, all_properties(~cellfun(@isempty,all_properties))) nan(1,5)]';
% [range_all(:,2)] = [cellfun(@max, all_properties(~cellfun(@isempty,all_properties))) nan(1,5)]';

%nontouch units 
all_properties = {whisking_fr,non_whisking_fr,prop_touch,prop_whisking_touch};
non_touch = cellfun(@(x) ~strcmp(x.meta.touchProperties.responseType,'excited'),U); 
mean_nt = [cellfun(@(x) mean(x(non_touch)),all_properties) nan(1,8)]';
std_nt = [cellfun(@(x) std(x(non_touch)),all_properties) nan(1,8)]';
median_nt= [cellfun(@(x) median(x(non_touch)),all_properties) nan(1,8)]';
[range_nt(:,1)] = [cellfun(@(x) min(x(non_touch)),all_properties) nan(1,8)]';
[range_nt(:,2)] = [cellfun(@(x) max(x(non_touch)),all_properties) nan(1,8)]';

% touch units
% touch_properties = {whisking_fr(touchUnits),non_whisking_fr(touchUnits),prop_touch(touchUnits),prop_whisking_touch(touchUnits),...
%     onset_latency,resp_window_length,touch_response_spks,pResponse_touch};
% mean_touch = [cellfun(@mean, touch_properties) nan(1,1)]';
% std_touch = [cellfun(@std, touch_properties) nan(1,1)]';
% median_touch = [cellfun(@median, touch_properties) nan(1,1)]';
% range_touch(:,1) = [cellfun(@min, touch_properties) nan(1,1)]';
% range_touch(:,2) = [cellfun(@max, touch_properties) nan(1,1)]';

%tuned touch units
[~,ix_idx] = intersect(touchUnits,tuned_units);
OL_properties = {whisking_fr(tuned_units),non_whisking_fr(tuned_units),prop_touch(tuned_units),prop_whisking_touch(tuned_units),...
    onset_latency(ix_idx),resp_window_length(ix_idx),touch_response_spks(ix_idx),pResponse_touch(ix_idx),...
    pResponse_peak(ix_idx),pResponse_trough(ix_idx),Response_peak(ix_idx),Response_trough(ix_idx)};
mean_ol = [cellfun(@nanmean, OL_properties) ]';
std_ol = [cellfun(@nanstd, OL_properties) ]';
median_ol = [cellfun(@nanmedian, OL_properties) ]';
range_ol(:,1) = [cellfun(@min, OL_properties) ]';
range_ol(:,2) = [cellfun(@max, OL_properties) ]';

%nontuned touch units
non_tuned_idx = setdiff(1:numel(touchUnits),ix_idx);
non_OL_properties = {whisking_fr(non_tuned_units),non_whisking_fr(non_tuned_units),prop_touch(non_tuned_units),prop_whisking_touch(non_tuned_units),...
    onset_latency(non_tuned_idx),resp_window_length(non_tuned_idx),touch_response_spks(non_tuned_idx),pResponse_touch(non_tuned_idx),...
    pResponse_peak(non_tuned_idx),pResponse_trough(non_tuned_idx),Response_peak(non_tuned_idx),Response_trough(non_tuned_idx)};
mean_non_ol = [cellfun(@nanmean, non_OL_properties) ]';
std_non_ol = [cellfun(@nanstd, non_OL_properties) ]';
median_non_ol = [cellfun(@nanmedian, non_OL_properties) ]';
range_non_ol(:,1) = [cellfun(@min, non_OL_properties) ]';
range_non_ol(:,2) = [cellfun(@max, non_OL_properties) ]';


%% Table output
final_table = table(mean_nt,std_nt,median_nt,range_nt,... 
    mean_non_ol,std_non_ol,median_non_ol,range_non_ol,...
    mean_ol,std_ol,median_ol,range_ol);

id_properties = {'whisking (hz)','quiet (hz)','touch evoked','touch+w evoked',...
    'onset latency','response window','touch response spks','p(response) touch','p(response) peak','p(response} @ trough','response peak','response trough'}; 
final_table.Properties.RowNames = id_properties;

filename = 'spiking_characteristics_raw.xlsx';
cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
writetable(final_table,filename)

%% Histograms 
% A)
figure(20);clf
% subplot(2,4,1)
% histogram(whisking_fr,0:1:40,'facecolor',[.8 .8 .8])
% hold on; histogram(whisking_fr(touchUnits),0:1:40)
% hold on; histogram(whisking_fr(tuned_units),0:1:40)
% title('whisking firing rate')
% legend('all units','touch units','location touch units')
% 
% subplot(2,4,2);
% histogram(non_whisking_fr,0:1:30,'facecolor',[.8 .8 .8])
% hold on; histogram(non_whisking_fr(touchUnits),0:1:30)
% hold on; histogram(non_whisking_fr(tuned_units),0:1:30)
% title('non whisking firing rate')

subplot(2,6,[1 2])
scatter(non_whisking_fr,whisking_fr,'filled','markerfacecolor',[.8 .8 .8])
hold on;scatter(non_whisking_fr(touchUnits),whisking_fr(touchUnits),'filled','r')
hold on; scatter(non_whisking_fr(tuned_units),whisking_fr(tuned_units),'filled','g')
hold on; plot([0 50],[0 50],'--k')
set(gca,'xlim',[0 50],'ylim',[0 50])
axis square
xlabel('non whisking FR');ylabel('whisking FR')

subplot(2,6,3)
histogram(prop_touch,0:.05:1,'facecolor',[.8 .8 .8])
hold on; histogram(prop_touch(touchUnits),0:.05:1)
hold on; histogram(prop_touch(tuned_units),0:.05:1)
title('proportion spikes in touch window')
legend('all units','touch units','location touch units')
subplot(2,6,4)
histogram(prop_whisking_touch,0:.05:1,'facecolor',[.8 .8 .8])
hold on; histogram(prop_whisking_touch(touchUnits),0:.05:1)
hold on; histogram(prop_whisking_touch(tuned_units),0:.05:1)
title('proportion spikes in touch+whisking window')

[~,touch_tuned_idx] = intersect(touchUnits,tuned_units);
subplot(2,6,5)
hold on; histogram(onset_latency,0:1:35,'facecolor','r')
hold on; histogram(onset_latency(touch_tuned_idx),0:1:35,'facecolor','g')
title('onset latency (ms)')

subplot(2,6,6);
hold on; histogram(resp_window_length,0:1:40,'facecolor','r')
hold on; histogram(resp_window_length(touch_tuned_idx),0:1:40,'facecolor','g')
title('response window duration')

subplot(2,6,7)
hold on; histogram(touch_response_spks,0:.05:4,'facecolor','r')
hold on; histogram(touch_response_spks(touch_tuned_idx),0:.05:4,'facecolor','g')
title('spikes in touch response window')

subplot(2,6,8)
hold on; histogram(pResponse_touch,0:.05:1,'facecolor','r')
hold on; histogram(pResponse_touch(touch_tuned_idx),0:.05:1,'facecolor','g')
title('probability of touch response')


pResponse_ol_trough = pResponse_trough(ix_idx);
pResponse_ol_peak = pResponse_peak(ix_idx);
response_ol_peak = Response_peak(ix_idx);
response_ol_trough = Response_trough(ix_idx);
[~,p,~,stat] = ttest(response_ol_peak,response_ol_trough)

subplot(2,6,[9 10])
scatter(pResponse_ol_trough(pResponse_sigs<=.01),pResponse_ol_peak(pResponse_sigs<=.01),'filled','markerfacecolor','g')
hold on; scatter(pResponse_ol_trough(pResponse_sigs>.01),pResponse_ol_peak(pResponse_sigs>.01),'go')
x = mean(pResponse_ol_trough) ;
y = mean(pResponse_ol_peak) ; 
xerr = std(pResponse_ol_trough) ; 
yerr = std(pResponse_ol_peak); 
hold on; errorbar(x,y,yerr,yerr,xerr,xerr,'ko','capsize',0)
hold on; plot([0 1],[0 1],'--k')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.5:1,'ytick',0:.5:1)
axis square
xlabel('p(response) trough');ylabel('p(response) peak')
[~,p,~,stat] = ttest(pResponse_ol_peak,pResponse_ol_trough)


subplot(2,6,[11 12])
scatter(response_ol_trough(response_sigs<=.01),response_ol_peak(response_sigs<=.01),'filled','markerfacecolor','g')
hold on; scatter(response_ol_trough(response_sigs>.01),response_ol_peak(response_sigs>.01),'go')
x = nanmean(response_ol_trough) ;
y = nanmean(response_ol_peak) ; 
xerr = nanstd(response_ol_trough) ; 
yerr = nanstd(response_ol_peak); 
hold on; errorbar(x,y,yerr,yerr,xerr,xerr,'ko','capsize',0)
hold on; plot([0 150],[0 150],'--k')
set(gca,'xlim',[0 150],'ylim',[0 150],'xtick',0:50:150,'ytick',0:50:150)
axis square
xlabel('trough response');ylabel('peak response')

    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\table1\';
    fn = 'population_spiking_characteristics.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])

