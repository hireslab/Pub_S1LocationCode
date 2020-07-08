function [ff_bin] = poisson_sampling_justification(Uarray,tune_struct)

selectedCells = find(cellfun(@(x) x.is_tuned==1,tune_struct));

pcheck = zeros(numel(selectedCells),2);
ff_group = [];
ff_bin = cell(1,numel(selectedCells));
for rec = 1:length(selectedCells)
    %stimulus and response variables definitions
    array = Uarray{selectedCells(rec)};
    touch_response_window = array.meta.touchProperties.responseWindow(1):array.meta.touchProperties.responseWindow(2);
    [tVar] = atTouch_sorter(array,touch_response_window);
    response_period = range(array.meta.touchProperties.responseWindow)+1; 
    
    %raw data
    y = sum(tVar.allTouches.R_ntk,2); %firing rate in response window
    x = normalize_var(tVar.allTouches.S_ctk(:,end),-1,1); %normalized pole pos (-1 far , 1 close)   
    pcheck(rec,:) = [mean(y) var(y)];
    
    response_spks = cellfun(@(x) (x./1000).*response_period,tune_struct{selectedCells(rec)}.stim_response.raw_response,'uniformoutput',0);
    ff_bin{rec}.fano_factor = cellfun(@(x) var(x)./mean(x),response_spks);
    ff_bin{rec}.num_spks = cellfun(@(x) mean(x),response_spks);
    
    %line of best fit and regressing out 
%     tc_x = tune_struct{selectedCells(rec)}.stim_response.values(:,1);
%     tc_y = tune_struct{selectedCells(rec)}.stim_response.values(:,2);
%     
    %cant use this because if we do the mean of residuals it is always 0,
    %assuming it is line of best fit. 
%     f = fit(tc_x,tc_y,'smoothingspline','SmoothingParam',.999);
%     figure(8);clf
%     hold on; plot(f,tc_x,tc_y);
%     pause
%     predicted = f(x);
%     residuals = y - predicted;
%     figure(9);clf
%     scatter(x,predicted)
%     hold on; scatter(x,y)
%     pause
%     pcheck(rec,:) = [mean(residuals) var(residuals)];
    
    ff_group.fano_factor(rec)  = var(y)./mean(y);
    ff_group.num_spks(rec)  = mean(y);
end

all_bin_ff = cell2mat(cellfun(@(x) x.fano_factor,ff_bin,'uniformoutput',0)');
all_bin_fr = cell2mat(cellfun(@(x) x.num_spks,ff_bin,'uniformoutput',0)');

mean_bin_ff = cell2mat(cellfun(@(x) nanmean(x.fano_factor),ff_bin,'uniformoutput',0)');
mean_bin_fr = cell2mat(cellfun(@(x) nanmean(x.num_spks),ff_bin,'uniformoutput',0)');

mega_mean_FF = mean(mean_bin_ff);
mega_SEM_FF = std(mean_bin_ff);
mega_mean_FR = mean(mean_bin_fr);
mega_SEM_FR = std(mean_bin_fr);

figure(80);clf
scatter(all_bin_fr,all_bin_ff,'k.')
hold on; scatter(mean_bin_fr,mean_bin_ff,'filled','r')%
hold on; errorbar(mega_mean_FR,mega_mean_FF,mega_SEM_FF,mega_SEM_FF,mega_SEM_FR,mega_SEM_FR,'ko','capsize',0)
% hold on; scatter(ff_group.num_spks,ff_group.fano_factor,'filled','b')%FF not accounting for pole position
set(gca,'ylim',[0 3],'ytick',0:.5:3,'xlim',[0 4],'xtick',0:.5:4)
xlabel('Mean spike count');ylabel('Fano factor')
axis square
legend('FF, location binned response','FF, all "location regressed" response')
% legend('FF, location binned response','FF, all "location regressed" response','FF, all response')
% 
% figure(48);clf
% scatter(pcheck(:,1),pcheck(:,2),'k')
% hold on; plot([0 max(pcheck,[],'all')],[0 max(pcheck,[],'all')],'--k')
% set(gca,'xlim',[0 max(pcheck,[],'all')],'ylim',[0 max(pcheck,[],'all')])
% axis square
% xlabel('mean spikes in response window')
% ylabel('variance spikes in response window')
% title(['fano factor = ' num2str(mean(fano_factor_group)) '+/-' num2str(std(fano_factor_group))])