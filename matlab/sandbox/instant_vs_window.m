%% BUILD
hilbert_feature = {'angle','phase','midpoint','amplitude','velocity'};
capture_window = {'instant','lag','lag_window'};
naming = {'instant','lag','window'};
for k = 3
    for g = 2:5
        wStruct= whisk_location_quant_vsupp(U,1:length(U),hilbert_feature{g},'off',capture_window{k});
        cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo')
        save(['whisk_' hilbert_feature{g} '_' naming{k}],'wStruct')
    end
end
%% OR LOAD
clear whisk_struct
hilbert_feature = 'angle';
mdl_versions = {'instant','lag','window'};

for b = 1:numel(mdl_versions)
    fileName = ['whisk_' hilbert_feature '_' mdl_versions{b}];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\tmp_structs\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\tmp_structs\' fileName '.mat'])
        whisk_struct.(mdl_versions{b}) = wStruct;
    end
end

% itune = find(cellfun(@(x) x.is_tuned==1,whisk_struct.instant));
% ltune = find(cellfun(@(x) x.is_tuned==1,whisk_struct.lag));
wtune = find(cellfun(@(x) x.is_tuned==1,whisk_struct.window));
both_tuned=intersect(intersect(itune,wtune,'stable'),ltune,'stable');
%% intersect venn diagram
    A = [numel(wtune) numel(itune) numel(ltune)]; 
    I = [numel(intersect(wtune,itune)) numel(intersect(wtune,ltune)) numel(intersect(itune,ltune)) numel(both_tuned)];
    figure(1230);clf
    h2 = venn(A,I,'ErrMinMode','TotalError','FaceColor',{'r','c','g'});
    axis image,  title ('Overlap of whisking tune')
    set(gca,'xtick',[],'ytick',[])
%% PLOT TRIALS
% cellNum = datasample(intersect(itune,wtune),1);
[r] = numSubplots(numel(both_tuned));
figure(5480);clf
for g =1:numel(both_tuned)
    cellNum = both_tuned(g);
    
    iraw_x = cellfun(@median,(whisk_struct.instant{cellNum}.stim_response.raw_stim));
    iraw_y = cellfun(@mean,(whisk_struct.instant{cellNum}.stim_response.raw_response));
    ix = whisk_struct.instant{cellNum}.stim_response.bars_stim;
    iy = whisk_struct.instant{cellNum}.stim_response.bars_fit.mean;
    
    raw_x = cellfun(@median,(whisk_struct.window{cellNum}.stim_response.raw_stim));
    raw_y = cellfun(@mean,(whisk_struct.window{cellNum}.stim_response.raw_response));
    x = whisk_struct.window{cellNum}.stim_response.bars_stim;
    y = whisk_struct.window{cellNum}.stim_response.bars_fit.mean;
    
    lraw_x = cellfun(@median,(whisk_struct.lag{cellNum}.stim_response.raw_stim));
    lraw_y = cellfun(@mean,(whisk_struct.lag{cellNum}.stim_response.raw_response));
    lx = whisk_struct.lag{cellNum}.stim_response.bars_stim;
    ly = whisk_struct.lag{cellNum}.stim_response.bars_fit.mean;
    
    subplot(r(1),r(2),g)
    hold on; plot(ix,iy,'c')
    hold on; plot(x,y,'r')
    hold on; plot(lx,ly,'g')
    
        hold on; scatter(iraw_x,iraw_y,'c.');alpha(.25)
        hold on; scatter(raw_x,raw_y,'r.'); alpha(.25)
        hold on; scatter(lraw_x,lraw_y,'g.'); alpha(.25)
        
    if g == numel(both_tuned)
        legend('instantaneous','lagged window','lag')
    end
    %     xlabel('whisker angle')
    %     ylabel('firing rate')
    
    if strcmp(hilbert_feature,'phase')
        set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'xlim',[-4 4])
    end
end

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
% fn = ['instant_vs_windowed_curves.eps'];
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])


%% PREFERENCE comparison
i_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.instant(both_tuned),'uniformoutput',0)');
l_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.lag(both_tuned),'uniformoutput',0)');
w_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.window(both_tuned),'uniformoutput',0)');

figure(413);clf
subplot(1,2,1)
plot(repmat(1:3,numel(both_tuned),1)',[i_pw(:,1) l_pw(:,1) w_pw(:,1)]','ko-')
set(gca,'xlim',[.5 3.5],'xtick',1:3,'xticklabel',mdl_versions);
ylabel('angle preference')

subplot(1,2,2)
sd_pref = std([i_pw(:,1) l_pw(:,1) w_pw(:,1)],[],2);
histogram(sd_pref,0:1:50,'normalization','probability','facecolor','k','facealpha',1)
title('SD distribution between three models')
ylabel('proportion of co-tuned units')
xlabel('SD angle preference')

%% Residual SS and Log Likelihood
residual_ss = nan(numel(both_tuned),numel(mdl_versions));
for g = 1:numel(mdl_versions)
    for b = 1:numel(both_tuned)
        cellNum = both_tuned(b);
        x = whisk_struct.(mdl_versions{g}){cellNum}.stim_response.bars_fit.x;
        y = whisk_struct.(mdl_versions{g}){cellNum}.stim_response.bars_fit.mean;
        xq = cellfun(@median,(whisk_struct.(mdl_versions{g}){cellNum}.stim_response.raw_stim));
        pred_yq = interp1(x,y,xq);
        true_yq = cellfun(@mean,(whisk_struct.(mdl_versions{g}){cellNum}.stim_response.raw_response));
        residual_ss(b,g) = nansum((true_yq-pred_yq).^2);
    end
end
[~,p_t(1)] = ttest(residual_ss(:,1),residual_ss(:,2));
[~,p_t(2)] =ttest(residual_ss(:,2),residual_ss(:,3));
[~,p_t(3)] =ttest(residual_ss(:,1),residual_ss(:,3));

i_ll=cellfun(@(x) mean([x.stim_response.bars_fit.models.logLikelihood]),whisk_struct.instant(both_tuned));
l_ll=cellfun(@(x) mean([x.stim_response.bars_fit.models.logLikelihood]),whisk_struct.lag(both_tuned));
w_ll=cellfun(@(x) mean([x.stim_response.bars_fit.models.logLikelihood]),whisk_struct.window(both_tuned));

figure(480);clf
subplot(1,2,1)
plot(repmat(1:3,numel(both_tuned),1),[i_ll' l_ll' w_ll'],'ko')
hold on;errorbar(1:3,mean([i_ll' l_ll' w_ll']),std([i_ll' l_ll' w_ll'])./sqrt(numel(both_tuned)),'ro')
set(gca,'xlim',[0 4],'xtick',1:3,'xticklabel',mdl_versions,'yscale','log');
[p,~,~] = anova1([i_ll' l_ll' w_ll'],[],'off');
ylabel('bars fitting log likelihood')
title(['p-val = ' num2str(p)])

subplot(1,2,2)
plot(repmat(1:3,numel(both_tuned),1),residual_ss,'ko')
hold on;errorbar(1:3,mean(residual_ss),std(residual_ss)./sqrt(numel(both_tuned)),'ro')
set(gca,'xlim',[0 4],'xtick',1:3,'xticklabel',mdl_versions,'yscale','log');
[p,~,~] = anova1(residual_ss,[],'off');
ylabel('bars fitting SS residuals')
title(['p-val 1|2 , 2|3 , 1|3 = ' num2str(p_t)])


%% INTERSECTION OF PEAKS
i_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.instant(both_tuned),'uniformoutput',0)') ;
w_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.window(both_tuned),'uniformoutput',0)');

figure(480);clf
errorbar(i_pw(:,1),w_pw(:,1),w_pw(:,2),w_pw(:,3),i_pw(:,2),i_pw(:,3),'ko','capsize',0)
if strcmp(hilbert_feature,'angle')
    hold on; plot([-30 80],[-30 80],'--k')
    set(gca,'xlim',[-30 80],'ylim',[-30 80],'xtick',-40:20:80,'ytick',-40:20:80)
elseif strcmp(hilbert_feature,'phase')
    hold on; plot([-pi pi],[-pi pi],'--k')
    set(gca,'xlim',[-4 4],'ylim',[-4 4],'xtick',-pi:pi:pi,'ytick',-pi:pi:pi,...
        'xticklabel',{'-\pi',0,'\pi'},'yticklabel',{'-\pi',0,'\pi'})
end
xlabel('instantaneous preference')
ylabel('lagged window preference')
axis square

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = ['instant_vs_windowed_ix.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
%% OL TUNING VS WINDOW 
