%% GLOBALS
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig7\';

%% scatter of phase x hilbert and responses
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,touchCells,'pole','off');
location_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
%% link between phase and amp vs angle
for rec = datasample(1:numel(location_units),1)
    %stimulus and response variables definitions
    array = U{location_units(rec)};
    touch_response_window = array.meta.touchProperties.responseWindow;
    [tVar] = atTouch_sorter(array,-25:50);
    
    phase = tVar.allTouches.S_ctk(:,5);
    midpoint = tVar.allTouches.S_ctk(:,4);
    amp = tVar.allTouches.S_ctk(:,3);
    angle = tVar.allTouches.S_ctk(:,1);
    
    figure(480);clf
    subplot(1,2,1);
    scatter3(phase,amp,angle,'k.')
    xlabel('phase');ylabel('amp');zlabel('angle')
    axis square
    
    subplot(1,2,2);
    scatter3(phase,midpoint,angle,'k.')
    xlabel('phase');ylabel('midpoint');zlabel('angle')
    axis square
end

%%
% location_unit 21 may be a good example
for rec = 4
    %stimulus and response variables definitions
    array = U{location_units(rec)};
    touch_response_window = array.meta.touchProperties.responseWindow;
    [tVar] = atTouch_sorter(array,-25:50);
    
    phase = tVar.allTouches.S_ctk(:,5);
    response = mean(tVar.allTouches.R_ntk(:,touch_response_window(1):touch_response_window(2)),2);
    unique_color_codes = flipud(unique(normalize_var(response,0,.9)));
    [~,idx] = sort(response);
    [unique_responses,~,unique_idx] = unique(response);
    
    figure(1230);clf
    feat_numbers = [1 3 4];
    y_names = {'angle','amp','midpoint'};
    for d = 1:3
        feature = tVar.allTouches.S_ctk(:,feat_numbers(d));
        subplot(1,3,d)
        for g = 1:numel(unique_responses)
            hold on; scatter(phase(unique_idx==g),feature(unique_idx==g),'filled','markerfacecolor',[1 1 1] .* unique_color_codes(g))
        end
        set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'})
        xlabel('phase at touch');
        ylabel([y_names(d) ' at touch'])
    end
    suptitle(['cell number : ' num2str(location_units(rec))])
    
    
end

fn = 'feature_scatter_response_2.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% phase tuning at low amp vs high amp

%amp distribution
amp_range = 1:40;
midpoint_range = -20:50;
amp_cdf = zeros(numel(amp_range),numel(location_units));
midpoint_cdf = zeros(numel(midpoint_range),numel(location_units));
phase_amp_corr = zeros(numel(location_units),1);
phase_mp_corr = zeros(numel(location_units),1);

for rec = 1:length(location_units)
    array = U{location_units(rec)};
    touch_response_window = array.meta.touchProperties.responseWindow;
    [tVar] = atTouch_sorter(array,-25:50);
    
    angle = tVar.allTouches.S_ctk(:,1); 
    phase = tVar.allTouches.S_ctk(:,5);
    amp = tVar.allTouches.S_ctk(:,3);
    amp_counts = histcounts(amp,min(amp_range)-.5:1:max(amp_range)+.5);
    amp_cdf(:,rec) = cumsum(amp_counts) ./ sum(amp_counts);
    
    midpoint = tVar.allTouches.S_ctk(:,4);
    midpoint_counts = histcounts(midpoint,min(midpoint_range)-.5:1:max(midpoint_range)+.5);
    midpoint_cdf(:,rec) = cumsum(midpoint_counts) ./ sum(midpoint_counts);
    
    phase_angle_corr(rec) = corr(phase,amp,'rows','complete');
    phase_amp_corr(rec) = corr(phase,amp,'rows','complete');
    phase_mp_corr(rec) = corr(phase,midpoint,'rows','complete');
end
figure(123);clf
subplot(1,3,1);
plot(amp_range,amp_cdf,'color',[.8 .8 .8])
hold on; plot(amp_range,mean(amp_cdf,2),'k')
set(gca,'xtick',0:10:40,'ytick',0:.5:1)
title('amplitude')
axis square

subplot(1,3,2);
plot(midpoint_range,midpoint_cdf,'color',[.8 .8 .8])
hold on; plot(midpoint_range,mean(midpoint_cdf,2),'k')
set(gca,'xtick',-20:20:50,'ytick',0:.5:1,'xlim',[-20 50])
title('midpoint')
axis square

subplot(1,3,3);
scatter(ones(numel(location_units),1),phase_amp_corr,'filled','markerfacecolor',[.8 .8 .8])
hold on; errorbar(1,mean(phase_amp_corr),std(phase_amp_corr),'ko')
hold on; scatter(ones(numel(location_units),1).*2,phase_mp_corr,'filled','markerfacecolor',[.8 .8 .8])
hold on; errorbar(2,mean(phase_mp_corr),std(phase_mp_corr),'ko')
axis square
set(gca,'xlim',[.5 2.5],'xtick',1:2,'xticklabel',{'amp','mp'})

%
fn = 'slow_feature_distribution.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%%
amp_threshold = 11;
midpoint_threshold = 11;
rc = numSubplots(numel(location_units));
figure(669);clf
figure(670);clf
figure(679);clf
figure(680);clf

for rec = 1:length(location_units)
    %stimulus and response variables definitions
    array = U{location_units(rec)};
    touch_response_window = array.meta.touchProperties.responseWindow;
    [tVar] = atTouch_sorter(array,touch_response_window(1):touch_response_window(2));
    
    phase = tVar.allTouches.S_ctk(:,5);
    amp = tVar.allTouches.S_ctk(:,3);
    midpoint = tVar.allTouches.S_ctk(:,4);
    response = mean(tVar.allTouches.R_ntk,2);
    low_amp = amp<prctile(amp,50);
    high_amp = amp>=prctile(amp,50);
    low_mp = midpoint<prctile(midpoint,50);
    high_mp = midpoint>=prctile(midpoint,50);
    
    %samples
    low_amp_samples = histcounts(phase(low_amp),linspace(-3.15,3.15,13));
    low_amp_cdf = cumsum(low_amp_samples)./sum(low_amp_samples);
    high_amp_samples = histcounts(phase(high_amp),linspace(-3.15,3.15,13));
    high_amp_cdf = cumsum(high_amp_samples)./sum(high_amp_samples);
    
    low_mp_samples = histcounts(phase(low_mp),linspace(-3.15,3.15,13));
    low_mp_cdf = cumsum(low_mp_samples)./sum(low_mp_samples);
    high_mp_samples = histcounts(phase(high_mp),linspace(-3.15,3.15,13));
    high_mp_cdf = cumsum(high_mp_samples)./sum(high_mp_samples);
    
    %response at phase
    [low_sorted] = binslin(phase(low_amp),response(low_amp) .* 1000,'equalX',13,-pi,pi);%amplitude
    [high_sorted] = binslin(phase(high_amp),response(high_amp) .* 1000,'equalX',13,-pi,pi);%amplitude
    [low_sorted_midpoint] = binslin(phase(low_mp),response(low_mp) .* 1000,'equalX',13,-pi,pi); %midpoint
    [high_sorted_midpoint] = binslin(phase(high_mp),response(high_mp) .* 1000,'equalX',13,-pi,pi);%midpoint
    
    high_mean = cellfun(@(x) nanmean(x),high_sorted);
    high_sem = cellfun(@(x) nanstd(x)./sqrt(numel(x)),high_sorted);
    low_mean = cellfun(@(x) nanmean(x),low_sorted);
    low_sem = cellfun(@(x) nanstd(x)./sqrt(numel(x)),low_sorted);
    
    high_mean_midpoint = cellfun(@(x) nanmean(x),high_sorted_midpoint);
    high_sem_midpoint = cellfun(@(x) nanstd(x)./sqrt(numel(x)),high_sorted_midpoint);
    low_mean_midpoint = cellfun(@(x) nanmean(x),low_sorted_midpoint);
    low_sem_midpoint = cellfun(@(x) nanstd(x)./sqrt(numel(x)),low_sorted_midpoint);
    
    %     figure(669); %plotting sampling at low vs high amp
    %     subplot(rc(1),rc(2),rec)
    %     hold on;bar(linspace(-pi,pi,12),low_amp_samples,'b','facealpha',1)
    %     hold on;bar(linspace(-pi,pi,12),high_amp_samples,'r','facealpha',1)
    %     yyaxis right
    %     hold on; plot(linspace(-pi,pi,12),low_amp_cdf,'b--');
    %     hold on; plot(linspace(-pi,pi,12),high_amp_cdf,'r--');
    %     set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi',0,'\pi'})
    
    
    figure(670); %plotting phase tuning at touch for low vs high amp
    subplot(rc(1),rc(2),rec)
    hold on; errorbar(linspace(-pi,pi,12),low_mean,low_sem,'bo-')
    hold on;errorbar(linspace(-pi,pi,12),high_mean,high_sem,'ro-')
    set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi',0,'\pi'})
    
    
    %     figure(679); %plotting sampling at low vs high amp
    %     subplot(rc(1),rc(2),rec)
    %     hold on;bar(linspace(-pi,pi,12),low_mp_samples,'b','facealpha',.5)
    %     hold on;bar(linspace(-pi,pi,12),high_mp_samples,'r','facealpha',.5)
    %     yyaxis right
    %     hold on; plot(linspace(-pi,pi,12),low_mp_cdf,'b--');
    %     hold on; plot(linspace(-pi,pi,12),high_mp_cdf,'r--');
    %     set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi',0,'\pi'})
    
    
    figure(680); %plotting phase tuning at touch for low vs high amp
    subplot(rc(1),rc(2),rec)
    hold on; errorbar(linspace(-pi,pi,12),low_mean_midpoint,low_sem_midpoint,'bo-')
    hold on;errorbar(linspace(-pi,pi,12),high_mean_midpoint,high_sem_midpoint,'ro-')
    set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi',0,'\pi'})
    
    
    
    
end

% figure(669)
% suptitle('amplitude')
% fn = 'phase_low_vs_high_amp_samp.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

figure(670)
suptitle('amplitude')
fn = 'phase_tuning_low_vs_high_amp.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% figure(679)
% suptitle('midpoint')
% fn = 'phase_low_vs_high_midpoint_samp.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

figure(680)
suptitle('midpoint')
fn = 'phase_tuning_low_vs_high_midpoint.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

