%% GLOBALS
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';

%% STIMULUS CORRELATION
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
corr_mat = zeros(7, 7);
for rec = 1:length(touchCells)
    %stimulus and response variables definitions
    array = U{touchCells(rec)};
    [tVar] = atTouch_sorter(array,-25:50);
    it_features = tVar.allTouches.S_ctk(:,[7 1 5 3 4]);
    dt_features = tVar.allTouches.dtS_ctk(:,2:3);
    all_features = [it_features dt_features];
    all_features(any(isnan(all_features),2),:) = [];
    
    corr_mat = corr_mat + corr(all_features);
end

figure(430);clf
imagesc(abs(corr_mat./numel(touchCells)))
set(gca,'xtick',1:7,'xticklabel',{'pole','angle','phase','amp','midpoint','dkappa','dtheta'})
axis square
colormap gray
caxis([0 1])

fn = 'stimulus_correlation.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
%% TOUCH and FIRING RATE
% U = defTouchResponse(U,.95,'off');
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
mod_idx_touch = cellfun(@(x) x.meta.touchProperties.mod_idx_relative,U(touchCells));
firing_rate = log10(cellfun(@(x) mean(x.R_ntk(:))*1000,U));

%% POLE + IT FEATURES
pole_tuned = object_location_quantification(U,touchCells,'pole','off');
angle_tuned = object_location_quantification(U,touchCells,'angle','off');
phase_tuned = object_location_quantification(U,touchCells,'phase','off');
amp_tuned = object_location_quantification(U,touchCells,'amplitude','off');
mp_tuned = object_location_quantification(U,touchCells,'midpoint','off');

tunedCells = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
nontunedCells = setdiff(touchCells,tunedCells);

mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_tuned(touchCells));
mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,angle_tuned(touchCells));
mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,phase_tuned(touchCells));
mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,amp_tuned(touchCells));
mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,mp_tuned(touchCells));

peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_tuned(touchCells)) * -1;
peak_idx_angle = cellfun(@(x) x.calculations.tune_peak,angle_tuned(touchCells));
peak_idx_phase = cellfun(@(x) x.calculations.tune_peak,phase_tuned(touchCells));
peak_idx_amp = cellfun(@(x) x.calculations.tune_peak,amp_tuned(touchCells));
peak_idx_midpoint = cellfun(@(x) x.calculations.tune_peak,mp_tuned(touchCells));

% HSV PLOTTING

% [~,sort_idx] = sort(peak_idx_pole);
%
% peaks = {peak_idx_pole,peak_idx_angle,peak_idx_phase,peak_idx_amp,peak_idx_midpoint};
% mod_idx = {mod_idx_pole,mod_idx_angle,mod_idx_phase,mod_idx_amp,mod_idx_midpoint};
% final_image = nan(numel(peaks),numel(selectedCells),3);
% for b = 1:numel(peaks)
%     hues = normalize_var(peaks{b}(sort_idx),.7,1);
%     hsv = [hues' mod_idx{b}(sort_idx)' ones(numel(hues),1)];
%     rgb_values = hsv2rgb(hsv);
%     final_image(b,:,1) = rgb_values(:,1);
%     final_image(b,:,2) = rgb_values(:,2);
%     final_image(b,:,3) = rgb_values(:,3);
% end
%
% figure(580);clf
% imshow(final_image)
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'hsv_touch.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
%% EMERGENCE OF TUNING THAT GOES BEYOND PHASE (pushing boundary laid by Curtis and Kleinfeld 2009) 
% pushes boundary because Curtis and Kleinfeld show ONLY phase at touch
% modulation and relatively low angle at touch modulation. THUS L5b must be
% performing some deeper level of computation/merging of other features to reveal angle at touch
% tuning. 
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,touchCells,'pole','off');
angle_tuned = object_location_quantification(U,touchCells,'angle','off');
phase_tuned = object_location_quantification(U,touchCells,'phase','off');
amp_tuned = object_location_quantification(U,touchCells,'amplitude','off');
mp_tuned = object_location_quantification(U,touchCells,'midpoint','off');

phase_tuned_cells = ismember(touchCells,find(cellfun(@(x) x.is_tuned==1,phase_tuned)));

mod_depth_pole = cellfun(@(x) x.calculations.mod_depth,pole_tuned(touchCells));
mod_depth_angle = cellfun(@(x) x.calculations.mod_depth,angle_tuned(touchCells));
mod_depth_phase = cellfun(@(x) x.calculations.mod_depth,phase_tuned(touchCells));
mod_depth_amp = cellfun(@(x) x.calculations.mod_depth,amp_tuned(touchCells));
mod_depth_mp = cellfun(@(x) x.calculations.mod_depth,mp_tuned(touchCells));

mod_idx_abs_pole = cellfun(@(x) x.calculations.mod_idx_abs,angle_tuned(touchCells));
mod_idx_abs_angle = cellfun(@(x) x.calculations.mod_idx_abs,angle_tuned(touchCells));
mod_idx_abs_phase = cellfun(@(x) x.calculations.mod_idx_abs,phase_tuned(touchCells));
mod_idx_abs_amp = cellfun(@(x) x.calculations.mod_idx_abs,amp_tuned(touchCells));
mod_idx_abs_mp = cellfun(@(x) x.calculations.mod_idx_abs,mp_tuned(touchCells));

abs_against = {mod_idx_abs_angle,mod_idx_abs_phase,mod_idx_abs_amp,mod_idx_abs_mp};
depth_against = {mod_depth_angle,mod_depth_phase,mod_depth_amp,mod_depth_mp};
labels = {'angle','amp','mp'};
figure(480);clf
for g = 1:length(abs_against)
    subplot(2,3,g)
    scatter(mod_depth_pole,depth_against{g},'k')
    hold on; scatter(mod_depth_pole(phase_tuned_cells),depth_against{g}(phase_tuned_cells),'filled','k')
    hold on; plot([0 10],[0 10],'--k')
    set(gca,'xlim',[0 4],'ylim',[0 4])
    axis square
    xlabel('mod depth phase')
    ylabel(['mod depth ' labels{g}])
    [~,p_rel] = ttest(mod_depth_pole,depth_against{g});
    title(num2str(p_rel))
    
    subplot(2,3,g+3)
    scatter(mod_idx_abs_phase,abs_against{g},'k')
    hold on; scatter(mod_idx_abs_phase(phase_tuned_cells),abs_against{g}(phase_tuned_cells),'filled','k')
    hold on; plot([.1 150],[.1 150],'--k')
    set(gca,'xlim',[1 150],'ylim',[1 150],'xscale','log','yscale','log')
    axis square
    xlabel('mod idx abs phase')
    ylabel(['mod idx abs ' labels{g}])
    [~,p_abs] = ttest(mod_idx_abs_phase,abs_against{g});
    title(num2str(p_abs));
end

fn = 'phase_vs_feat_modulation.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


group = {mod_depth_angle,mod_depth_amp,mod_depth_mp};

figure(452);clf
for g = 1:length(group)
hold on; scatter(group{g}./mod_depth_phase,ones(numel(group{g}),1).*g,'filled','markerfacecolor',[.8 .8 .8])
hold on; errorbar(nanmean(group{g}./mod_depth_phase),g,nanstd(group{g}./mod_depth_phase),'horizontal','ro')
end
hold on; plot([1 1],[0 4],'k--')
set(gca,'ytick',1:3,'yticklabel',{'angle','amp','midpoint'},'ylim',[.5 3.5])
xlabel('modulation ratio over phase')

fn = 'mod_ratio.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% DYNAMIC FEATURES

dkappa_tuned = dynamic_touch_quantification(U,touchCells,'dkappa','off');
dtheta_tuned = dynamic_touch_quantification(U,touchCells,'dtheta','off');

mod_idx_dk = cellfun(@(x) x.calculations.mod_idx_relative,dkappa_tuned(touchCells));
mod_idx_dt = cellfun(@(x) x.calculations.mod_idx_relative,dtheta_tuned(touchCells));


%% ADAPTATION
[adaptation] = adaptation_quantification(U,touchCells,'on');
% [adaptation] = adaptation_quantification(U,tunedCells,'off');
mod_idx_adaptation = cellfun(@(x) x.calculations.mod_idx_relative,adaptation(touchCells));

% TOUCH ORDER HEAT
population_heat_tuned = normalize_var(cell2mat(cellfun(@(x) x.lh,adaptation(tunedCells),'uniformoutput',0)')',0,1);
population_heat_nontuned = normalize_var(cell2mat(cellfun(@(x) x.lh,adaptation(nontunedCells),'uniformoutput',0)')',0,1);
[~,idx] = sort(population_heat_tuned(1,:));
[~,idx_non] = sort(population_heat_nontuned(1,:));
figure(81);clf
subplot(1,2,1)
imagesc(population_heat_tuned(:,fliplr(idx)))
caxis([0 1])
title('location tuned')
subplot(1,2,2);
imagesc(population_heat_nontuned(:,fliplr(idx_non)))
caxis([0 1])
title('non-location tuned')
colormap gray
xlabel('cell number')
ylabel('touch order')
colorbar
set(gca,'ydir','reverse')
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'adaptation_map.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

%% PROTRACTION VS RETRACTION

alpha_value = 0.05;
dir_mod = touch_directional_selectivity(U,touchCells,alpha_value,'off');
mod_idx_directional = cellfun(@(x) x.mod_idx_relative,dir_mod(touchCells));

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'pro_vs_ret_unity_scatter.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])


%% PLOTTING TUNED UNITS vs UNTUNED

touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tunedCells = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
nontunedCells = setdiff(touchCells,tunedCells);

tuned_touch_idx = find(ismember(touchCells,tunedCells));
nontuned_touch_idx = find(ismember(touchCells,nontunedCells));

[~,idx] = sort(mod_idx_pole);
tunedidx = ismember(idx,tuned_touch_idx);
nontunedidx = ismember(idx,nontuned_touch_idx);

tune_map_1 = [mod_idx_touch ; mod_idx_adaptation];
tune_map_it = [mod_idx_pole ; mod_idx_angle ; mod_idx_phase ; mod_idx_amp ; mod_idx_midpoint];
tune_map_dt = [mod_idx_dk ; mod_idx_dt];
tune_map_fr = firing_rate(touchCells);

figure(81);clf
subplot(4,2,1);
imagesc(tune_map_1(:,idx(tunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'touch','adaptation'})
title('tuned')
caxis([0 1])
subplot(4,2,3);
imagesc(tune_map_it(:,idx(tunedidx)))
set(gca,'ytick',1:6,'yticklabel',{'pole','angle','phase','amp','midpoint'})
caxis([0 1])
subplot(4,2,5);
imagesc(tune_map_dt(:,idx(tunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'max dkappa','max dtheta'})
caxis([0 1])
subplot(4,2,7);
imagesc(tune_map_fr(idx(tunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'log firing rate'})
caxis([-2 2])

subplot(4,2,2);
imagesc(tune_map_1(:,idx(nontunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'touch','adaptation'})
title('nontuned')
caxis([0 1])
subplot(4,2,4);
imagesc(tune_map_it(:,idx(nontunedidx)))
set(gca,'ytick',1:6,'yticklabel',{'pole','angle','phase','amp','midpoint'})
caxis([0 1])
subplot(4,2,6);
imagesc(tune_map_dt(:,idx(nontunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'max dkappa','max dtheta'})
caxis([0 1])
subplot(4,2,8);
imagesc(tune_map_fr(idx(nontunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'log firing rate'})
caxis([-2 2])

colormap(gray)
colorbar

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map_split.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%PROTRACTION RETRACTION TUNING
figure(91);clf
tune_mod_idx_directional = mod_idx_directional(tune_idx);
nontune_mod_idx_directional = mod_idx_directional(non_tune_idx);
subplot(1,2,1);
imagesc(tune_mod_idx_directional(tune_idx))
caxis([-1 1])
subplot(1,2,2);
imagesc(nontune_mod_idx_directional(non_tune_idx))
caxis([-1 1])

stretch_resolution = 500 ;
stretch_map = redbluecmap;
new_map = nan(stretch_resolution,3);
for j = 1:3
    new_map(:,j) = interp1(linspace(1,stretch_resolution,length(stretch_map)),stretch_map(:,j),1:stretch_resolution);
end
colormap(new_map)
caxis([-1 1])
colorbar
set(gca,'ytick',1:5,'yticklabel',{'direction'})

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map_direction_split.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])






%% WHISKING FEATURE TUNING

%whisk x quiet
masks = cellfun(@(x) maskBuilder(x),U,'uniformoutput',0);
whisking_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.whisking .*y.touch,U,masks,'uniformoutput',0);
quiet_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.quiet .*y.touch,U,masks,'uniformoutput',0);
whisking_tp = cellfun(@(x,y) nansum(nansum(y.whisking .*y.touch)),U,masks);
quiet_tp = cellfun(@(x,y) nansum(nansum(y.quiet .*y.touch)),U,masks);

[~,p] = cellfun(@(x,y) ttest2(x(:),y(:)),whisking_spks_mat,quiet_spks_mat);
fr_whisk = (cellfun(@(x) nansum(x(:)),whisking_spks_mat)./whisking_tp)*1000;
fr_quiet = (cellfun(@(x) nansum(x(:)),quiet_spks_mat)./quiet_tp)*1000;

red_dots = intersect(find(p<.05),find(fr_whisk>fr_quiet));
blue_dots = intersect(find(p<.05),find(fr_whisk<fr_quiet));
gray_dots = setdiff(1:numel(masks),[red_dots blue_dots]);

figure(480);clf
scatter(fr_quiet(red_dots),fr_whisk(red_dots),'filled','r')
hold on; scatter(fr_quiet(blue_dots),fr_whisk(blue_dots),'filled','b')
hold on; scatter(fr_quiet(gray_dots),fr_whisk(gray_dots),'filled','markerfacecolor',[.8 .8 .8])
axis square
hold on; plot([0 max([fr_quiet fr_whisk])],[0 max([fr_quiet fr_whisk])],'--k')
set(gca,'xlim',[0 max([fr_quiet fr_whisk])],'ylim',[0 max([fr_quiet fr_whisk])])
xlabel('quiet FR');ylabel('whisking FR')
title(['red=' num2str(numel(red_dots)) ' blue=' num2str(numel(blue_dots)) ' gray=' num2str(numel(gray_dots))])


saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
fn = 'whisk_quiet_unity.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%%
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
nontouchCells = find(cellfun(@(x) ~strcmp(x.meta.touchProperties.responseType,'excited'),U));
tunedCells = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
nontunedCells = setdiff(touchCells,tunedCells);

selectedCells = 1:length(U);
% pole_whisk = whisking_location_quantification(U,selectedCells,'pole','off');
angle_whisk = whisking_location_quantification(U,selectedCells,'angle','off');
midpoint_whisk = whisking_location_quantification(U,selectedCells,'midpoint','off');
amp_whisk = whisking_location_quantification(U,selectedCells,'amplitude','off');
phase_whisk = whisking_location_quantification(U,selectedCells,'phase','off');

% w_mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_whisk(selectedCells));
w_mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,angle_whisk(selectedCells));
w_mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,phase_whisk(selectedCells));
w_mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,amp_whisk(selectedCells));
w_mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,midpoint_whisk(selectedCells));

% w_peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_whisk(selectedCells)) * -1;
w_peak_idx_angle = cellfun(@(x) x.calculations.tune_peak,angle_whisk(selectedCells));
w_peak_idx_phase = cellfun(@(x) x.calculations.tune_peak,phase_whisk(selectedCells));
w_peak_idx_amp = cellfun(@(x) x.calculations.tune_peak,amp_whisk(selectedCells));
w_peak_idx_midpoint = cellfun(@(x) x.calculations.tune_peak,midpoint_whisk(selectedCells));


% HSV PLOTTING WHISK
% [~,sort_idx] = sort(w_peak_idx_pole);
%
% peaks = {w_peak_idx_pole,w_peak_idx_angle,w_peak_idx_phase,w_peak_idx_amp,w_peak_idx_midpoint};
% mod_idx = {w_mod_idx_pole,w_mod_idx_angle,w_mod_idx_phase,w_mod_idx_amp,w_mod_idx_midpoint};
% final_image = nan(numel(peaks),numel(selectedCells),3);
% for b = 1:numel(peaks)
%     hues = normalize_var(peaks{b}(sort_idx),.7,1);
%     hsv = [hues' mod_idx{b}(sort_idx)' ones(numel(hues),1)];
%     rgb_values = hsv2rgb(hsv);
%     final_image(b,:,1) = rgb_values(:,1);
%     final_image(b,:,2) = rgb_values(:,2);
%     final_image(b,:,3) = rgb_values(:,3);
% end
% %
% figure(580);clf
% imshow(final_image)
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'hsv_whisk.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
% 
% % HEAT GRAY WHISK TOUCH TUNING SORTING 
% [~,idx] = sort(w_mod_idx_angle);
% nontouchidx = ismember(idx,nontouchCells);
% tunedidx = ismember(idx,tunedCells);
% nontunedidx = ismember(idx,nontunedCells);
% 
% tune_map = [w_mod_idx_angle ; w_mod_idx_phase ; w_mod_idx_amp ; w_mod_idx_midpoint];
% fr_map  =  firing_rate;
% figure(19);clf
% subplot(2,3,1)
% imagesc(tune_map(:,idx(nontouchidx)))
% set(gca,'ytick',1:6,'yticklabel',{'angle','phase','amp','midpoint'})
% caxis([0 1])
% subplot(2,3,2)
% imagesc(tune_map(:,idx(tunedidx)))
% caxis([0 1])
% subplot(2,3,3)
% imagesc(tune_map(:,idx(nontunedidx)))
% caxis([0 1])
% colormap(gray)
% colorbar
% 
% subplot(2,3,4)
% imagesc(fr_map(:,idx(nontouchidx)))
% set(gca,'ytick',1,'yticklabel',{'log firing rate'})
% caxis([0 1])
% subplot(2,3,5)
% imagesc(fr_map(:,idx(tunedidx)))
% caxis([0 1])
% subplot(2,3,6)
% imagesc(fr_map(:,idx(nontunedidx)))
% caxis([-2 2])
% colormap(gray)
% colorbar
% 
% fn = 'gray_whisk_touch.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])


% HEAT GRAY WHISK SORTING
[~,idx] = sort(w_mod_idx_angle);
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,angle_whisk);
whisktunedidx = ismember(idx,find(whisk_angle_tuned));
nonwhisktunedidx = ismember(idx,find(~whisk_angle_tuned));
tune_map = [w_mod_idx_angle ; w_mod_idx_phase ; w_mod_idx_amp ; w_mod_idx_midpoint];
fr_map  =  firing_rate;

figure(29);clf
subplot(2,2,1)
imagesc(tune_map(:,idx(whisktunedidx)))
set(gca,'ytick',1:6,'yticklabel',{'angle','phase','amp','midpoint'})
caxis([0 1])
subplot(2,2,2)
imagesc(tune_map(:,idx(nonwhisktunedidx)))
set(gca,'ytick',1:6,'yticklabel',{'angle','phase','amp','midpoint'})
caxis([0 1])
subplot(2,2,3)
imagesc(fr_map(:,idx(whisktunedidx)))
set(gca,'ytick',1,'yticklabel',{'log firing rate'})
caxis([-2 2])
subplot(2,2,4)
imagesc(fr_map(:,idx(nonwhisktunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'log firing rate'})
caxis([-2 2])
colormap(gray)

fn = 'gray_whisk_sort.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% figure(144);clf
% imagesc(abs(corr(tune_map')))
% colormap(gray)
% colorbar
% caxis([0 1])
% axis square
% set(gca,'xtick',[],'ytick',[])
%
% fn = 'mod_depth_correlation_whisk.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])





