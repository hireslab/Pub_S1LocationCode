%% Load whisking structures 
clear whisk_struct
hilbert_feature = {'angle','phase','midpoint','amplitude','velocity'};

for b = 1:numel(hilbert_feature)
    fileName = ['whisk_' hilbert_feature{b} '_window'];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'])
        whisk_struct.(hilbert_feature{b}) = wStruct;
   end
end
tuned_units = cellfun(@(x) x.is_tuned,whisk_struct.angle)==1;

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';

%% build whisking structures 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil

wStruct= whisking_location_quantification(U,1:length(U),'velocity','off');

%% whisk x quiet (A) 
masks = cellfun(@(x) maskBuilder(x),U,'uniformoutput',0);
whisking_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.whisking .*y.touch,U,masks,'uniformoutput',0);
quiet_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.quiet .*y.touch,U,masks,'uniformoutput',0);
whisking_tp = cellfun(@(x,y) nansum(nansum(y.whisking .*y.touch)),U,masks);
quiet_tp = cellfun(@(x,y) nansum(nansum(y.quiet .*y.touch)),U,masks);


for b = 1:numel(whisking_spks_mat)
    w = whisking_spks_mat{b}(:);
    q = quiet_spks_mat{b}(:);
    n1 = nansum(w==1); N1 = nansum(w==0);
    n2 = nansum(q==1); N2 = nansum(q==0);

    %generate sample data using proportions in data since matlab chi square
    %requires vectors to be equal lengths.
    x1 = [repmat('a',N1,1); repmat('b',N2,1)];
    x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
    %chi squared calcluation
    [tbl,chi2stat,p(b)] = crosstab(x1,x2); %

end
    
% [~,p] = cellfun(@(x,y) ttest2(x(:),y(:)),whisking_spks_mat,quiet_spks_mat);
fr_whisk = (cellfun(@(x) nansum(x(:)),whisking_spks_mat)./whisking_tp)*1000;
fr_quiet = (cellfun(@(x) nansum(x(:)),quiet_spks_mat)./quiet_tp)*1000;

red_dots = intersect(find(p<.05),find(fr_whisk>fr_quiet));
blue_dots = intersect(find(p<.05),find(fr_whisk<fr_quiet));
gray_dots = setdiff(1:numel(masks),[red_dots blue_dots]);

figure(480);clf
scatter(fr_quiet(p<.05),fr_whisk(p<.05),'filled','k')
hold on; scatter(fr_quiet(gray_dots),fr_whisk(gray_dots),'ko')
axis square
hold on; plot([0.01 100],[0.01 100],'--k')
set(gca,'xlim',[0 100],'ylim',[0 100],'xscale','log','yscale','log')
xlabel('quiet FR');ylabel('whisking FR')
title(['red=' num2str(numel(red_dots)) ' blue=' num2str(numel(blue_dots)) ' n.s.=' num2str(numel(gray_dots))])


saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'whisk_quiet_unity.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% PLos Reviewer Edits (1.5)
low_fr = fr_whisk(fr_quiet< 5) - fr_quiet(fr_quiet<5);
high_fr = fr_whisk(fr_quiet>= 5) - fr_quiet(fr_quiet>=5);

[~,low_p] = ttest(fr_whisk(fr_quiet< 5),fr_quiet(fr_quiet<5));
[~,high_p] = ttest(fr_whisk(fr_quiet>= 5) - fr_quiet(fr_quiet>=5));


%% angle and phase tuning curves (B)
angle_tuned = find(cellfun(@(x) x.is_tuned==1,whisk_struct.angle));
rc = numSubplots(numel(angle_tuned));
tmp = cellfun(@(x) range(x.stim_response.bars_fit.mean),whisk_struct.phase(angle_tuned));
[~,idx] = sort(tmp);

figure(480);clf
check_cells = idx(end-10:end);

for b = 1:numel(check_cells)
    cellNum = angle_tuned(check_cells(b));
    ix = whisk_struct.angle{cellNum}.stim_response.bars_stim;
    iy = whisk_struct.angle{cellNum}.stim_response.bars_fit.mean;
    i_err = whisk_struct.angle{cellNum}.stim_response.bars_fit.confBands(:,2) - iy;
    x = cellfun(@median,(whisk_struct.angle{cellNum}.stim_response.raw_stim));
    y = cellfun(@mean,(whisk_struct.angle{cellNum}.stim_response.raw_response));
    
    ix_p = whisk_struct.phase{cellNum}.stim_response.bars_stim;
    iy_p = whisk_struct.phase{cellNum}.stim_response.bars_fit.mean;
    i_err_p = whisk_struct.phase{cellNum}.stim_response.bars_fit.confBands(:,2) - iy_p;
    x_p = cellfun(@median,(whisk_struct.phase{cellNum}.stim_response.raw_stim));
    y_p = cellfun(@mean,(whisk_struct.phase{cellNum}.stim_response.raw_response));
    
    
    subplot(2,numel(check_cells),b)
    hold on; bar(x,y,'k')
    shadedErrorBar(ix,iy,i_err,'c')

    subplot(2,numel(check_cells),numel(check_cells)+b)
    hold on;bar(x_p,y_p,'k')
    shadedErrorBar(ix_p,iy_p,i_err_p,'r')
    set(gca,'xlim',[-pi pi],'xtick',[])
end

fn = ['matched_tuning_curves.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Angle and phase heat map (C/D)
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.angle);
whisk_phase_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.phase);

bars_io = cellfun(@(x) [x.stim_response.bars_stim' x.stim_response.bars_fit.mean ],whisk_struct.angle(whisk_angle_tuned),'uniformoutput',0);
whisk_sampled = round([min(cellfun(@(x) min(x(:,1)),bars_io)):max(cellfun(@(x) max(x(:,1)),bars_io))]);
interped_resp = cellfun(@(x) interp1(x(:,1),x(:,2),whisk_sampled),bars_io,'uniformoutput',0);
% norm_resp = cellfun(@(x) normalize_var(x,0,1),interped_resp,'uniformoutput',0);
norm_resp = cellfun(@(x) norm_new(x),interped_resp,'uniformoutput',0);

[~,maxidx] = cellfun(@(x) max(x),norm_resp);
[~,sorted_idx] = sort(maxidx);

sort_whisk_tuning = cell2mat(norm_resp(sorted_idx)');

figure(45);clf
pcolor(cell2mat(norm_resp(sorted_idx)')); %flipud because pcolor flips compared to imagesc
% set(gca,'xtick',1:10:numel(whisk_sampled),'xticklabel',whisk_sampled(1:10:end))
set(gca,'xtick',4:10:numel(whisk_sampled),'xticklabel',whisk_sampled(4:10:end))
colormap turbo
colorbar
caxis([0 1])
xlabel('whisker position');

%squished heatmap for phase and angle
whisk_cells_idx = {find(whisk_angle_tuned),  find(whisk_phase_tuned)};

hvar_names = {'angle','phase'};
whisk_structs = {whisk_struct.angle,whisk_struct.phase};
clear whisk_heat
figure(1210);clf
for e = 1:numel(whisk_structs)
    
    for d = 1:numel(whisk_cells_idx{e})
        curr_w = whisk_structs{e}{whisk_cells_idx{e}(d)}.stim_response.values;
        
        %clean nan rows
        curr_w = curr_w(~any(isnan(curr_w),2),:);
        
        if strcmp(hvar_names{e},'pole')
            whisk_x = -2:.1:2;
        elseif strcmp(hvar_names{e},'phase')
            whisk_x = linspace(-pi,pi,21);
        elseif strcmp(hvar_names{e},'angle')
            whisk_x = -30:80;
        else
            whisk_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
        end
        whisk_heat{e}{d} = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
    end
    
    %this is the part for squishing the map for angle; 
    if strcmp(hvar_names{e},'angle')
%         whisk_stretch_bins = median(cellfun(@(x) sum(~isnan(x)),whisk_heat{e}));
        whisk_stretch_bins = 20; 
        raw_whisk = cellfun(@(x) x(~isnan(x)),whisk_heat{e},'uniformoutput',0);
        new_whisk = cellfun(@(x) interp1(linspace(1,whisk_stretch_bins,numel(x)),x,1:whisk_stretch_bins),raw_whisk,'uniformoutput',0);
        whisk_heat{e} = new_whisk;
    end
    
    full_map = norm_new(cell2mat(whisk_heat{e}')');
    full_map = [nan(1,size(full_map,2)) ; full_map ; nan(1,size(full_map,2))]; %pad edges since pcolor trims off edges. 
    [~,idx] = max(full_map);
    [~,sort_idx] = sort(idx);
    subplot(1,2,e)
    pcolor(full_map(:,sort_idx)');
    colormap turbo
    title(hvar_names{e})
    set(gca,'xtick',[]);
end

fn = 'whisker_position_pop.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% angle and phase pie chart (E)
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.angle);
whisk_phase_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.phase);
co_tuned = intersect(find(whisk_angle_tuned),find(whisk_phase_tuned)); 
angle_only = setdiff(find(whisk_angle_tuned),co_tuned);
phase_only = setdiff(find(whisk_phase_tuned),co_tuned);
untuned = setdiff(1:length(U),[co_tuned angle_only phase_only]);

figure(580620);clf
pie([numel(untuned) numel(co_tuned) numel(angle_only) numel(phase_only)],...
    {'not tuned','co-tuned','angle only','phase only'})

fn = 'proportion_angle_phase_pie.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% Normalized phase and angle tuning curves (F examples)
tuned_units = (double(cellfun(@(x) x.is_tuned,whisk_struct.angle)==1) + double(cellfun(@(x) x.is_tuned,whisk_struct.phase)==1))>0;

% phase_mod = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.phase(tuned_units));
% angle_mod = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.angle(tuned_units));
% relative_mod = (angle_mod - phase_mod) ./ (phase_mod+angle_mod); %negative value = phase more than angle
phase_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(tuned_units));
angle_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(tuned_units));
relative_mod = (angle_mod_abs - phase_mod_abs) ./ (phase_mod_abs+angle_mod_abs); %negative value = phase more than angle

[s_rel_mod,idx]=sort(relative_mod);

phase_tc = cellfun(@(x) x.stim_response.bars_fit.mean',whisk_struct.phase(tuned_units),'uniformoutput',0);
stretch_bins =  mean(cellfun(@numel,phase_tc)); 
angle_tc = cellfun(@(x) x.stim_response.bars_fit.mean,whisk_struct.angle(tuned_units),'uniformoutput',0);
d_angle_tc = cellfun(@(x) interp1(linspace(1,stretch_bins,numel(x)),x,1:stretch_bins),angle_tc,'uniformoutput',0);

phase_tc = phase_tc(idx);
d_angle_tc = d_angle_tc(idx);

figure(58012);clf
rc = numSubplots(numel(phase_tc));
for b = 1:numel(phase_tc)
    maxmin = [min([d_angle_tc{b} phase_tc{b}]) max([d_angle_tc{b} phase_tc{b}])];
    subplot(rc(1),rc(2),b)
    plot(phase_tc{b},'r')
    hold on; plot(d_angle_tc{b},'b')
    hold on; plot([stretch_bins/2 stretch_bins/2],maxmin,'--k')
    axis tight
    
    title(num2str(s_rel_mod(b)))
    set(gca,'xtick',[],'ytick',maxmin,'ylim',maxmin,'yticklabel',round(maxmin))
end
legend('phase','angle')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'phase_x_angle_curves.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% phase X angle modulation comparison (F/G)
tuned_units = (double(cellfun(@(x) x.is_tuned,whisk_struct.angle)==1) + double(cellfun(@(x) x.is_tuned,whisk_struct.phase)==1))>0;

% calculate angle modulation 
phase_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(tuned_units));
angle_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(tuned_units));
relative_mod = (angle_mod_abs - phase_mod_abs) ./ (phase_mod_abs+angle_mod_abs); %negative value = phase more than angle

% normalize angle and phase preferences
phase_tc = cellfun(@(x) x.stim_response.bars_fit.mean',whisk_struct.phase(tuned_units),'uniformoutput',0);
stretch_bins =  mean(cellfun(@numel,phase_tc)); 
angle_tc = cellfun(@(x) x.stim_response.bars_fit.mean,whisk_struct.angle(tuned_units),'uniformoutput',0);
d_angle_tc = cellfun(@(x) interp1(linspace(1,stretch_bins,numel(x)),x,1:stretch_bins),angle_tc,'uniformoutput',0);
[~,idx] = cellfun(@max,d_angle_tc);
angle_pref = normalize_var(idx,0,1);
phase_pref = cellfun(@(x) x.calculations.tune_peak,whisk_struct.phase(tuned_units));

% define parameters for color 
num_brew_elem = 30; 
sel_map = round(normalize_var([relative_mod -1 1],1,num_brew_elem));
div_cmap = cbrewer('div','RdBu',num_brew_elem);
sel_map = sel_map(1:end-2);

% plot 
figure(88);clf
subplot(2,4,[1 2 5 6])
scatter(phase_pref',angle_pref',100,div_cmap(sel_map,:),'filled');
set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'ylim',[-.1 1.1])
axis square
xlabel('phase preference')
ylabel('norm angle preference')

sorted_angle = binslin(angle_pref',relative_mod','equalE',5,0,1);
sorted_phase = binslin(phase_pref',relative_mod','equalE',5,-pi,pi);

subplot(2,4,[3 4])
bar(linspace(0,1,4),cellfun(@mean,sorted_angle),'k')
hold on;errorbar(linspace(0,1,4),cellfun(@mean,sorted_angle),cellfun(@(x) std(x)./sqrt(numel(x)),sorted_angle),'k.')
set(gca,'ylim',[-.25 .75],'ytick',-.25:.25:.75)
subplot(2,4,[7 8])
bar(linspace(-pi,pi,4),cellfun(@mean,sorted_phase),'k')
hold on;errorbar(linspace(-pi,pi,4),cellfun(@mean,sorted_phase),cellfun(@(x) std(x)./sqrt(numel(x)),sorted_phase),'k.')
set(gca,'ylim',[-.25 .75],'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'ytick',-.25:.25:.75)

fn = 'phase_x_angle_x_preference.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% inherit cotuned from main_whisk_touch
ix = logical(tuned_units .* cotuned);
nix = logical(tuned_units .* ~cotuned);

nix_phase = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(nix));
nix_angle = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(nix)); 

ix_phase = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(ix));
ix_angle = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(ix));

figure(57);clf
% ABSOLUTE MODULATION 
scatter(phase_mod_abs,angle_mod_abs,'k')
hold on; scatter(ix_phase, ix_angle,'filled','r')
hold on; scatter(nix_phase, nix_angle, 'filled','b')
hold on; plot([0 30], [0 30],'k--')
hold on; errorbar(mean(phase_mod_abs),mean(angle_mod_abs),...
    std(angle_mod_abs)./sqrt(sum(tuned_units)),std(angle_mod_abs)./sqrt(sum(tuned_units)),...
    std(phase_mod_abs)./sqrt(sum(tuned_units)),std(phase_mod_abs)./sqrt(sum(tuned_units))...
    ,'ko','capsize',0)

hold on; errorbar(mean(ix_phase),mean(ix_angle),...
    std(ix_angle)./sqrt(sum(tuned_units)),std(ix_angle)./sqrt(sum(tuned_units)),...
    std(ix_phase)./sqrt(sum(tuned_units)),std(ix_phase)./sqrt(sum(tuned_units))...
    ,'ro','capsize',0)

hold on; errorbar(mean(nix_phase),mean(nix_angle),...
    std(nix_angle)./sqrt(sum(tuned_units)),std(nix_angle)./sqrt(sum(tuned_units)),...
    std(nix_phase)./sqrt(sum(tuned_units)),std(nix_phase)./sqrt(sum(tuned_units))...
    ,'bo','capsize',0)

set(gca,'xlim',[0 30],'ylim',[0 30])
axis square
xlabel('Phase absolute modulation')
ylabel('Angle absolute modulation')
[~,p,~,stats] = ttest(phase_mod_abs,angle_mod_abs);
title(['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])



% MODULATION DEPTH
% scatter(phase_mod,angle_mod,'k')
% hold on; plot([0 1], [0 1],'k--')
% hold on; errorbar(mean(phase_mod),mean(angle_mod),...
%     std(angle_mod)./sqrt(sum(tuned_units)),std(angle_mod)./sqrt(sum(tuned_units)),...
%     std(phase_mod)./sqrt(sum(tuned_units)),std(phase_mod)./sqrt(sum(tuned_units))...
%     ,'ro','capsize',0)
% set(gca,'xlim',[0 1],'ylim',[0 1])
% xlabel('phase mod depth')
% ylabel('angle mod depth')
% axis square
% [~,p,~,stats] = ttest(phase_mod,angle_mod);
% title(['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])

fn = 'phase_x_angle.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% Heat of modulation depth comparison (H/I/J)
selectedCells = 1:length(wStruct);

% w_mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_whisk(selectedCells));
w_mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.angle(selectedCells));
w_mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.phase(selectedCells));
w_mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.amplitude(selectedCells));
w_mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.midpoint(selectedCells));
w_mod_idx_velocity = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.velocity(selectedCells));

mod_idx_collated = [w_mod_idx_angle;w_mod_idx_phase;w_mod_idx_amp;w_mod_idx_midpoint;w_mod_idx_velocity];

[~,p,~,stats] = ttest(mod_idx_collated(1,:),mod_idx_collated(4,:)); 

angle_abs= cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(selectedCells));
phase_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(selectedCells));
amp_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.amplitude(selectedCells));
midpoint_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.midpoint(selectedCells));
velocity_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.velocity(selectedCells));

% w_peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_whisk(selectedCells)) * -1;
w_peak_idx_angle = cellfun(@(x) x.calculations.tune_peak,whisk_struct.angle(selectedCells));
w_peak_idx_phase = cellfun(@(x) x.calculations.tune_peak,whisk_struct.phase(selectedCells));
w_peak_idx_amp = cellfun(@(x) x.calculations.tune_peak,whisk_struct.amplitude(selectedCells));
w_peak_idx_midpoint = cellfun(@(x) x.calculations.tune_peak,whisk_struct.midpoint(selectedCells));
w_peak_idx_velocity = cellfun(@(x) x.calculations.tune_peak,whisk_struct.velocity(selectedCells));

% HEAT GRAY WHISK SORTING
[~,idx] = sort(w_mod_idx_angle);
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.angle);
whisktunedidx = ismember(idx,find(whisk_angle_tuned));
nonwhisktunedidx = ismember(idx,find(~whisk_angle_tuned));
tune_map = [w_mod_idx_angle ; w_mod_idx_phase ; w_mod_idx_amp ; w_mod_idx_midpoint ; w_mod_idx_velocity];
fr_map  = log10(cellfun(@(x) mean(x.R_ntk(:))*1000,U));

figure(29);clf
subplot(2,2,1)
imagesc(tune_map(:,idx(whisktunedidx)))
set(gca,'ytick',1:6,'yticklabel',fields(whisk_struct))
caxis([0 1])
subplot(2,2,2)
imagesc(tune_map(:,idx(nonwhisktunedidx)))
set(gca,'ytick',1:6,'yticklabel',fields(whisk_struct))
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

figure(3480);clf
subplot(1,2,1)
imagesc(abs(corr([tune_map(:,idx(whisktunedidx)) ; fr_map(:,idx(whisktunedidx))]')));
caxis([.5 1])
colorbar
set(gca,'xtick',1:10,'xticklabel',{'angle','phase','amp','midpoint','velocity','firing rate'},'ytick',[]);xtickangle(45)
axis square
subplot(1,2,2)
imagesc(abs(corr([tune_map(:,idx(nonwhisktunedidx)) ; fr_map(:,idx(nonwhisktunedidx))]')));
set(gca,'xtick',1:10,'xticklabel',{'angle','phase','amp','midpoint','velocity','firing rate'},'ytick',[]);xtickangle(45)
xtickangle(45)
colormap gray
colorbar
caxis([0 1])
axis square

fn = 'modulation_correlation.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

tuned_units = find(cellfun(@(x) x.is_tuned,whisk_struct.angle));

a_p = angle_abs(tuned_units) - phase_abs(tuned_units);
a_a = angle_abs(tuned_units) - amp_abs(tuned_units);
a_m = angle_abs(tuned_units) - midpoint_abs(tuned_units);
a_v = angle_abs(tuned_units) - velocity_abs(tuned_units);

figure(384);clf
hold on; scatter(ones(1,numel(a_p)),a_p,'k')
hold on; scatter(ones(1,numel(a_p)).*2,a_a,'k')
hold on; scatter(ones(1,numel(a_p)).*3,a_m,'k')
hold on; scatter(ones(1,numel(a_p)).*4,a_v,'k')
hold on; errorbar(1,mean(a_p),std(a_p)./sqrt(numel(a_p)),'ro')
hold on; errorbar(2,mean(a_a),std(a_a)./sqrt(numel(a_a)),'ro')
hold on; errorbar(3,mean(a_m),std(a_m)./sqrt(numel(a_m)),'ro')
hold on; errorbar(4,mean(a_v),std(a_v)./sqrt(numel(a_v)),'ro')

[~,p,~,stats] = ttest(angle_abs(tuned_units),phase_abs(tuned_units));
text(1,-8,['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])
[~,p,~,stats] = ttest(angle_abs(tuned_units),amp_abs(tuned_units));
text(1,-6,['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])
[~,p,~,stats] = ttest(angle_abs(tuned_units),midpoint_abs(tuned_units));
text(1,-4,['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])
[~,p,~,stats] = ttest(angle_abs(tuned_units),velocity_abs(tuned_units));
text(1,-2,['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])

set(gca,'ylim',[-10 20],'xlim',[0 5],'xtick',1:4,'xticklabel',{'phase','amp','midpoint','velocity'})
fn = 'diff_abs_mod.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% HSV PLOTTING WHISK
% [~,sort_idx] = sort(w_peak_idx_angle);
% 
% peaks = {w_peak_idx_angle,w_peak_idx_phase,w_peak_idx_amp,w_peak_idx_midpoint,w_peak_idx_velocity};
% mod_idx = {w_mod_idx_angle,w_mod_idx_phase,w_mod_idx_amp,w_mod_idx_midpoint,w_mod_idx_velocity};
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

%% example whisker trace (SF)
for k = 33
    currentCell = U{k};
    
    %      trialSample = [datasample(1:currentCell.k,5)];
    %      trialSample = [134 136 datasample(1:currentCell.k,1)]; %CELL 83
    %         trialSample = [99 136 datasample(1:currentCell.k,3)]; %% CELL 22
    trialSample = [119 98 43 datasample(1:currentCell.k,2)]; %% CELL 33 EXAMPLE 119 SO GOOD!CELL TUNED TO ANGLE 10 . Dam, it's also phase tuned. LULz
    %         trialSample = [60 1 datasample(1:currentCell.k,1)]; %% CELL 74 maybe 96 too
    figure(2310);clf
    for g = 1:length(trialSample)
        
        whisk = currentCell.S_ctk(1,:,trialSample(g));
        midpoint = currentCell.S_ctk(4,:,trialSample(g));
        amp = currentCell.S_ctk(3,:,trialSample(g));
        amplitude = whisk;
        amplitude(amp<5) = nan;
        touchesOn = [find(currentCell.S_ctk(9,:,trialSample(g))==1) find(currentCell.S_ctk(12,:,trialSample(g))==1)];
        touchesOff = [find(currentCell.S_ctk(10,:,trialSample(g))==1) find(currentCell.S_ctk(13,:,trialSample(g))==1)];
        short_touches = (touchesOff-touchesOn)<=5;
        fill_values = touchesOff(short_touches)+5;
        touchesOff(short_touches) = fill_values;
        
        
        if ~isempty(touchesOn)
            for p = 1:numel(touchesOn)
                amplitude(touchesOn(p)+1:touchesOff(p)+30) = nan;
            end
        end
        
        spikes = find(currentCell.R_ntk(1,:,trialSample(g))==1);
        subplot(length(trialSample),1,g)
        plot(1:currentCell.t,whisk,'k')
        hold on;plot(1:currentCell.t,amplitude,'c')
        
        hold on; scatter(spikes,whisk(spikes),'ko','filled','markerfacecolor',[.8 .8 .8])
        hold on; scatter(spikes,amplitude(spikes),'ko','filled')
        
        if ~isempty(touchesOn) && ~isempty(spikes)
            for f = 1:length(touchesOn)
                tp = touchesOn(f):touchesOff(f);
                hold on; plot(tp,whisk(tp),'r-');
            end
        end
        
        title(['cell num ' num2str(k) ' at trial num ' num2str(trialSample(g))])
        set(gca,'ylim',[-30 70])
    end

end
fn = 'example_traces .eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% firing rate X depth of recording (SF)
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
whisk_cells =  cellfun(@(x) x.is_tuned,whisk_struct.angle)==1;
jc_silent_cell = [766 819 895 631 776 815 910 871 844 902   941   840   888   748   732   940   686   944   950   933]; %Phils  from 902

figure(480);clf
subplot(3,1,[1 2])
scatter( cellfun(@(y) y.meta.depth,U),cellfun(@(y) mean(y.R_ntk(:))*1000, U),'k','filled')
hold on;scatter( cellfun(@(y) y.meta.depth,U(whisk_cells)),cellfun(@(y) mean(y.R_ntk(:))*1000, U(whisk_cells)),'c','filled')
hold on; scatter(jc_silent_cell,ones(length(jc_silent_cell),1)*-1,[],'filled','markerfacecolor',[.6 .6 .6]);
hold on; plot([700 700],[0 30],'--k')
hold on; plot([900 900],[0 30],'--k')
title(['n = response (' num2str(numel(U)) ') + silent (' num2str(numel(jc_silent_cell)) ')' ])
set(gca,'ytick',0:10:30,'xtick',600:100:1000,'xlim',[600 1000])
ylabel('firing rate (hz)');
xlabel('depth (um)')


hold on; subplot(3,1,3)
histogram(cellfun(@(y) y.meta.depth,U(setdiff(1:length(U),whisk_cells))),600:25:1000,'facecolor','k')
hold on; histogram(cellfun(@(y) y.meta.depth,U(whisk_cells)),600:25:1000,'facecolor','c')
hold on;histogram(jc_silent_cell,600:25:1000,'facecolor',[.6 .6 .6])
set(gca,'xtick',600:100:1000,'xlim',[600 1000])

fn = 'scatter_depth_firingrate.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% CDF of sampling (SF)
built = find(cellfun(@(x) isfield(x,'stim_response'),whisk_struct.angle));
bars_built = cellfun(@(x) isfield(x.stim_response,'bars_fit'),whisk_struct.angle(built));
full_build = built(bars_built);

masks = cellfun(@(x) maskBuilder(x),U(full_build),'uniformoutput',0);

thetas = cellfun(@(x,y) squeeze(x.S_ctk(1,:,:)) .* y.whisking,U(full_build),masks,'uniformoutput',0);

hist_thetas = cellfun(@(x) histcounts(x(:),-25.5:80.5),thetas,'uniformoutput',0);

cdf_sums = cell2mat(cellfun(@(x) cumsum(x),hist_thetas','uniformoutput',0))';
cdf_mat = cdf_sums./(cdf_sums(end,:));

mean_time_per_trial_whisking = cellfun(@(x) nansum(x.whisking(:))./size(x.whisking,2)./1000,masks); % in seconds
num_bins_per_cell = cellfun(@(x) numel(x.stim_response.raw_stim),whisk_struct.angle(full_build));

figure(23);clf
subplot(2,3,[1 2 4 5])
plot(repmat(-25:80,size(cdf_mat,2),1)',cdf_mat,'color',[.9 .9 .9])
hold on; shadedErrorBar(-25:80,nanmean(cdf_mat,2),nanstd(cdf_mat,[],2)./sqrt(numel(full_build)))
axis square
set(gca,'ytick',0:.2:1,'xtick',-20:20:80)
xlabel('whisker angle during whisking')
ylabel('cumulative proportion of whisks')

subplot(2,3,3);
scatter(ones(1,numel(full_build)),mean_time_per_trial_whisking,'k','markeredgecolor',[.9 .9 .9])
hold on; errorbar(1,mean(mean_time_per_trial_whisking),std(mean_time_per_trial_whisking),'ko')
set(gca,'ylim',[0 3],'ytick',0:1:3,'xtick',[])
ylabel('time spent whisking per trial (s)');
axis square
title(['mean=' num2str(mean(mean_time_per_trial_whisking)) ' , SD=' num2str(std(mean_time_per_trial_whisking))])

subplot(2,3,6);
scatter(ones(1,numel(full_build)),num_bins_per_cell,'k','markeredgecolor',[.9 .9 .9])
hold on; errorbar(1,mean(num_bins_per_cell),std(num_bins_per_cell),'ko')
set(gca,'ylim',[0 80],'ytick',0:25:100,'xtick',[])
ylabel('number of whisking bins');
axis square
title(['mean=' num2str(mean(num_bins_per_cell)) ' , SD=' num2str(std(num_bins_per_cell))])

fn = 'whisk_sampling.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% phase map (SF)
hilbert_feature = 'phase';

tuned_units = cellfun(@(x) x.is_tuned,whisk_struct.phase)==1;

preferred_phase = cellfun(@(x) x.calculations.tune_peak,whisk_struct.phase(tuned_units));
phase_mod = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.phase(tuned_units));
retractions = sign(preferred_phase)<0;
conversions = zeros(1,numel(preferred_phase));
conversions(retractions)=180;
deg_per_pie_unit =  180/pi;

preferred_degs = (abs(preferred_phase) .* deg_per_pie_unit) + conversions;

figure(380);clf
polarplot(preferred_phase,phase_mod,'o');

fn = 'phase_map.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% ASIDE: phase x angle heat 
tuned_units = find((double(cellfun(@(x) x.is_tuned,whisk_struct.angle)==1) + double(cellfun(@(x) x.is_tuned,whisk_struct.phase)==1))>0);
rc = numSubplots(numel(tuned_units));
bins = 20;
figure(520);clf
for b = 1:numel(tuned_units)
    array = U{tuned_units(b)};
    spikes = squeeze(array.R_ntk);
    
    %whisking mask
    timePostTouchToTrim = 30;
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touchOnIdx = touchOnIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchOffIdx = touchOffIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(touchOnIdx)
        touchEx_mask(touchOnIdx(i):touchOffIdx(i)+timePostTouchToTrim) = NaN; %added 30 to add time from touch offset
    end
    touchEx_mask(1:100,:) = 1; %since bleedover from end of trials before, tihs ensure we keep end
    amplitude = squeeze(array.S_ctk(3,:,:));
    whisking = nan(size(squeeze(array.S_ctk(1,:,:))));
    whisking(amplitude>5)=1;
    whisking_mask = whisking .* touchEx_mask;
    

    angle = squeeze(array.S_ctk(1,:,:));
    phase = squeeze(array.S_ctk(5,:,:));

   	angle_feature = angle(whisking_mask==1);
    phase_feature = phase(whisking_mask==1);
    
    filtered_spikes =spikes(whisking_mask==1);
    
    [phase_sorted] = binslin(phase_feature,[angle_feature filtered_spikes],'equalE',bins+1,-pi,pi);
    [angle_sorted,ang_sorted_by] = cellfun(@(x) binslin(x(:,1),x(:,2),'equalE',bins+1,prctile(angle_feature,1),prctile(angle_feature,99)),phase_sorted,'uniformoutput',0);
    
    angX = round(nanmedian(cell2mat(cellfun(@(x) cellfun(@(y) nanmedian(y),x),ang_sorted_by,'uniformoutput',0)'),2));
    countsX = cell2mat(cellfun(@(x) cellfun(@(y) sum(y),x),angle_sorted,'uniformoutput',0)');
    numSamplesX = cell2mat(cellfun(@(x) cellfun(@(y) numel(y),x),angle_sorted,'uniformoutput',0)');
    
   figure(520);
   subplot(rc(1),rc(2),b);
   plot_data = (countsX./numSamplesX).*1000;
   imagesc(plot_data)
   set(gca,'xtick',[1 (1+bins)/2 bins],'xticklabel',{'-\pi',0,'\pi'},'ydir','normal',...
       'ytick',1:3:bins,'yticklabel',angX(1:3:bins))
%    xlabel('phase free whisking')
%    ylabel('angle free whisking')
   axis square
   colormap turbo
   colorbar
   caxis([0 prctile(countsX(:)./numSamplesX(:).*1000,99)])
end
