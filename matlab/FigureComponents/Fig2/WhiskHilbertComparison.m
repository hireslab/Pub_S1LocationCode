function WhiskHilbertComparison(U, whisk_struct)

selectedCells = 1:length(whisk_struct.angle);

% calculate relative modulation index for each cell for heatmap
w_mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.angle(selectedCells));
w_mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.phase(selectedCells));
w_mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.amplitude(selectedCells));
w_mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.midpoint(selectedCells));
w_mod_idx_velocity = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.velocity(selectedCells));

% calculate absolute modulation index for tuned units in comparing
% modulation
angle_abs= cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(selectedCells));
phase_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(selectedCells));
amp_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.amplitude(selectedCells));
midpoint_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.midpoint(selectedCells));
velocity_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.velocity(selectedCells));

tuned_units = find(cellfun(@(x) x.is_tuned,whisk_struct.angle));
a_p = angle_abs(tuned_units) - phase_abs(tuned_units);
a_a = angle_abs(tuned_units) - amp_abs(tuned_units);
a_m = angle_abs(tuned_units) - midpoint_abs(tuned_units);
a_v = angle_abs(tuned_units) - velocity_abs(tuned_units);

% plot gray heatmap
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

% fn = 'gray_whisk_sort.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

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

% fn = 'modulation_correlation.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

% plot relative modulation to angle
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

% fn = 'diff_abs_mod.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])