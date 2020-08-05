function PlotPhaseAngleModulation(whisk_struct)
%% part 1: phase and angle curves
tuned_units = (double(cellfun(@(x) x.is_tuned,whisk_struct.angle)==1) + double(cellfun(@(x) x.is_tuned,whisk_struct.phase)==1))>0;
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
suptitle('co-tuned curves modulation')

%% part 2: plot modulation unity
tuned_units = (double(cellfun(@(x) x.is_tuned,whisk_struct.angle)==1) + double(cellfun(@(x) x.is_tuned,whisk_struct.phase)==1))>0;

% calculate angle modulation 
phase_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(tuned_units));
angle_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(tuned_units));

% inherit cotuned from Fig6
% ix = logical(tuned_units .* cotuned);
% nix = logical(tuned_units .* ~cotuned);
% 
% nix_phase = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(nix));
% nix_angle = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(nix)); 
% 
% ix_phase = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.phase(ix));
% ix_angle = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.angle(ix));

% plotting absolute modulation 
figure(57);clf
scatter(phase_mod_abs,angle_mod_abs,'k')
% hold on; scatter(ix_phase, ix_angle,'filled','r')
% hold on; scatter(nix_phase, nix_angle, 'filled','b')
hold on; plot([0 30], [0 30],'k--')
hold on; errorbar(mean(phase_mod_abs),mean(angle_mod_abs),...
    std(angle_mod_abs)./sqrt(sum(tuned_units)),std(angle_mod_abs)./sqrt(sum(tuned_units)),...
    std(phase_mod_abs)./sqrt(sum(tuned_units)),std(phase_mod_abs)./sqrt(sum(tuned_units))...
    ,'ko','capsize',0)
% hold on; errorbar(mean(ix_phase),mean(ix_angle),...
%     std(ix_angle)./sqrt(sum(tuned_units)),std(ix_angle)./sqrt(sum(tuned_units)),...
%     std(ix_phase)./sqrt(sum(tuned_units)),std(ix_phase)./sqrt(sum(tuned_units))...
%     ,'ro','capsize',0)
% hold on; errorbar(mean(nix_phase),mean(nix_angle),...
%     std(nix_angle)./sqrt(sum(tuned_units)),std(nix_angle)./sqrt(sum(tuned_units)),...
%     std(nix_phase)./sqrt(sum(tuned_units)),std(nix_phase)./sqrt(sum(tuned_units))...
%     ,'bo','capsize',0)
set(gca,'xlim',[0 30],'ylim',[0 30])
axis square
xlabel('Phase absolute modulation')
ylabel('Angle absolute modulation')
[~,p,~,stats] = ttest(phase_mod_abs,angle_mod_abs);
title(['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])
