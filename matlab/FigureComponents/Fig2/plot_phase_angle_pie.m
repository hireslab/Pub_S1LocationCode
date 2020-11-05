function plot_phase_angle_pie(whisk_struct)

whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.angle);
whisk_phase_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.phase);
co_tuned = intersect(find(whisk_angle_tuned),find(whisk_phase_tuned)); 
angle_only = setdiff(find(whisk_angle_tuned),co_tuned);
phase_only = setdiff(find(whisk_phase_tuned),co_tuned);
untuned = setdiff(1:length(whisk_struct.angle),[co_tuned angle_only phase_only]);

figure(580620);clf
pie([numel(untuned) numel(co_tuned) numel(angle_only) numel(phase_only)],...
    {'not tuned','co-tuned','angle only','phase only'})

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'proportion_angle_phase_pie.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
