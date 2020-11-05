function modulation_width(touch_struct)

tuned_units = find(cellfun(@(x) x.is_tuned==1,touch_struct.pole));

peak_response = cellfun(@(x) x.calculations.tune_peak,touch_struct.pole(tuned_units));
interp_norm_y = cell(1,length(peak_response));
for g = 1:length(peak_response)
    sr = touch_struct.pole{tuned_units(g)}.stim_response;
    centered_x = sr.values(:,1) - peak_response(g) ;
    norm_y = norm_new(sr.values(:,2));
    
    interp_centered_x = -2:.1:2;
    raw_x = round(centered_x,2);
    [~,idx] = unique(raw_x);
    interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
    figure(30);
    hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
end

pop_mean = nanmean(cell2mat(interp_norm_y'));
pop_sem = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
hold on; shadedErrorBar(interp_centered_x,pop_mean,pop_sem,'k')

set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
xlabel('distance from peak (mm)')
axis square

% fn = 'all_modulation_width.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
