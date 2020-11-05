function modulation_width_comparison(U,touch_struct)
tuned_units = cellfun(@(x) x.is_tuned==1,touch_struct.pole);
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
trained = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

naive_tuned_units = intersect(find(tuned_units),find(naive));
expert_tuned_units = intersect(find(tuned_units),find(trained));

tuned_units_list = {naive_tuned_units,expert_tuned_units};

figure(30);clf
pop_mean = cell(1,2);
pop_sem = cell(1,2);
for b=1:length(tuned_units_list)
    peak_response = cellfun(@(x) x.calculations.tune_peak,touch_struct.pole(tuned_units_list{b}));
    for g = 1:length(peak_response)
        sr = touch_struct.pole{tuned_units_list{b}(g)}.stim_response;
        centered_x = sr.values(:,1) - peak_response(g) ;
        norm_y = normalize_var(sr.values(:,2),0,1);
        
        interp_centered_x = -2:.1:2;
        raw_x = round(centered_x,2);
        [~,idx] = unique(raw_x);
        interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
        figure(30);subplot(1,3,b)
        hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
    end
    
    pop_mean{b} = nanmean(cell2mat(interp_norm_y'));
    pop_sem{b} = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
    hold on; shadedErrorBar(interp_centered_x,pop_mean{b},pop_sem{b},'k')
    
    set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
    xlabel('distance from peak (mm)')
    axis square
end

subplot(1,3,3);
shadedErrorBar(interp_centered_x,pop_mean{1},pop_sem{1},'k')
hold on;  shadedErrorBar(interp_centered_x,pop_mean{2},pop_sem{2},'b')
set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
xlabel('distance from peak (mm)')
axis square

% figure(30);
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
% fn = 'naive_expert_modulation_width.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
