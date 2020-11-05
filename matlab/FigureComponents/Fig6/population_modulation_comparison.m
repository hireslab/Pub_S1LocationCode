function population_modulation_comparison(touch_struct,whisk_struct,hilbert_var)

tUnits = cellfun(@(x) x.is_tuned==1,touch_struct.(hilbert_var));
wUnits = cellfun(@(x) x.is_tuned==1,whisk_struct.(hilbert_var));

co_tuned = find(tUnits.*wUnits);
touch_only = setdiff(find(tUnits),co_tuned);
whisk_only = setdiff(find(wUnits),co_tuned);
untuned = setdiff(1:length(touch_struct.(hilbert_var)),[co_tuned, touch_only, whisk_only]);

touch_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,touch_struct.(hilbert_var));
whisk_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.(hilbert_var));

touch_mod_abs_tune = cellfun(@(x) x.calculations.mod_idx_abs,touch_struct.(hilbert_var)(co_tuned));
whisk_mod_abs_tune = cellfun(@(x) x.calculations.mod_idx_abs,whisk_struct.(hilbert_var)(co_tuned));

figure(8540);clf
hold on;scatter(whisk_mod_abs(untuned),touch_mod_abs(untuned),'ko')
hold on;scatter(whisk_mod_abs(touch_only),touch_mod_abs(touch_only),'filled','ro')
hold on;scatter(whisk_mod_abs(whisk_only),touch_mod_abs(whisk_only),'filled','co')
hold on;scatter(whisk_mod_abs_tune,touch_mod_abs_tune,'filled','ko')
hold on; errorbar(mean(whisk_mod_abs(untuned)),mean(touch_mod_abs(untuned)),...
    std(touch_mod_abs(untuned))./sqrt(numel((untuned))),std(touch_mod_abs(untuned))./sqrt(numel((touch_only))),...
    std(whisk_mod_abs(untuned))./sqrt(numel((untuned))),std(whisk_mod_abs(untuned))./sqrt(numel((touch_only)))...
    ,'go','capsize',0)
hold on; errorbar(mean(whisk_mod_abs_tune),mean(touch_mod_abs_tune),...
    std(touch_mod_abs_tune)./sqrt(numel(co_tuned)),std(touch_mod_abs_tune)./sqrt(numel(co_tuned)),...
    std(whisk_mod_abs_tune)./sqrt(numel(co_tuned)),std(whisk_mod_abs_tune)./sqrt(numel(co_tuned))...
    ,'ko','capsize',0)
hold on; errorbar(mean(whisk_mod_abs(touch_only)),mean(touch_mod_abs(touch_only)),...
    std(touch_mod_abs(touch_only))./sqrt(numel((touch_only))),std(touch_mod_abs(touch_only))./sqrt(numel((touch_only))),...
    std(whisk_mod_abs(touch_only))./sqrt(numel((touch_only))),std(whisk_mod_abs(touch_only))./sqrt(numel((touch_only)))...
    ,'ro','capsize',0)
hold on; errorbar(mean(whisk_mod_abs(whisk_only)),mean(touch_mod_abs(whisk_only)),...
    std(touch_mod_abs(whisk_only))./sqrt(numel((whisk_only))),std(touch_mod_abs(whisk_only))./sqrt(numel((whisk_only))),...
    std(whisk_mod_abs(whisk_only))./sqrt(numel((whisk_only))),std(whisk_mod_abs(whisk_only))./sqrt(numel((whisk_only)))...
    ,'co','capsize',0)


hold on; plot([.1 100],[.1 100],'--k')
set(gca,'xlim',[.1 100],'ylim',[.1 100],...
'yscale','log','xscale','log')
axis square
xlabel('whisk abs mod (spks/s)')
ylabel('touch abs mod (spks/s)')

[~,p,~,stats] = ttest(whisk_mod_abs,touch_mod_abs);
text(5, 5,['all built : ' 'p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])
[~,p,~,stats] = ttest(whisk_mod_abs_tune,touch_mod_abs_tune);
text(5, 4,['co tuned : ' 'p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
% fn = [hilbert_var '_touch_x_whisk_absmod.eps'];
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])