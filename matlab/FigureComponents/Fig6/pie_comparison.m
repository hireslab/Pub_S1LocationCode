function pie_comparison(touch_struct,whisk_struct,hilbert_var)

tUnits = cellfun(@(x) x.is_tuned==1,touch_struct.(hilbert_var));
wUnits = cellfun(@(x) x.is_tuned==1,whisk_struct.(hilbert_var));

co_tuned = numel(find(tUnits.*wUnits));
touch_only = numel(setdiff(find(tUnits),find(tUnits.*wUnits)));
whisk_only = numel(setdiff(find(wUnits),find(tUnits.*wUnits)));
untuned = length(touch_struct.(hilbert_var)) - co_tuned - touch_only - whisk_only; 

% bootstrapping
samples = 1000; 
whisk_tuning = find(binornd(1,mean(wUnits),1,samples)==1);
touch_tuning = find(binornd(1,mean(tUnits),1,samples)==1);

co_tuned_bs = intersect(whisk_tuning,touch_tuning);
touch_only_bs = numel(setdiff(touch_tuning,co_tuned_bs));
whisk_only_bs = numel(setdiff(whisk_tuning,co_tuned_bs));
untuned_bs = samples - numel(co_tuned_bs) - touch_only_bs - whisk_only_bs;

% plotting
figure(340);clf
subplot(1,2,1)
pie([touch_only co_tuned  whisk_only untuned]);
subplot(1,2,2)
pie([touch_only_bs numel(co_tuned_bs)  whisk_only_bs untuned_bs]);
title ('bootstrap')
legend({'touch','both','whisk','none'},'location','southeast') 

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
% fn = ['B_intersect_proportions.eps'];
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
