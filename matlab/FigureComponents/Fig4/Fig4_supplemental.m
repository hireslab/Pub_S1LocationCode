function Fig4_supplemental(U)

% num touches
touch_nums = cellfun(@(x) sum(squeeze(x.S_ctk(9,:,:)==1) + squeeze(x.S_ctk(13,:,:)==1)),U,'uniformoutput',0);
whisk_mat = cellfun(@(x) squeeze(x.S_ctk(3,:,:)>5),U,'uniformoutput',0);

w_prop_expert = cellfun(@(x) mean(x(:)),whisk_mat(trained));
w_prop_naive = cellfun(@(x) mean(x(:)),whisk_mat(naive));

expert_mean = [cellfun(@(x) nanmean(x),touch_nums(trained))];
naive_mean = [cellfun(@(x) nanmean(x),touch_nums(naive))];

figure(48);clf
subplot(1,2,1)
histogram(expert_mean,0:1:20,'normalization','probability','facealpha',1,'facecolor','y')
hold on; histogram(naive_mean,0:1:20,'normalization','probability','facecolor',[.8 .8 .8])
set(gca,'xtick',0:5:20,'xlim',[0 20],'ylim',[0 .25],'ytick',0:.25:1)
xlabel('touch counts')
ylabel('proportion of units')
[~,ksp_counts] = kstest2(expert_mean,naive_mean);
disp(['touch distribution KS p = ' num2str(ksp_counts)])
title(num2str(ksp_counts));

figure(48);subplot(1,2,2)
histogram(w_prop_expert,0:.05:1,'normalization','probability','facealpha',1,'facecolor','y')
hold on; histogram(w_prop_naive,0:.05:1,'normalization','probability','facecolor',[.8 .8 .8])
set(gca,'xtick',0:.5:1,'xlim',[0 1],'ylim',[0 .5],'ytick',0:.25:1)
xlabel('proportion of time whisking')
ylabel('proportion of units')
[~,ksp_whisk] = kstest2(w_prop_expert,w_prop_naive);
disp(['touch distribution KS p = ' num2str(ksp_whisk)])
title(num2str(p_whisk));

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
% fn = 'naive_expert_stimulus.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])