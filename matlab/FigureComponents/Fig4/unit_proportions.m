function unit_proportions(U,touch_struct)
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
trained = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

num_naive_touch = sum(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U(naive)));
num_expert_touch = sum(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U(trained)));
num_naive_OL = sum((cellfun(@(x) x.is_tuned==1,touch_struct.pole(naive)))) ;
num_expert_OL = sum((cellfun(@(x) x.is_tuned==1,touch_struct.pole(trained)))) ;

naive_proportion_touch = num_naive_touch ./ sum(naive);
expert_proportion_touch = num_expert_touch ./ sum(trained);
naive_proportion_OL  = num_naive_OL ./ sum(naive);
expert_proportion_OL  = num_expert_OL ./ sum(trained);

x = table([num_naive_touch - num_naive_OL; num_naive_OL],[ num_expert_touch - num_expert_OL ;num_expert_OL]);
[h,p,stats] = fishertest(x);

figure(99);clf
subplot(1,2,1)
hold on; bar(1:2,[naive_proportion_touch expert_proportion_touch],'facecolor',[.8 .8 .8])
hold on; bar(1:2,[naive_proportion_OL expert_proportion_OL],'facecolor','k')
ylabel('proportion of units')
set(gca,'xtick',1:2,'xticklabel',{['naive n=' num2str(sum(naive))],['expert n=' num2str(sum(trained))]},'ytick',0:.25:1,'ylim',[0 1])
legend('touch','location tuned','location','northwest')

subplot(1,2,2)
hold on; bar(1:2,[num_naive_OL./num_naive_touch, num_expert_OL./num_expert_touch],'facecolor','k')

ylabel('proportion of units')
set(gca,'xtick',1:2,'xticklabel',{['naive n=' num2str(num_naive_OL)],['expert n=' num2str(num_expert_OL)]},'ytick',0:.25:1,'ylim',[0 1])
title('proportion of touch that is OL')

% figure(99)
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
% fn = 'unit_distribution_bar.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])