function session_accuracy(U)

trained = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
task_list = {naive,trained};

figure(41);clf
for k = 1:length(task_list)
    accuracy =  cellfun(@(x) mean(x.meta.trialCorrect),U(task_list{k}));
    hold on; scatter(ones(1,sum(task_list{k})).*k,accuracy,'markeredgecolor',[.7 .7 .7]);
    hold on; errorbar(1*k,mean(accuracy),std(accuracy),'ko')
end
hold on; plot([0 3],[.5 .5],'--k')
set(gca,'ylim',[0 1],'ytick',0:.25:1,'xtick',[1,2],...
    'xticklabel',{'naive','expert'},'xlim',[0 3])
% 
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
% fn = 'naive_expert_performance.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])