function tuning_preference_scatter(touch_struct,whisk_struct,hilbert_var)

tUnits = cellfun(@(x) x.is_tuned==1,touch_struct.(hilbert_var));
wUnits = cellfun(@(x) x.is_tuned==1,whisk_struct.(hilbert_var));

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

touch_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],touch_struct.(hilbert_var)(tUnits),'uniformoutput',0)') ;
whisking_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.(hilbert_var)(wUnits),'uniformoutput',0)');

%scatter of whisking (Y) vs touch (X)
figure(3850);clf
if strcmp(hilbert_var,'pole')
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*2.5,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal')%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*2.5,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical')%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko')
    
    set(gca,'xlim',[-3 3],'ylim',[-3 3],'xdir','reverse','ydir','reverse',...
        'xtick',-1:1:1,'ytick',-1:1:1)
    hold on; plot([-1 1],[-1 1],'--k')
    
    figure(3851);clf
    subplot(2,1,1)
    histogram(touch_pw(:,1),-3:.20:3,'facecolor','b','facealpha',1)
    set(gca,'xdir','reverse','xlim',[-3 3])
    
    subplot(2,1,2);
    histogram(whisking_pw(:,1),-3:.20:3,'facecolor','c','facealpha',1)
    set(gca,'xdir','reverse','xlim',[-3 3],'ytick',0:2:6,'ylim',[0 6])
elseif strcmp(hilbert_var,'phase')
    subplot(4,1,[1:2])
%     hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-4,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
%     hold on; errorbar(ones(1,length(touch_nonIX_idx))*-4,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[-4 4],'ylim',[-4 4],...
        'xtick',-pi:pi:pi,'ytick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'yticklabel',{'-\pi',0,'\pi'})
    hold on; plot([-4 4],[-4 4],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(touch_nonIX_idx,1),-pi:pi/4:pi,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(whisk_nonIX_idx,1),-pi:pi/4:pi,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
elseif strcmp(hilbert_var,'angle')
    subplot(4,1,[1:2])
    %     hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-30,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    %     hold on; errorbar(ones(1,length(touch_nonIX_idx))*-30,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[-30 70],'ylim',[-30 70],...
        'xtick',-40:20:80,'ytick',-40:20:80)
    hold on; plot([-40 80],[-40 80],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(touch_nonIX_idx,1),-30:10:80,'facecolor','b','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-30 70],'ylim',[0 .5])
    subplot(4,1,4);
    histogram(whisking_pw(whisk_nonIX_idx,1),-430:10:80,'facecolor','c','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-30 70],'ylim',[0 .5])
    
elseif strcmp(hilbert_var,'amplitude')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*0,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*0,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[0 40],'ylim',[0 40],...
        'xtick',0:10:40,'ytick',0:10:40)
    hold on; plot([0 40],[0 40],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),0:5:40,'facecolor','b','facealpha',1,'normalization','probability')
    set(gca,'xlim',[0 40],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1,'normalization','probability')
    set(gca,'xlim',[0 40],'ylim',[0 20])
elseif strcmp(hilbert_var,'midpoint')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-20,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*-20,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[-20 60],'ylim',[-20 60],...
        'xtick',-20:20:60,'ytick',-20:20:60)
    hold on; plot([-20 60],[-20 60],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),-20:10:60,'facecolor','b','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-20 60],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-20 60],'ylim',[0 20])
end

figure(3850);subplot(4,1,[1 2])
axis square
xlabel('whisk tune peak');ylabel('touch tune peak')
title(['whisk=' num2str(numel(whisk_nonIX_idx)) ', touch=' num2str(numel(touch_nonIX_idx)) ', both=' num2str(numel(touch_ix_idx))])

% figure(3850);
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
% fn = [hilbert_var '_intersect.eps'];
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])


% Distance from unity line (D stats)  
t_w_vect = [touch_pw(touch_ix_idx,1) whisking_pw(whisk_ix_idx,1)];
if strcmp(hilbert_var,'angle')
    diagonal_vect = repmat([-30:80]',1,2);
elseif strcmp(hilbert_var,'phase')
    diagonal_vect = repmat((-pi:pi/20:pi)',1,2);
end

euclid_dist_all = pdist2(t_w_vect,diagonal_vect);
min_ed = min(euclid_dist_all,[],2);
[mean(min_ed) std(min_ed)]
[~,p,~,stats] = ttest(min_ed)