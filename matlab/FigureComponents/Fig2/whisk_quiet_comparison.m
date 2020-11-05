function whisk_quiet_comparison(U)

% block out certain periods of whisk and quiet 
masks = cellfun(@(x) maskBuilder(x),U,'uniformoutput',0);
whisking_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.whisking .*y.touch,U,masks,'uniformoutput',0);
quiet_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.quiet .*y.touch,U,masks,'uniformoutput',0);
whisking_tp = cellfun(@(x,y) nansum(nansum(y.whisking .*y.touch)),U,masks);
quiet_tp = cellfun(@(x,y) nansum(nansum(y.quiet .*y.touch)),U,masks);

for b = 1:numel(whisking_spks_mat)
    w = whisking_spks_mat{b}(:);
    q = quiet_spks_mat{b}(:);
    n1 = nansum(w==1); N1 = nansum(w==0);
    n2 = nansum(q==1); N2 = nansum(q==0);

    %generate sample data using proportions in data since matlab chi square
    %requires vectors to be equal lengths.
    x1 = [repmat('a',N1,1); repmat('b',N2,1)];
    x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
    
    %chi squared calculation
    [~,~,p(b)] = crosstab(x1,x2); %

end

fr_whisk = (cellfun(@(x) nansum(x(:)),whisking_spks_mat)./whisking_tp)*1000;
fr_quiet = (cellfun(@(x) nansum(x(:)),quiet_spks_mat)./quiet_tp)*1000;

red_dots = intersect(find(p<.05),find(fr_whisk>fr_quiet));
blue_dots = intersect(find(p<.05),find(fr_whisk<fr_quiet));
gray_dots = setdiff(1:numel(masks),[red_dots blue_dots]);

figure(480);clf
scatter(fr_quiet(p<.05),fr_whisk(p<.05),'filled','k')
hold on; scatter(fr_quiet(gray_dots),fr_whisk(gray_dots),'ko')
axis square
hold on; plot([0.01 100],[0.01 100],'--k')
set(gca,'xlim',[0 100],'ylim',[0 100],'xscale','log','yscale','log')
xlabel('quiet FR');ylabel('whisking FR')

fn = '2A.eps';
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% PLos Reviewer Edits (1.5)
low_fr = fr_whisk(fr_quiet< 5) - fr_quiet(fr_quiet<5);
high_fr = fr_whisk(fr_quiet>= 5) - fr_quiet(fr_quiet>=5);

[~,low_p] = ttest(fr_whisk(fr_quiet< 5),fr_quiet(fr_quiet<5));
[~,high_p] = ttest(fr_whisk(fr_quiet>= 5) - fr_quiet(fr_quiet>=5));