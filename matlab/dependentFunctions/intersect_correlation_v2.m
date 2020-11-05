function intersect_correlation_v2(tStruct,wStruct,hilbertVar)
% Build heatmap of equal tuning sized tuning
touch_units = cellfun(@(x) x.is_tuned==1,tStruct.(hilbertVar));
whisk_units = cellfun(@(x) x.is_tuned==1,wStruct.(hilbertVar));
whisk_or_touch = find(sum(double(touch_units)+double(whisk_units),1)>0);

touch_heat = cell(1,numel(whisk_or_touch));
whisk_heat = cell(1,numel(whisk_or_touch));

for g = 1:numel(whisk_or_touch)
    try
        curr_t = tStruct.(hilbertVar){whisk_or_touch(g)}.stim_response.values;
        curr_w = wStruct.(hilbertVar){whisk_or_touch(g)}.stim_response.values;
        %clean nan rows
        curr_w = curr_w(~any(isnan(curr_w),2),:);
        curr_t = curr_t(~any(isnan(curr_t),2),:);
        
        if strcmp(hilbertVar,'pole')
            touch_x = -1:.1:1;
            whisk_x = -2:.1:2;
        elseif strcmp(hilbertVar,'phase')
            touch_x = linspace(-pi,pi,21);
            whisk_x = linspace(-pi,pi,21);
        elseif strcmp(hilbertVar,'angle')
            touch_x = -30:80;
            whisk_x = -30:80;
        else
            touch_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
            whisk_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
        end
        touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
        whisk_heat{g} = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
    end
end

% curve matching correlation
%clean out cells that are not either whisk or touch tuned
to_toss = unique([find(cellfun(@isempty,touch_heat)) find(cellfun(@isempty,whisk_heat))]);
touch_heat(to_toss) = [];
whisk_heat(to_toss) = [];

%1) CORRELATION VIA BIN MATCHING
matching_bins = cellfun(@(x,y) intersect(find(~isnan(x)),find(~isnan(y))),touch_heat,whisk_heat,'uniformoutput',0);
matched_touch = cellfun(@(x,y) x(y),touch_heat,matching_bins,'uniformoutput',0);
matched_whisk = cellfun(@(x,y) x(y),whisk_heat,matching_bins,'uniformoutput',0);
%NEW ADDITION FOR ANDREW, normalize values from 0:1 200227
matched_touch = cellfun(@(x) normalize_var(x,0,1),matched_touch,'uniformoutput',0);
matched_whisk = cellfun(@(x) normalize_var(x,0,1),matched_whisk,'uniformoutput',0);

%2) CORRELATION VIA STRETCHING OR USE FOR SHUFFLE
match_stretch = mean(cellfun(@numel,matching_bins));
match_stretch_touch=cellfun(@(x) interp1(linspace(1,match_stretch,numel(x)),x,1:match_stretch),matched_touch,'uniformoutput',0);
match_stretch_whisk=cellfun(@(x) interp1(linspace(1,match_stretch,numel(x)),x,1:match_stretch),matched_whisk,'uniformoutput',0);

%real data correlation
real_corr = cellfun(@(x,y) corr(x',y'),matched_touch,matched_whisk);

%shuffle data correlation
num_shuffles = 1000;
shuff_corr = cell(1,1000);
for k = 1:num_shuffles
    shuff_whisk = match_stretch_whisk(randperm(length(match_stretch_whisk)));
    shuff_touch = match_stretch_touch(randperm(length(match_stretch_touch)));
    shuff_corr{k} =  cellfun(@(x,y) corr(x',y'),shuff_touch,shuff_whisk);
end
[~,p] = kstest2(real_corr,cell2mat(shuff_corr));


figure(53);clf
histogram(cell2mat(shuff_corr),-1:.2:1,'facecolor','b','normalization','probability')
hold on; histogram(real_corr,-1:.2:1,'facecolor','k','normalization','probability')
legend('shuff','data')
suptitle(['KS p = ' num2str(p)])
axis square

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
% fn = [hilbert_var '_correlations.eps'];
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

% correlation between tuning curves with whisking preferences. submission 1
% plos biology edit. 200609
[~,tid] = intersect(whisk_or_touch,find(touch_units));
[~,wid] = intersect(whisk_or_touch,find(whisk_units));
tw_ix = intersect(tid,wid); % find whisk and touch units

whisk_pref = cellfun(@(x) x.calculations.tune_peak,wStruct.(hilbertVar)(whisk_or_touch));
whisk_pref_tw = whisk_pref(tw_ix);
real_corr_tw = real_corr(tw_ix);
mdl = fitlm(whisk_pref_tw,real_corr_tw);

figure(3940);clf
scatter(whisk_pref_tw,real_corr_tw,'ko')
hold on; plot(-40:80,mdl.predict((-40:80)'),'k')
hold on; plot([-30 70],[0 0],'--k')
set(gca,'ytick',-1:.5:1,'xtick',-40:20:80,'xlim',[-30 70],'ylim',[-1.1 1.1])
axis square
xlabel('whisk angle preference')
ylabel('whisk and touch tuning curve correlation')
title(['corr = ' num2str(corr(whisk_pref_tw',real_corr_tw')) ', Rsq = ' num2str(mdl.Rsquared.Ordinary) ', n = ' num2str(length(real_corr_tw))])






