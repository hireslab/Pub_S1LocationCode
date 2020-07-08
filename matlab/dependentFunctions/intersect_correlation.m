function intersect_correlation(tStruct,wStruct,hilbertVar)
% Build heatmap of equal tuning sized tuning
touch_units = cellfun(@(x) x.is_tuned==1,tStruct);
whisk_units = cellfun(@(x) x.is_tuned==1,wStruct);
whisk_or_touch = find(sum(double(touch_units)+double(whisk_units),1)>0);
whisk_and_touch = intersect(find(touch_units),find(whisk_units));

[~,single_tuned] = setdiff(whisk_or_touch,whisk_and_touch);
[~,both_tuned] = setdiff(1:length(whisk_or_touch),single_tuned);

touch_heat = cell(1,numel(whisk_or_touch));
whisk_heat = cell(1,numel(whisk_or_touch));

for g = 1:numel(whisk_or_touch)
    try
        curr_t = tStruct{whisk_or_touch(g)}.stim_response.values;
        curr_w = wStruct{whisk_or_touch(g)}.stim_response.values;
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

%% curve matching correlation
titles = {'single','both'};
which_tuned = {single_tuned,both_tuned};

for b = 1:numel(which_tuned)
    ix_touch = touch_heat(which_tuned{b});
    ix_whisk = whisk_heat(which_tuned{b});
    
    to_toss = unique([find(cellfun(@isempty,ix_touch)) find(cellfun(@isempty,ix_whisk))]);
    ix_touch(to_toss) = [];
    ix_whisk(to_toss) = [];
    
    %corr for matching bits
    matching_bins = cellfun(@(x,y) intersect(find(~isnan(x)),find(~isnan(y))),ix_touch,ix_whisk,'uniformoutput',0);
    match_stretch = mean(cellfun(@numel,matching_bins));
    matched_touch = cellfun(@(x,y) x(y),ix_touch,matching_bins,'uniformoutput',0);
    matched_whisk = cellfun(@(x,y) x(y),ix_whisk,matching_bins,'uniformoutput',0);
    
    
    stretch_bins = median(cellfun(@(x) sum(~isnan(x)),[ix_touch ix_whisk]));
    raw_touch = cellfun(@(x) x(~isnan(x)),ix_touch,'uniformoutput',0);
    new_touch = cellfun(@(x) interp1(linspace(1,stretch_bins,numel(x)),x,1:stretch_bins),raw_touch,'uniformoutput',0);
    raw_whisk = cellfun(@(x) x(~isnan(x)),ix_whisk,'uniformoutput',0);
    new_whisk = cellfun(@(x) interp1(linspace(1,stretch_bins,numel(x)),x,1:stretch_bins),raw_whisk,'uniformoutput',0);
    match_stretch_touch=cellfun(@(x) interp1(linspace(1,match_stretch,numel(x)),x,1:match_stretch),matched_touch,'uniformoutput',0);
    match_stretch_whisk=cellfun(@(x) interp1(linspace(1,match_stretch,numel(x)),x,1:match_stretch),matched_whisk,'uniformoutput',0);
    
    
    match_corr = cellfun(@(x,y,z) corr(x(z)',y(z)'),ix_touch,ix_whisk,matching_bins);
    
    shuff_times = 1000;
    for f = 1:shuff_times
        shuff_touch = new_touch(randperm(numel(match_stretch_touch)));
        shuff_whisk = new_whisk(randperm(numel(match_stretch_whisk)));
        shuff_corr{f} = cellfun(@(x,y) corr(x',y'),shuff_touch,shuff_whisk);
    end
    [~,p] = kstest2(match_corr,cell2mat(shuff_corr));
    
    figure(53+b);clf
    subplot(1,2,1)
    histogram(match_corr,-1:.2:1,'facecolor','k')
    xlabel('correlation coefficient')
    title('intersecting locations')
    axis square
    subplot(1,2,2)
    histogram(cell2mat(shuff_corr),-1:.2:1,'facecolor','b','normalization','probability')
    hold on; histogram(match_corr,-1:.2:1,'facecolor','k','normalization','probability')
    legend('shuff','data')
    suptitle(['KS p = ' num2str(p)])
    axis square
    
    title(titles{b});
end
