function population_heatmap_builder(tStruct,wStruct,hilbertVar)
% this function plots heatmaps of whisking and touch tuning
%
% requires products of object_location_quantification and
% whisking_location_quantification first before running.
%
% uses touch tuning struct and whisk tuning struct for POLE position tuning
% to build population heatmaps

%% Plot heatmap of tuning across all units
tuned_units = cellfun(@(x) x.is_tuned==1,tStruct);
whisk_units = cellfun(@(x) x.is_tuned==1,wStruct);

sel_tstructs = tStruct(tuned_units);
sel_wstructs = wStruct(whisk_units);

touch_heat = cell(1,sum(tuned_units));
whisk_heat = cell(1,sum(whisk_units));

touch_ix_tmp = intersect(find(tuned_units),find(whisk_units));
whisk_ix_tmp = intersect(find(whisk_units),find(tuned_units));

[~,touch_ix_idx] = intersect(find(tuned_units),touch_ix_tmp);
[~,whisk_ix_idx] = intersect(find(whisk_units),whisk_ix_tmp);

%building heatmap for object location tuned touch units
for g = 1:numel(sel_tstructs)
    curr_t = sel_tstructs{g}.stim_response.values;
    
    %clean nan rows
    curr_t = curr_t(~any(isnan(curr_t),2),:);
    
    if strcmp(hilbertVar,'pole')
        touch_x = -1:.1:1;
    elseif strcmp(hilbertVar,'phase')
        touch_x = linspace(-pi,pi,21);
    elseif strcmp(hilbertVar,'angle')
        touch_x = -30:80;
    else
        touch_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
    end
    touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
    
end

%building heatmap for object location tuned whisk units
for d = 1:numel(sel_wstructs)
    curr_w = sel_wstructs{d}.stim_response.values;
    
    %clean nan rows
    curr_w = curr_w(~any(isnan(curr_w),2),:);
    
    if strcmp(hilbertVar,'pole')
        whisk_x = -2:.1:2;
    elseif strcmp(hilbertVar,'phase')
        whisk_x = linspace(-pi,pi,21);
    elseif strcmp(hilbertVar,'angle')
        whisk_x = -30:80;
    else
        whisk_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
    end
    whisk_heat{d} = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
end
%% heatmap of all units
figure(50);clf
title_names = {'all touch sorted','all whisk','both touch, touch sorted',...
    'both whisk, whisk sorted', 'both touch, whisk sorted', 'both whisk, touch sorted'};

unsorted_heat = norm_new(cell2mat(touch_heat')');
[~,t_max_idx] = max(unsorted_heat,[],1);
[~,t_idx] = sort(t_max_idx);
data{1} = unsorted_heat(:,t_idx)';

unsorted__whisk_heat = norm_new(cell2mat(whisk_heat')');
[~,w_max_idx] = max(unsorted__whisk_heat,[],1);
[~,w_idx] = sort(w_max_idx);
data{2} = unsorted__whisk_heat(:,w_idx)';

unsorted_heat_touch = norm_new(cell2mat(touch_heat(touch_ix_idx)')');
[~,t_max_idx] = max(unsorted_heat_touch,[],1);
[~,t_idx_ix] = sort(t_max_idx);
data{3} = unsorted_heat_touch(:,t_idx_ix)';

unsorted_whisk_heat = norm_new(cell2mat(whisk_heat(whisk_ix_idx)')');
[~,w_max_idx] = max(unsorted_whisk_heat,[],1);
[~,w_idx_ix] = sort(w_max_idx);
data{4} = unsorted_whisk_heat(:,w_idx_ix)';

unsorted_heat_touch = norm_new(cell2mat(touch_heat(touch_ix_idx)')');
data{5} = unsorted_heat_touch(:,w_idx_ix)';

unsorted_whisk_heat = norm_new(cell2mat(whisk_heat(whisk_ix_idx)')');
data{6} = unsorted_whisk_heat(:,t_idx_ix)';

for p = 1:numel(data)
    subplot(3,2,p)
    [nr,nc] = size(data{p});
    pcolor([data{p} nan(nr,1); nan(1,nc+1)]);
    caxis([0 1])
    shading flat;
    colormap turbo
    if strcmp(hilbertVar,'pole')
        set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
    elseif strcmp(hilbertVar,'phase')
        set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
    elseif strcmp(hilbertVar,'angle')
        set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
    end
    
    title(title_names{p})
end

%% squished heatmap
% touch_stretch_bins = median(cellfun(@(x) sum(~isnan(x)),touch_heat));
% whisk_stretch_bins = median(cellfun(@(x) sum(~isnan(x)),whisk_heat));
touch_stretch_bins = 50;
whisk_stretch_bins = 50;

raw_curves = cellfun(@(x) x(~isnan(x)),touch_heat,'uniformoutput',0);
new_curves = cellfun(@(x) interp1(linspace(1,touch_stretch_bins,numel(x)),x,1:touch_stretch_bins),raw_curves,'uniformoutput',0);
raw_whisk = cellfun(@(x) x(~isnan(x)),whisk_heat,'uniformoutput',0);
new_whisk = cellfun(@(x) interp1(linspace(1,whisk_stretch_bins,numel(x)),x,1:whisk_stretch_bins),raw_whisk,'uniformoutput',0);

figure(51);clf

unsorted_heat = norm_new(cell2mat(new_curves')');
[~,t_max_idx] = max(unsorted_heat,[],1);
[~,t_idx] = sort(t_max_idx);
data_squish{1} = unsorted_heat(:,t_idx)';

unsorted__whisk_heat = norm_new(cell2mat(new_whisk')');
[~,w_max_idx] = max(unsorted__whisk_heat,[],1);
[~,w_idx] = sort(w_max_idx);
data_squish{2} = unsorted__whisk_heat(:,w_idx)';

unsorted_heat_touch = norm_new(cell2mat(new_curves(touch_ix_idx)')');
[~,t_max_idx] = max(unsorted_heat_touch,[],1);
[~,t_idx_ix] = sort(t_max_idx);
data_squish{3} = unsorted_heat_touch(:,t_idx_ix)';

unsorted_whisk_heat = norm_new(cell2mat(new_whisk(whisk_ix_idx)')');
[~,w_max_idx] = max(unsorted_whisk_heat,[],1);
[~,w_idx_ix] = sort(w_max_idx);
data_squish{4} = unsorted_whisk_heat(:,w_idx_ix)';

unsorted_heat_touch = norm_new(cell2mat(new_curves(touch_ix_idx)')');
data_squish{5} = unsorted_heat_touch(:,w_idx_ix)';

unsorted_whisk_heat = norm_new(cell2mat(new_whisk(whisk_ix_idx)')');
data_squish{6} = unsorted_whisk_heat(:,t_idx_ix)';

for p = 1:numel(data_squish)
    subplot(3,2,p)
    [nr,nc] = size(data_squish{p});
    pcolor([data_squish{p} nan(nr,1); nan(1,nc+1)]);
    caxis([0 1])
    shading flat;
    colormap turbo
    if mod(p,2)==0
        set(gca,'xlim',[1 whisk_stretch_bins],'xtick',[])
    else
        set(gca,'xlim',[1 touch_stretch_bins],'xtick',[])
    end
    title(title_names{p})
end

%% histograms
[~,touch_peaks] = max(data{1},[],2);
[~,whisk_peaks] = max(data{2},[],2);
[~,touch_peaks_squish] = max(data_squish{1},[],2);
[~,whisk_peaks_squish] = max(data_squish{2},[],2);
figure(52);clf
subplot(2,2,1);histogram(touch_x(touch_peaks),linspace(touch_x(1),touch_x(end),8))
set(gca,'ylim',[0 25])
title('touch')
subplot(2,2,3);histogram(whisk_x(whisk_peaks),linspace(whisk_x(1),whisk_x(end),8),'facecolor','c')
set(gca,'ylim',[0 25])
title('whisk')
subplot(2,2,2);histogram(touch_peaks_squish,linspace(1,touch_stretch_bins,8))
set(gca,'ylim',[0 25])
title('touch squish')
subplot(2,2,4);histogram(whisk_peaks_squish,linspace(1,whisk_stretch_bins,8),'facecolor','c')
set(gca,'ylim',[0 25])
title('whisk squish')

