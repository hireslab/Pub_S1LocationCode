function heatmap_builder(touch_struct,variable)

tuned_structs = touch_struct.(variable)(cellfun(@(x) x.is_tuned==1,touch_struct.(variable)));
touch_heat = cell(1,numel(tuned_structs));

for g = 1:numel(tuned_structs)
    curr_t = tuned_structs{g}.stim_response.values;
    curr_t = curr_t(~any(isnan(curr_t),2),:);%clean nan rows
    [~,u_idx] = unique(curr_t(:,1)); %catch non-unique x-values
    curr_t = curr_t(u_idx,:);
    if strcmp(variable,'pole')
        touch_x = -1:.1:1;
    elseif strcmp(variable,'angle')
        touch_x = linspace(-30,60,21);
    end
    
    touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
end
unsorted_heat = norm_new(cell2mat(touch_heat')');
[~,t_max_idx] = max(unsorted_heat,[],1);
[~,t_idx] = sort(t_max_idx);
data = unsorted_heat(:,t_idx)';
%set unsampled heatmap to nan;
[nr,nc] = size(data);
figure(15);clf
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
title(['location tuned heat map (n=' num2str(size(data,1)) ')'])
colorbar
% 
% fn = 'heat_map_tuning.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
