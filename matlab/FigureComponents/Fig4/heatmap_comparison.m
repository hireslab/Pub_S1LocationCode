function heatmap_comparison(U,touch_struct)

tuned_units = cellfun(@(x) x.is_tuned==1,touch_struct.pole);
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
trained = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

naive_tuned_units = intersect(find(tuned_units),find(naive));
expert_tuned_units = intersect(find(tuned_units),find(trained));

tuned_units = {naive_tuned_units,expert_tuned_units};

figure(80);clf
touch_map_all=[];
for b=1:length(tuned_units)
    sel_tstructs = touch_struct.pole(tuned_units{b});
    
    touch_heat = cell(1,length(sel_tstructs));
    %building heatmap for object location tuned touch units
    for g = 1:numel(sel_tstructs)
        curr_t = sel_tstructs{g}.stim_response.values;
        curr_t = curr_t(~any(isnan(curr_t),2),:);%clean nan rows
        [~,u_idx] = unique(curr_t(:,1)); %catch non-unique x-values
        curr_t = curr_t(u_idx,:);
        touch_x = -1:.1:1;
        touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
    end
    
    %plotting all touch object location tuned units
    subplot(2,2,b)
    unsorted_heat = norm_new(cell2mat(touch_heat')');
    [~,t_max_idx] = max(unsorted_heat,[],1);
    [~,t_idx] = sort(t_max_idx);
    data = unsorted_heat(:,t_idx)';
    
    %set unsampled heatmap to nan;
    [nr,nc] = size(data);
    pcolor([nan(nr,1) data nan(nr,1); nan(1,nc+2)]);
    shading flat;
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)+4],'xticklabel',-1:1:1,'ydir','reverse')
    axis square
    
    subplot(2,2,4)
    [~,idx] = max(data,[],2);
    x_conv = linspace(1,-1,21);
    hold on; histogram(x_conv(idx),-1.1:.5:1.1,'normalization','probability')
    set(gca,'ylim',[0 .6],'ytick',0:.3:.6)
    legend('naive','expert')
    
    pref{b} = x_conv(idx);
    
    touch_map_all = [touch_map_all ; data];
end
colormap turbo
subplot(2,2,3)
[nr,nc] = size(touch_map_all);
pcolor([nan(nr,1) touch_map_all nan(nr,1); nan(1,nc+2)]);
shading flat;
set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)+2],'xticklabel',-1:1:1,'ydir','reverse')
axis square

% figure(80);
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
% fn = 'naive_expert_heat.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])