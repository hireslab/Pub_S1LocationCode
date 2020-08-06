function PlotPhaseAngleHeat(whisk_struct)

% find tuned units to either phase or angle
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.angle);
whisk_phase_tuned= cellfun(@(x) x.is_tuned==1,whisk_struct.phase);
whisk_cells_idx = {find(whisk_angle_tuned),  find(whisk_phase_tuned)};

% not normalized heatmap for angle 
bars_io = cellfun(@(x) [x.stim_response.bars_stim' x.stim_response.bars_fit.mean ],whisk_struct.angle(whisk_angle_tuned),'uniformoutput',0);
whisk_sampled = round([min(cellfun(@(x) min(x(:,1)),bars_io)):max(cellfun(@(x) max(x(:,1)),bars_io))]);
interped_resp = cellfun(@(x) interp1(x(:,1),x(:,2),whisk_sampled),bars_io,'uniformoutput',0);
norm_resp = cellfun(@(x) norm_new(x),interped_resp,'uniformoutput',0);

[~,maxidx] = cellfun(@(x) max(x),norm_resp);
[~,sorted_idx] = sort(maxidx);

figure(45);clf
pcolor(cell2mat(norm_resp(sorted_idx)')); %flipud because pcolor flips compared to imagesc
set(gca,'xtick',4:10:numel(whisk_sampled),'xticklabel',whisk_sampled(4:10:end))
colormap turbo
colorbar
caxis([0 1])
xlabel('whisker position');
title('SF 2: unshifted heat map')

% space normalized heatmap for phase and angle
hvar_names = {'angle','phase'};
whisk_structs = {whisk_struct.angle,whisk_struct.phase};
clear whisk_heat
figure(46);clf
for e = 1:numel(whisk_structs)
    
    for d = 1:numel(whisk_cells_idx{e})
        curr_w = whisk_structs{e}{whisk_cells_idx{e}(d)}.stim_response.values;
        %clean nan rows
        curr_w = curr_w(~any(isnan(curr_w),2),:);
        
        if strcmp(hvar_names{e},'pole')
            whisk_x = -2:.1:2;
        elseif strcmp(hvar_names{e},'phase')
            whisk_x = linspace(-pi,pi,21);
        elseif strcmp(hvar_names{e},'angle')
            whisk_x = -30:80;
        else
            whisk_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
        end
        whisk_heat{e}{d} = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
    end
    
    %this is the part for squishing the map for angle; 
    if strcmp(hvar_names{e},'angle')
%         whisk_stretch_bins = median(cellfun(@(x) sum(~isnan(x)),whisk_heat{e}));
        whisk_stretch_bins = 20; 
        raw_whisk = cellfun(@(x) x(~isnan(x)),whisk_heat{e},'uniformoutput',0);
        new_whisk = cellfun(@(x) interp1(linspace(1,whisk_stretch_bins,numel(x)),x,1:whisk_stretch_bins),raw_whisk,'uniformoutput',0);
        whisk_heat{e} = new_whisk;
    end
    
    full_map = norm_new(cell2mat(whisk_heat{e}')');
    full_map = [nan(1,size(full_map,2)) ; full_map ; nan(1,size(full_map,2))]; %pad edges since pcolor trims off edges. 
    [~,idx] = max(full_map);
    [~,sort_idx] = sort(idx);
    subplot(1,2,e)
    pcolor(full_map(:,sort_idx)');
    colormap turbo
    title(hvar_names{e})
    set(gca,'xtick',[]);
end
suptitle('space normalized heatmap')

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'whisker_position_pop.eps'
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])