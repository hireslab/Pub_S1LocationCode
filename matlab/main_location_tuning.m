%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells


%% Top level parameters and definitions
% U = defTouchResponse(U,.95,'off');
clearvars -except U
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';

pole_tuned = object_location_quantification(U,selectedCells,'pole','off'); %for old see object_location_v1.0
% pole_tuned = object_location_quantification(U,1:length(U),'pole','off'); %for old see object_location_v1.0

%% Raster + PSTH one cell (A) 
for i = 29
    %raster
    motors = normalize_var(U{i}.meta.motorPosition,1,-1);
    spikes = squeeze(U{i}.R_ntk);
    [~,sidx] = sort(motors);
    sidx = fliplr(sidx);
    
    figure(38);clf
    for g = 1:length(motors)
        subplot(1,2,1)
        spikeIdx = find(spikes(:,sidx(g)));
        hold on; scatter(spikeIdx,ones(numel(spikeIdx),1).*g,'k.')
    end
    set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
    
    
    touchOn = [find(U{i}.S_ctk(9,:,:)==1)  ;find(U{i}.S_ctk(12,:,:)==1)];
    touchOff = [find(U{i}.S_ctk(10,:,:)==1)  ;find(U{i}.S_ctk(13,:,:)==1)];
    touch_matrix = nan(size(spikes));
    for g = 1:length(touchOn)
        touch_matrix(touchOn(g):touchOff(g)) = 1;
    end
    
    touch_matrix = touch_matrix(:,sidx)';
    subplot(1,2,2)
    pcolor(touch_matrix)
    set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
    
    fn = 'example_raster.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
    
    %chunked psth
    sorted_spike_mat = spikes(:,sidx);
    chunk_one = sidx(1:round(numel(sidx)/3));
    chunk_two = sidx(round(numel(sidx)/3)+1: round(numel(sidx)/3)+1 + round(numel(sidx)/3));
    chunk_three = setdiff(sidx,[chunk_one chunk_two]);
    all_chunks = {chunk_three, chunk_two,chunk_one};
    
    figure(49);clf
    colors = [.3 .3 .3];
    for b = 1:length(all_chunks)
        hold on;plot(smooth(mean(spikes(:,all_chunks{b}),2)*1000,100),'color',colors.*b)
        set(gca,'ylim',[0 40],'xtick',0:1000:4000,'ytick',0:10:40)
        title(['motor pos = ' num2str(mean(motors(all_chunks{b})))])
    end
    title(['far:close ' num2str(cellfun(@(x) mean(motors(x)),all_chunks))])
    
    fn = 'example_psth.eps';
    export_fig([saveDir, fn], '-depsc ', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
end
%% firing rate X depth of recording (B)
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
location_cells = cellfun(@(x) x.is_tuned==1,pole_tuned);

jc_silent_cell = [766 819 895 631 776 815 910 871 844 902  941   840   888   748   732   940   686   944   950   933]; %Phils  from 902

figure(480);clf
subplot(3,1,[1 2])
scatter( cellfun(@(y) y.meta.depth,U),cellfun(@(y) mean(y.R_ntk(:))*1000, U),'k','filled')
hold on;scatter( cellfun(@(y) y.meta.depth,U(location_cells)),cellfun(@(y) mean(y.R_ntk(:))*1000, U(location_cells)),'g','filled')
hold on; scatter(jc_silent_cell,ones(length(jc_silent_cell),1)*-1,[],'c','filled');
hold on; plot([700 700],[0 30],'--k')
hold on; plot([900 900],[0 30],'--k')
title(['n = response (' num2str(numel(U)) ') + silent (' num2str(numel(jc_silent_cell)) ')' ])
set(gca,'ytick',0:10:30,'xtick',600:100:1000,'xlim',[600 1000])
ylabel('firing rate (hz)');
xlabel('depth (um)')


hold on; subplot(3,1,3)
histogram(cellfun(@(y) y.meta.depth,U(setdiff(1:length(U),location_cells))),600:25:1000,'facecolor','k')
hold on; histogram(cellfun(@(y) y.meta.depth,U(location_cells)),600:25:1000,'facecolor','g')
hold on;histogram(jc_silent_cell,600:25:1000,'facecolor','c')
set(gca,'xtick',600:100:1000,'xlim',[600 1000])

fn = 'scatter_depth_firingrate.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% 
% figure(48);clf
% subplot(1,2,1);
% plot(angle_tuned{51}.stim_response.bars_fit.x,angle_tuned{51}.stim_response.bars_fit.mean)
% subplot(1,2,2);
% plot(angle_tuned{115}.stim_response.bars_fit.x,angle_tuned{115}.stim_response.bars_fit.mean)

%% touch psth by quartiles of far, close and near (C)
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tuned_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));

touch_window = -25:50;
chunks = 3;
units_to_plot = [11 ,30, 6];
figure(20);clf
for g = 1:length(units_to_plot)
    selected_unit = tuned_units(units_to_plot(g));
    
    motors = normalize_var(U{selected_unit}.meta.motorPosition,1,-1);
    leftover = mod(length(motors),chunks);
    new_motors = datasample(motors,numel(motors)-leftover,'Replace',false);
    [s_motors,sidx] = sort(new_motors);
    all_chunks = reshape(s_motors,numel(new_motors)./chunks,[]);
    
    
    [tVar] = atTouch_sorter(U{selected_unit},touch_window);
    touch_motors = normalize_var(tVar.allTouches.S_ctk(:,end),-1,1);
    
    chunked_idx = cell(1,size(all_chunks,2));
    chunked_responses = cell(1,size(all_chunks,2));
    figure(20);
    colors = [.33 .33 .33];
    subplot(3,1,g)
    for b = 1:size(all_chunks,2)
        [rsmall,rbig] = bounds(all_chunks(:,b));
        chunked_idx = intersect(find(touch_motors<rbig),find(touch_motors>rsmall));
        chunked_responses{b} = tVar.allTouches.R_ntk(chunked_idx,:);
        
        hold on; plot(touch_window,smooth(nanmean(chunked_responses{b}).*1000,10),'color',colors.*b)
        
    end
    
    set(gca,'xtick',-25:25:50,'xlim',[-25 50])
end

fn = 'touch_psth_by_location.eps';
export_fig([saveDir, fn], '-depsc ', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% heatmap for object location tuned touch units (D)
variable = 'pole';
tuned_structs = pole_tuned(cellfun(@(x) x.is_tuned==1,pole_tuned));
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

fn = 'heat_map_tuning.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% modulation width (E)
% pole_tuned = object_location_quantification(U,selectedCells,'pole','off'); %for old see object_location_v1.0
tuned_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));

peak_response = cellfun(@(x) x.calculations.tune_peak,pole_tuned(tuned_units));
interp_norm_y = cell(1,length(peak_response));
for g = 1:length(peak_response)
    sr = pole_tuned{tuned_units(g)}.stim_response;
    centered_x = sr.values(:,1) - peak_response(g) ;
    norm_y = norm_new(sr.values(:,2));
    
    interp_centered_x = -2:.1:2;
    raw_x = round(centered_x,2);
    [~,idx] = unique(raw_x);
    interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
    figure(30);
    hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
end

pop_mean = nanmean(cell2mat(interp_norm_y'));
pop_sem = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
hold on; shadedErrorBar(interp_centered_x,pop_mean,pop_sem,'k')

set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
xlabel('distance from peak (mm)')
axis square

figure(30);
fn = 'all_modulation_width.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% Follicle at first touch vs later (SF)
pxpermm = 33;
clear fdist_map
for i = 1:numel(U)
    for b = 1:U{i}.k
        touchOn = [find(U{i}.S_ctk(9,:,b)==1) find(U{i}.S_ctk(12,:,b)==1)];
        if numel(touchOn)>1
            fx = U{i}.whisker.follicleX{b};
            fy = U{i}.whisker.follicleY{b};
            touchOn = touchOn(touchOn < numel(fx));
            fxy = [fx(touchOn)' fy(touchOn)'];
            fdist_map_x{i}{b} = fxy(2:end,1) - fxy(1,1);
            fdist_map_y{i}{b} = fxy(2:end,2) - fxy(1,2);
%             fdist = pdist2(fxy,fxy);
%             fdist_map{i}{b} = diag(fdist,1)';
        end
    end
end

num_touches_per_trial = cellfun(@(x) cellfun(@(y) numel(y),x), fdist_map_x,'uniformoutput',0); 
touch_num_cdf = cellfun(@(x) cumsum(histc(x,0:50)./numel(x)),num_touches_per_trial,'uniformoutput',0);
indiv = cell2mat(touch_num_cdf');
pop = mean(indiv);
pop_num_sd = std(indiv);
pop__num_sem = std(indiv)./sqrt(numel(U)); 

figure(80);clf
subplot(1,3,1)
plot(0:50,indiv,'color',[.8 .8 .8])
hold on;shadedErrorBar(0:50,pop,pop_num_sd,'r');
set(gca,'ylim',[0 1],'ytick',0:.25:1)
xlabel('number of touches per trial')
ylabel('proportion of all trials')

all_mean_fd = cellfun(@(x) nanmean(cell2nanmat(x),2), fdist_map_x,'uniformoutput',0);
pop_mean = nanmean(cell2nanmat(all_mean_fd),2) * -1;
pop_sd = nanstd(cell2nanmat(all_mean_fd),[],2);
pop_sem = nanstd(cell2nanmat(all_mean_fd),[],2) ./ sqrt(sum(~isnan(cell2nanmat(all_mean_fd)),2));
subplot(1,3,2)
errorbar(2:numel(pop_mean)+1,pop_mean./pxpermm,pop_sd./pxpermm,'-ko')
set(gca,'xlim',[0 20],'ylim',[-.5 .5],'ytick',-.5:.25:.5)
ylabel('anterior --- posterior')
xlabel('touch number')

all_mean_fd = cellfun(@(x) nanmean(cell2nanmat(x),2), fdist_map_y,'uniformoutput',0);
pop_mean = nanmean(cell2nanmat(all_mean_fd),2);
pop_sd = nanstd(cell2nanmat(all_mean_fd),[],2);
pop_sem = nanstd(cell2nanmat(all_mean_fd),[],2) ./ sqrt(sum(~isnan(cell2nanmat(all_mean_fd)),2));
subplot(1,3,3)
errorbar(2:numel(pop_mean)+1,pop_mean./pxpermm,pop_sd./pxpermm,'-ko')
set(gca,'xlim',[0 20],'ylim',[-.5 .5],'ytick',-.5:.25:.5)
ylabel('medial --- lateral')
xlabel('touch number')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
fn = 'scatter_depth_firingrate.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])







