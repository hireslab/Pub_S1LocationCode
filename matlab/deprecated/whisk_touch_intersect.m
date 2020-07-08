%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%%
hilbertVar = 'angle';

selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
%     fn = 'touch_all.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])

tStruct = object_location_quantification(U,selectedCells,hilbertVar,'off');
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
%     fn = 'touch_location_all.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])

fileName = ['whisk_' hilbertVar '_tune'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');
    
end

if strcmp(hilbertVar,'pole')
    population_heatmap_builder(tStruct,wStruct,hilbertVar)
    
    %     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
    %     fn = 'population_location.eps';
    %     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    %     fix_eps_fonts([saveDir, fn])
else
    disp('not building out population heatmaps. function not optimized for other variables')
end
%% scatter of tuning preference of whisk and touch

tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

touch_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],tStruct(tUnits),'uniformoutput',0)') ;
whisking_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],wStruct(wUnits),'uniformoutput',0)');

%scatter of whisking (Y) vs touch (X)
figure(3850);clf
if strcmp(hilbertVar,'pole')
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
elseif strcmp(hilbertVar,'phase')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-4,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*-4,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)

    set(gca,'xlim',[-4 4],'ylim',[-4 4],...
        'xtick',-pi:pi:pi,'ytick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'yticklabel',{'-\pi',0,'\pi'})
    hold on; plot([-4 4],[-4 4],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),-pi:pi/4:pi,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-pi:pi/4:pi,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
elseif strcmp(hilbertVar,'angle')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-30,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*-30,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)

    set(gca,'xlim',[-40 80],'ylim',[-40 80],...
        'xtick',-40:20:80,'ytick',-40:20:80)
    hold on; plot([-40 80],[-40 80],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),-40:10:80,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-40 80],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-40 80],'ylim',[0 20])
end

figure(3850);subplot(4,1,[1 2])
lm = fitlm(touch_pw(touch_ix_idx,1),whisking_pw(whisk_ix_idx,1));
predicts = lm.predict;
[s_vals,sort_idx] = sort(touch_pw(touch_ix_idx,1));
hold on; plot(s_vals,predicts(sort_idx))

legend('whisk tuned only','touch tuned only','both tuned')
axis square
xlabel('whisk tune peak');ylabel('touch tune peak')
title(['whisk=' num2str(numel(whisk_nonIX_idx)) ', touch=' num2str(numel(touch_nonIX_idx)) ', both=' num2str(numel(touch_ix_idx))])

figure(3850);
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = [hilbertVar '_intersect.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% % scatter of absolute modulation values
% touch_abs_mod = cell2mat(cellfun(@(x) x.calculations.mod_idx_abs,tStruct(tUnits),'uniformoutput',0)') ;
% whisk_abs_mod = cell2mat(cellfun(@(x) x.calculations.mod_idx_abs,wStruct(wUnits),'uniformoutput',0)') ;
% 
% [min_bound,max_bound] = bounds([touch_abs_mod(touch_ix_idx); whisk_abs_mod(whisk_ix_idx)]);
% 
% figure(2410);clf
% hold on; scatter(whisk_abs_mod(whisk_nonIX_idx),zeros(1,length(whisk_nonIX_idx)),'c');
% hold on; scatter(zeros(1,length(touch_nonIX_idx)),touch_abs_mod(touch_nonIX_idx),'b');
% hold on; scatter(whisk_abs_mod(whisk_ix_idx),touch_abs_mod(touch_ix_idx),'k')
% hold on; plot([0,100],[0,100],'--k')
% set(gca,'xlim',[0 100],'ylim',[0 100],'ytick',0:25:100,'xtick',0:25:100)
% axis square
% xlabel('whisk absolute modulation');
% ylabel('touch absolute modulation');
% 
% figure(2411);clf
% subplot(2,1,1)
% histogram(touch_abs_mod,0:5:100,'facecolor','b','facealpha',1)
% set(gca,'xlim',[0 100],'ylim',[0 10],'xtick',0:25:100,'ytick',0:5:10)
% 
% subplot(2,1,2)
% histogram(whisk_abs_mod,0:5:100,'facecolor','c','facealpha',1)
% set(gca,'xlim',[0 100],'ylim',[0 10],'xtick',0:25:100,'ytick',0:5:10)
%% intersection of whisking and touch
tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

whiskTuned = find(wUnits);
touchTuned = find(tUnits);
touch_whisk_tuned = intersect(find(tUnits),find(wUnits));

% touch_OL = cellfun(@(x) x.is_tuned==1,tStruct);
touch_OL = logical(ones(1,length(U)));
rc = numSubplots(sum(touch_OL));

sel_tstructs = tStruct(touch_OL);
sel_wstructs = wStruct(touch_OL);

figure(100);clf
figure(101);clf
whisk_touch_pair = cell(1,sum(touch_OL));
touch_diff_pair = cell(1,sum(touch_OL));
for g = 1:sum(touch_OL)
    if isfield(sel_wstructs{g},'stim_response') && isfield(sel_tstructs{g},'stim_response')
        curr_w = sel_wstructs{g}.stim_response.values;
        curr_t = sel_tstructs{g}.stim_response.values;
        
        curr_w = curr_w(~any(isnan(curr_w),2),:); %clean nan rows
        curr_t = curr_t(~any(isnan(curr_t),2),:);
        
        if strcmp(hilbertVar,'pole')
            whisk_x = round(round(min(curr_w(:,1)),1):.1:round(max(curr_w(:,1)),1),1);
            touch_x = round(round(min(curr_t(:,1)),1):.1:round(max(curr_t(:,1)),1),1);
        elseif strcmp(hilbertVar,'phase')
            whisk_x = linspace(-pi,pi,21);
            touch_x = linspace(-pi,pi,21);
        else
            whisk_x = round(round(min(curr_w(:,1))):1:round(max(curr_w(:,1))));
            touch_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
        end
        whisk_response = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
        whisk_std = interp1(curr_w(:,1),curr_w(:,3),whisk_x);
        whisk_CI = interp1(curr_w(:,1),curr_w(:,4),whisk_x);
        
        touch_response = interp1(curr_t(:,1),curr_t(:,2),touch_x);
        touch_std = interp1(curr_t(:,1),curr_t(:,3),touch_x);
        touch_CI =  interp1(curr_t(:,1),curr_t(:,4),touch_x);
        
        [~,~,whisk_idx] = intersect(touch_x,whisk_x);
        [overlap_x,~,touch_idx] = intersect(whisk_x,touch_x);
        
        %raw responses within touch ranges
        %     figure(99);subplot(rc(1),rc(2),g)
        %     shadedErrorBar(overlap_x,whisk_response(whisk_idx),whisk_CI(whisk_idx),'c')
        %     hold on; shadedErrorBar(overlap_x,touch_response(touch_idx),touch_CI(touch_idx),'b')
        %     if strcmp(hilbertVar,'pole')
        %         set(gca,'xlim',[-1 1],'xdir','reverse')
        %     elseif strcmp(hilbertVar,'phase')
        %         set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
        %     end
        
        %raw responses
        figure(100);subplot(rc(1),rc(2),g)
        shadedErrorBar(whisk_x,whisk_response,whisk_CI,'c')
        hold on; shadedErrorBar(touch_x,touch_response,touch_CI,'b')
        if strcmp(hilbertVar,'pole')
            set(gca,'xlim',[-1 2],'xdir','reverse')
            axis square
        elseif strcmp(hilbertVar,'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
        end
        
        if any(g == whiskTuned(whisk_nonIX_idx))
            title('whisk only')
        elseif any(g == touchTuned(touch_nonIX_idx))
            title('touch only')
        elseif any(g ==touch_whisk_tuned)
            title('T+W')
        end
        
        %normalized responses MAXMIN
%         m = min(whisk_response);
%         range = max(whisk_response) - m;
%         norm_whisk= (whisk_response - m) ./ range;
%         norm_whisk_CI = (whisk_CI - m) ./ range;
%         
%         m = min(touch_response);
%         range = max(touch_response) - m;
%         norm_touch= (touch_response - m) ./ range;
%         norm_touch_CI = (touch_CI - m) ./ range;
%         
        %zscoring
        mu = nanmean(whisk_response);
        sigma = nanstd(whisk_response);
        norm_whisk= (whisk_response - mu) ./ sigma;
        norm_whisk_CI = (whisk_CI - mu) ./ sigma;
        
        mu = nanmean(touch_response);
        sigma = nanstd(touch_response);
        norm_touch= (touch_response - mu) ./ sigma;
        norm_touch_CI = (touch_CI - mu) ./ sigma;
        

    
        figure(101);subplot(rc(1),rc(2),g)
            shadedErrorBar(whisk_x,norm_whisk,norm_whisk_CI,'c')
            hold on; shadedErrorBar(touch_x,norm_touch,norm_touch_CI,'b')
%         plot(whisk_x,norm_whisk,'c')
%         hold on; plot(touch_x,norm_touch,'b')
        if strcmp(hilbertVar,'pole')
            set(gca,'xlim',[-1 2],'xdir','reverse','ylim',[-3 3])
            axis square
        elseif strcmp(hilbertVar,'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
        end
        
        if any(g == whiskTuned(whisk_nonIX_idx))
            title('whisk only')
        elseif any(g == touchTuned(touch_nonIX_idx))
            title('touch only')
        elseif any(g ==touch_whisk_tuned)
            title('T+W')
        end
        
    end
    
    
    
end

%     figure(100);
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig5\';
%     fn = 'whisk_touch_tuning_curves.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])
% 
%     figure(101);
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig5\';
%     fn = 'whisk_touch_tuning_curves_normalized.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])



%% WHISK TUNING CURVES
whisk_units = find(wUnits);
rc = numSubplots(numel(whisk_units));
num_bins = 50;
figure(239);clf
for g = 1:numel(whisk_units)
    array = U{whisk_units(g)};
    
    spikes = squeeze(array.R_ntk);
    timePostTouchToTrim = 30;
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touchOnIdx = touchOnIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchOffIdx = touchOffIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(touchOnIdx)
        touchEx_mask(touchOnIdx(i):touchOffIdx(i)+timePostTouchToTrim) = NaN; %added 30 to add time from touch offset
    end
    touchEx_mask(1:100,:) = 1; %since bleedover from end of trials before, tihs ensure we keep end
    
    amplitude = squeeze(array.S_ctk(3,:,:));
    whisking = nan(size(squeeze(array.S_ctk(1,:,:))));
    whisking(amplitude>5)=1;
    
    whisking_mask = whisking .* touchEx_mask;
    
    
    if strcmp(hilbertVar,'phase')
        x_raw = squeeze(array.S_ctk(5,:,:)) .* whisking_mask;
    elseif strcmp(hilbertVar,'angle')
        x_raw = squeeze(array.S_ctk(1,:,:)) .* whisking_mask;
    end
    y_raw = spikes .* whisking_mask;
    
    x = x_raw(~any(isnan([x_raw(:) y_raw(:)]),2));
    y = y_raw(~any(isnan([x_raw(:) y_raw(:)]),2));
    
    if strcmp(hilbertVar,'phase')
        [sorted, sortedBy] = binslin(x,y,'equalE',num_bins+1,-pi,pi);
        y_bin = cellfun(@(x) nanmean(x),sorted) .* 1000;
        x_bin = linspace(-pi,pi,num_bins);
    elseif strcmp(hilbertVar,'angle')
        [sorted, sortedBy] = binslin(x,y,'equalE',num_bins+1,min(x),max(x));
        y_bin = cellfun(@(x) nanmean(x),sorted) .* 1000;
        x_bin = linspace(min(x),max(x),num_bins);
    end
    
    try
        f = fit(x_bin',y_bin,'gauss1');
        
        %     f = fit(x,y .* 1000,'gauss1');
        
        figure(239);subplot(rc(1),rc(2),g)
        scatter(x_bin,y_bin,'.k')
        yyaxis right
        hold on; plot(linspace(min(x),max(x),20),f(linspace(min(x),max(x),20)),'g');
        hold on; plot(wStruct{whisk_units(g)}.stim_response.values(:,1),wStruct{whisk_units(g)}.stim_response.values(:,2),'r')
        hold on; plot(wStruct{whisk_units(g)}.stim_response.bars_fit.x,wStruct{whisk_units(g)}.stim_response.bars_fit.mean,'b')
    end
end


