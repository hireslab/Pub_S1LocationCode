%% Load whisking structures 
clear whisk_struct
hilbert_feature = {'angle','phase','midpoint','amplitude','velocity'};

for b = 1:numel(hilbert_feature)
    fileName = ['whisk_' hilbert_feature{b} '_window'];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'])
        whisk_struct.(hilbert_feature{b}) = wStruct;
   end
end
tuned_units = cellfun(@(x) x.is_tuned,whisk_struct.angle)==1;

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';

%% Build whisking structures if needed
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil

wStruct= whisking_location_quantification(U,1:length(U),'velocity','off');

%% Figure 2A: whisking vs quiet comparison
WhiskQuietComparison(U)

%% Figure 2B: angle and phase tuning curves
PlotPhaseAngleCurves(whisk_struct)

fn = ['matched_tuning_curves.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Figure 2C: population heatmap of phase and angle tuning
PlotPhaseAngleHeat(whisk_struct)

fn = 'whisker_position_pop.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% Figure 2E: Proportion of units tuned to phase, angle, neither, or both
PlotPhaseAnglePie(whisk_struct)

fn = 'proportion_angle_phase_pie.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% Figure 2F/G: Phase and angle modulation comparison 
PlotPhaseAngleModulation(whisk_struct)

%% Figure 2H-J: Hilbert decomposition comparison 
WhiskHilbertComparison(U,whisk_struct)

%% example whisker trace (SF)
for k = 33
    currentCell = U{k};
    
    %      trialSample = [datasample(1:currentCell.k,5)];
    %      trialSample = [134 136 datasample(1:currentCell.k,1)]; %CELL 83
    %         trialSample = [99 136 datasample(1:currentCell.k,3)]; %% CELL 22
    trialSample = [119 98 43 datasample(1:currentCell.k,2)]; %% CELL 33 EXAMPLE 119 SO GOOD!CELL TUNED TO ANGLE 10 . Dam, it's also phase tuned. LULz
    %         trialSample = [60 1 datasample(1:currentCell.k,1)]; %% CELL 74 maybe 96 too
    figure(2310);clf
    for g = 1:length(trialSample)
        
        whisk = currentCell.S_ctk(1,:,trialSample(g));
        midpoint = currentCell.S_ctk(4,:,trialSample(g));
        amp = currentCell.S_ctk(3,:,trialSample(g));
        amplitude = whisk;
        amplitude(amp<5) = nan;
        touchesOn = [find(currentCell.S_ctk(9,:,trialSample(g))==1) find(currentCell.S_ctk(12,:,trialSample(g))==1)];
        touchesOff = [find(currentCell.S_ctk(10,:,trialSample(g))==1) find(currentCell.S_ctk(13,:,trialSample(g))==1)];
        short_touches = (touchesOff-touchesOn)<=5;
        fill_values = touchesOff(short_touches)+5;
        touchesOff(short_touches) = fill_values;
        
        
        if ~isempty(touchesOn)
            for p = 1:numel(touchesOn)
                amplitude(touchesOn(p)+1:touchesOff(p)+30) = nan;
            end
        end
        
        spikes = find(currentCell.R_ntk(1,:,trialSample(g))==1);
        subplot(length(trialSample),1,g)
        plot(1:currentCell.t,whisk,'k')
        hold on;plot(1:currentCell.t,amplitude,'c')
        
        hold on; scatter(spikes,whisk(spikes),'ko','filled','markerfacecolor',[.8 .8 .8])
        hold on; scatter(spikes,amplitude(spikes),'ko','filled')
        
        if ~isempty(touchesOn) && ~isempty(spikes)
            for f = 1:length(touchesOn)
                tp = touchesOn(f):touchesOff(f);
                hold on; plot(tp,whisk(tp),'r-');
            end
        end
        
        title(['cell num ' num2str(k) ' at trial num ' num2str(trialSample(g))])
        set(gca,'ylim',[-30 70])
    end

end
fn = 'example_traces .eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% firing rate X depth of recording (SF)
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
whisk_cells =  cellfun(@(x) x.is_tuned,whisk_struct.angle)==1;
jc_silent_cell = [766 819 895 631 776 815 910 871 844 902   941   840   888   748   732   940   686   944   950   933]; %Phils  from 902

figure(480);clf
subplot(3,1,[1 2])
scatter( cellfun(@(y) y.meta.depth,U),cellfun(@(y) mean(y.R_ntk(:))*1000, U),'k','filled')
hold on;scatter( cellfun(@(y) y.meta.depth,U(whisk_cells)),cellfun(@(y) mean(y.R_ntk(:))*1000, U(whisk_cells)),'c','filled')
hold on; scatter(jc_silent_cell,ones(length(jc_silent_cell),1)*-1,[],'filled','markerfacecolor',[.6 .6 .6]);
hold on; plot([700 700],[0 30],'--k')
hold on; plot([900 900],[0 30],'--k')
title(['n = response (' num2str(numel(U)) ') + silent (' num2str(numel(jc_silent_cell)) ')' ])
set(gca,'ytick',0:10:30,'xtick',600:100:1000,'xlim',[600 1000])
ylabel('firing rate (hz)');
xlabel('depth (um)')


hold on; subplot(3,1,3)
histogram(cellfun(@(y) y.meta.depth,U(setdiff(1:length(U),whisk_cells))),600:25:1000,'facecolor','k')
hold on; histogram(cellfun(@(y) y.meta.depth,U(whisk_cells)),600:25:1000,'facecolor','c')
hold on;histogram(jc_silent_cell,600:25:1000,'facecolor',[.6 .6 .6])
set(gca,'xtick',600:100:1000,'xlim',[600 1000])

fn = 'scatter_depth_firingrate.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% CDF of sampling (SF)
built = find(cellfun(@(x) isfield(x,'stim_response'),whisk_struct.angle));
bars_built = cellfun(@(x) isfield(x.stim_response,'bars_fit'),whisk_struct.angle(built));
full_build = built(bars_built);

masks = cellfun(@(x) maskBuilder(x),U(full_build),'uniformoutput',0);

thetas = cellfun(@(x,y) squeeze(x.S_ctk(1,:,:)) .* y.whisking,U(full_build),masks,'uniformoutput',0);

hist_thetas = cellfun(@(x) histcounts(x(:),-25.5:80.5),thetas,'uniformoutput',0);

cdf_sums = cell2mat(cellfun(@(x) cumsum(x),hist_thetas','uniformoutput',0))';
cdf_mat = cdf_sums./(cdf_sums(end,:));

mean_time_per_trial_whisking = cellfun(@(x) nansum(x.whisking(:))./size(x.whisking,2)./1000,masks); % in seconds
num_bins_per_cell = cellfun(@(x) numel(x.stim_response.raw_stim),whisk_struct.angle(full_build));

figure(23);clf
subplot(2,3,[1 2 4 5])
plot(repmat(-25:80,size(cdf_mat,2),1)',cdf_mat,'color',[.9 .9 .9])
hold on; shadedErrorBar(-25:80,nanmean(cdf_mat,2),nanstd(cdf_mat,[],2)./sqrt(numel(full_build)))
axis square
set(gca,'ytick',0:.2:1,'xtick',-20:20:80)
xlabel('whisker angle during whisking')
ylabel('cumulative proportion of whisks')

subplot(2,3,3);
scatter(ones(1,numel(full_build)),mean_time_per_trial_whisking,'k','markeredgecolor',[.9 .9 .9])
hold on; errorbar(1,mean(mean_time_per_trial_whisking),std(mean_time_per_trial_whisking),'ko')
set(gca,'ylim',[0 3],'ytick',0:1:3,'xtick',[])
ylabel('time spent whisking per trial (s)');
axis square
title(['mean=' num2str(mean(mean_time_per_trial_whisking)) ' , SD=' num2str(std(mean_time_per_trial_whisking))])

subplot(2,3,6);
scatter(ones(1,numel(full_build)),num_bins_per_cell,'k','markeredgecolor',[.9 .9 .9])
hold on; errorbar(1,mean(num_bins_per_cell),std(num_bins_per_cell),'ko')
set(gca,'ylim',[0 80],'ytick',0:25:100,'xtick',[])
ylabel('number of whisking bins');
axis square
title(['mean=' num2str(mean(num_bins_per_cell)) ' , SD=' num2str(std(num_bins_per_cell))])

fn = 'whisk_sampling.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% phase map (SF)
hilbert_feature = 'phase';

tuned_units = cellfun(@(x) x.is_tuned,whisk_struct.phase)==1;

preferred_phase = cellfun(@(x) x.calculations.tune_peak,whisk_struct.phase(tuned_units));
phase_mod = cellfun(@(x) x.calculations.mod_idx_relative,whisk_struct.phase(tuned_units));
retractions = sign(preferred_phase)<0;
conversions = zeros(1,numel(preferred_phase));
conversions(retractions)=180;
deg_per_pie_unit =  180/pi;

preferred_degs = (abs(preferred_phase) .* deg_per_pie_unit) + conversions;

figure(380);clf
polarplot(preferred_phase,phase_mod,'o');

fn = 'phase_map.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
