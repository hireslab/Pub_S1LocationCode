mouseNum = 'AH0716';
date = '170902';
cellName = 'JC1403';
cellNum = 'AAAB';
trialNum = '0194'

% mouseNum = 'AH0762';
% date = '171018';
% cellName = 'JC1505';
% cellNum = 'AAAB';
% trialNum = '0129'
% 
% mouseNum = 'AH0669';
% date = '170328';
% cellName = 'JC1239';
% cellNum = 'AAAA';
% trialNum = '0129'


% cd('C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays')
% load(['sweepArray_' mouseNum '_' date '_' cellName '_' cellNum])


cd(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\EPHUS\' mouseNum filesep cellName])

loadTrial = [cellName cellNum trialNum '.xsg'];
xsg_file = load([cellName cellNum trialNum '.xsg'], '-mat')


pulse_sets = xsg_file.header.stimulator.stimulator.pulseNameArray;
is_stim = any(strcmp('10_Hz_0',pulse_sets));

stim_params = xsg_file.header.stimulator.stimulator.pulseParameters{strcmp('10_Hz_0',pulse_sets)};

pulse_times = cumsum(repmat(stim_params.squarePulseTrainISI,stim_params.squarePulseTrainNumber,1)) - .1;
pulse_dur = pulse_times + stim_params.squarePulseTrainWidth;

pulse_on = ((stim_params.squarePulseTrainDelay + pulse_times) .* 10000) ;
pulse_off = ((stim_params.squarePulseTrainDelay + pulse_dur) .* 10000) ; 

figure(480);clf
if is_stim
    spike_trace = xsg_file.data.ephys.trace_1; 
    plot(spike_trace -mean(spike_trace) ,'k'); 
    for g = 1:length(pulse_on)
        hold on; plot([pulse_on(g) pulse_off(g)],[5 5],'b')
    end

    set(gca,'xtick',0:5000:40000,'xticklabel',0:.5:4,'xlim',[5000 11000])
end
% 
% 
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
    fn = 'sample_interneuron_zoom.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
    
%% load interneuron data struct
for i = 3
    % raster
    motors = normalize_var(U{i}.meta.motorPosition,1,-1);
    spikes = squeeze(U{i}.R_ntk);
    [~,sidx] = sort(motors);
    sidx = fliplr(sidx);
    
    figure(38);clf
    for g = 1:length(motors)
        spikeIdx = find(spikes(:,sidx(g)));
        hold on; scatter(spikeIdx,ones(numel(spikeIdx),1).*g,'k.')
    end
    set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
    axis square
end

    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
    fn = 'sample_interneuron_raster.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])