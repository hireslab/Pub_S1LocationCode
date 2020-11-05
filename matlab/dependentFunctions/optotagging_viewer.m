% mouseNum = 'AH0716';
% date = '170902';
% cellName = 'JC1403';
% cellNum = 'AAAB';
% trialNum = '0194'

mouseNum = 'AH0762';
date = '171018';
cellName = 'JC1505';
cellNum = 'AAAB';
trialNum = '0129'
% 
mouseNum = 'AH0669';
date = '170328';
cellName = 'JC1239';
cellNum = 'AAAA';
trialNum = '0129'


% cd('C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays')
% load(['sweepArray_' mouseNum '_' date '_' cellName '_' cellNum])


cd(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\EPHUS\' mouseNum filesep cellName])
spike_threshold = 3.5;
loadTrial = [cellName cellNum trialNum '.xsg'];
xsg_file = load([cellName cellNum trialNum '.xsg'], '-mat');

spike_trace = xsg_file.data.ephys.trace_1; 
filt_threshold = spike_trace-mean(spike_trace) > spike_threshold ;
spike_filtered = spike_trace .* filt_threshold;
spike_derivative = diff(spike_filtered);
spike_times_idx = numeric_pattern(sign(spike_derivative)', [1, -1]);
spike_times = zeros(1,length(spike_derivative));
spike_times(spike_times_idx) = 1;


pulse_sets = xsg_file.header.stimulator.stimulator.pulseNameArray;
is_stim = any(strcmp('10_Hz_0',pulse_sets));
stim_params = xsg_file.header.stimulator.stimulator.pulseParameters{strcmp('10_Hz_0',pulse_sets)};
pulse_times = cumsum(repmat(stim_params.squarePulseTrainISI,stim_params.squarePulseTrainNumber,1)) - .1;
pulse_dur = pulse_times + stim_params.squarePulseTrainWidth;
pulse_on = ((stim_params.squarePulseTrainDelay + pulse_times) .* 10000) ;
pulse_off = ((stim_params.squarePulseTrainDelay + pulse_dur) .* 10000) ; 

% plot stimulation trace
figure(480);clf
if is_stim
    plot_trace = spike_trace -mean(spike_trace);
    plot(plot_trace ,'k'); 
    hold on; scatter(spike_times_idx,plot_trace(spike_times_idx+1),'ro')
    for g = 1:length(pulse_on)
        hold on; plot([pulse_on(g) pulse_off(g)],[5 5],'b')
    end
    set(gca,'xtick',0:5000:40000,'xticklabel',0:.5:4,'xlim',[5000 11000])
end

% spike from pulse onset
start = -5;
stop = 25;
xtick_spacing = 5; 
sampling_window = round([(start*10): (stop*10)]+ pulse_on);
spike_window = plot_trace(sampling_window(:));
sampled_spikes = reshape(spike_window,size(sampling_window))';

figure(39);clf;subplot(2,1,1)
plot([(start*10): (stop*10)],sampled_spikes,...
    'color',[.8 .8 .8])
set(gca,'xtick',start*10:xtick_spacing*10:stop*10,...
    'xticklabel',start:xtick_spacing:stop)
hold on; plot([0 100],[5.5 5.5],'b')
xlabel('time from pulse onset (ms)')
ylabel('amplitude')
title('10 pulses and their responses')

% spike shape
spike_windows = round([-8:12] + spike_times_idx');
selected_spikes = plot_trace(spike_windows)';
subplot(2,1,2)
plot(selected_spikes,'color',[.8 .8 .8])
hold on; plot(mean(selected_spikes,2),'r')
ptt_ratio = max(selected_spikes) ./ abs(min(selected_spikes));
[~, max_idx] = max(selected_spikes);
[~, min_idx] = min(selected_spikes);
ptt_duration = mean((min_idx - max_idx) / 10 );

set(gca,'xtick',0:10:20,'xticklabel',0:1:2)
ylabel('amplitude')
xlabel('ms')

text(2, -.5, ['ptt duration = ' num2str(ptt_duration) 'ms'])
text(2, -1.5, ['ptt ratio = ' num2str(mean(ptt_ratio))])
suptitle('example interneuron')

% 
% 
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
%     fn = 'sample_interneuron_zoom.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])
    
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