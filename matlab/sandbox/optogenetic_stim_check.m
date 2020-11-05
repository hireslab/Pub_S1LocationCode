%% load data
T_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Raw\Excitatory\TArrays\';
S_directory = 'C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SweepArrays\';
ST_directory = 'C:\Users\jacheung\Dropbox\HLabBackup\Jon\DATA\SpikesData\SpikeArrays\';

% preallocate directory list and peak to trough ratio/duration vectors
T_array_list = dir([T_directory '/*.mat']);
ptt_ratio = zeros(1, length(T_array_list));
ptt_duration = zeros(1, length(T_array_list));

% plot parameters
rc = numSubplots(length(T_array_list));

for g = 1:length(T_array_list)
    disp(['calculating cell ' num2str(g) '/' num2str(length(T_array_list))])
    load([T_directory T_array_list(g).name])
    
    sweep_name = ['sweepArray_' T.mouseName '_' T.sessionName '_' T.cellNum '_' T.cellCode '.mat'];
    spikes_trials = ['spikes_trials_' T.mouseName '_' T.sessionName '_' T.cellNum '_' T.cellCode '.mat'];
    load([S_directory sweep_name])
    load([ST_directory spikes_trials])
    waveforms = s.get_spike_waveforms_all;

    % spike filter to grab indices of filtered spikes
    st_idx = cell(1,length(s.sweeps));
    for i = 1:length(s.sweeps)
        sts = spikes_trials.spikesTrials{i}.spikeTimes;
        raw_sts = s.sweeps{i}.spikeTimes;
        [~,st_idx{i}] = intersect(raw_sts,sts);
    end
    
    all_waveforms = cell2mat(cellfun(@(x,y) x.spikeWaveforms{2}(:,y),waveforms,st_idx,'uniformoutput',0));
    avg_waveform = mean(all_waveforms,2);
    peak_normalized_waves = (all_waveforms- mean(avg_waveform))';
    
    % calculate peak to trough ratio and duration
    [max_vals, max_idx] = max(peak_normalized_waves,[],2);
    [min_vals, min_idx] = min(peak_normalized_waves,[],2);
    ptt_ratio(g) = median(abs(max_vals ./ min_vals));
    ptt_duration(g) = median(min_idx - max_idx)/10;
    disp(['peak to trough ratio: ' num2str(ptt_ratio(g)) ' , duration: ' num2str(ptt_duration(g)) 'ms'])

    % plot
    figure(10); subplot(rc(1),rc(2),g)
    hold on; plot(mean(peak_normalized_waves),'r')
end
% 
% figure(10);clf
% plot(peak_normalized_waves,'color',[.9 .9 .9],'linewidth',.1)
% hold on; plot(apeak_normalized_waves,'r')
