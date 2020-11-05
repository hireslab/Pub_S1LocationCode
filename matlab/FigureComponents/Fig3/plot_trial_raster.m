function plot_trial_raster(U, cell_number)

%raster
motors = normalize_var(U{cell_number}.meta.motorPosition,1,-1);
spikes = squeeze(U{cell_number}.R_ntk);
[~,sidx] = sort(motors);
sidx = fliplr(sidx);

figure(38);clf
for g = 1:length(motors)
    subplot(1,2,1)
    spikeIdx = find(spikes(:,sidx(g)));
    hold on; scatter(spikeIdx,ones(numel(spikeIdx),1).*g,'k.')
end
set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
title('trial raster sorted by pole location')
touchOn = [find(U{cell_number}.S_ctk(9,:,:)==1)  ;find(U{cell_number}.S_ctk(12,:,:)==1)];
touchOff = [find(U{cell_number}.S_ctk(10,:,:)==1)  ;find(U{cell_number}.S_ctk(13,:,:)==1)];
touch_matrix = nan(size(spikes));
for g = 1:length(touchOn)
    touch_matrix(touchOn(g):touchOff(g)) = 1;
end

touch_matrix = touch_matrix(:,sidx)';
subplot(1,2,2)
pcolor(touch_matrix)
set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
title('touch windows')

% fn = 'example_raster.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

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

% fn = 'example_psth.eps';
% export_fig([saveDir, fn], '-depsc ', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
