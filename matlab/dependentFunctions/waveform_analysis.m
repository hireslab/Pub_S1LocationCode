%% script used to build SFig1 D/E: peak to trough waveform analysis. 

data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Raw\';

T_directory = [data_directory 'Excitatory\TArrays_Naive\'];
[n_excit_ptt_ratio, n_excit_ptt_duration] = peak_to_trough_metrics(T_directory);

T_directory = [data_directory 'Excitatory\TArrays_Trained\'];
[t_excit_ptt_ratio, t_excit_ptt_duration] = peak_to_trough_metrics(T_directory);

T_directory = [data_directory 'Excitatory\TArrays_Phil\'];
[p_excit_ptt_ratio, p_excit_ptt_duration] = peak_to_trough_metrics(T_directory);

T_directory = [data_directory 'Interneurons\TArrays\'];
[inhib_ptt_ratio, inhib_ptt_duration] = peak_to_trough_metrics(T_directory);
%%
excit_ptt_ratio = [n_excit_ptt_ratio t_excit_ptt_ratio p_excit_ptt_ratio];
excit_ptt_duration = [n_excit_ptt_duration t_excit_ptt_duration p_excit_ptt_duration];

figure(480);clf;
subplot(1,3,3)
scatter(excit_ptt_ratio,excit_ptt_duration,'b','filled')
hold on; scatter(inhib_ptt_ratio,inhib_ptt_duration,'r','filled')
set(gca,'xlim',[0 8],'ylim', [0 2])
xlabel('PTT ratio')
ylabel('PTT duration (ms)')

subplot(1,3,1)
bar([1,2],[median(excit_ptt_ratio) median(inhib_ptt_ratio)],'k','facealpha',1)
hold on; errorbar(1,median(excit_ptt_ratio),std(excit_ptt_ratio)./sqrt(numel(excit_ptt_ratio)),'b')
hold on; errorbar(2,median(inhib_ptt_ratio),std(inhib_ptt_ratio)./sqrt(numel(inhib_ptt_ratio)),'r')
set(gca,'xlim',[0,3],'xtick',[1,2], 'xticklabel',{'excit','inhib'})
ylabel('peak to trough ratio')

subplot(1,3,2)
bar([1,2],[median(excit_ptt_duration) median(inhib_ptt_duration)],'k','facealpha',1)
hold on; errorbar(1,median(excit_ptt_duration),std(excit_ptt_duration)./sqrt(numel(excit_ptt_duration)),'b')
hold on; errorbar(2,median(inhib_ptt_duration),std(inhib_ptt_duration)./sqrt(numel(inhib_ptt_duration)),'r')
set(gca,'xlim',[0,3],'xtick',[1,2], 'xticklabel',{'excit','inhib'})
ylabel('peak to trough duration (ms')

suptitle('peak to trough waveform analysis')

 fn = 'peak_trough.eps';
 export_fig(fn, '-depsc', '-painters', '-r1200', '-transparent')
 fix_eps_fonts([ fn])
