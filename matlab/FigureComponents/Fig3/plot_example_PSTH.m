function plot_example_PSTH(U ,touch_struct, units_to_plot)
tuned_units = find(cellfun(@(x) x.is_tuned==1,touch_struct.pole));

touch_window = -25:50;
chunks = 3;
figure(20);clf
colors = [.5 .75 .1];
for g = 1:length(units_to_plot)
    selected_unit = tuned_units(units_to_plot(g));
    
    motors = normalize_var(U{selected_unit}.meta.motorPosition,1,-1);
    leftover = mod(length(motors),chunks);
    new_motors = datasample(motors,numel(motors)-leftover,'Replace',false);
    [s_motors,~] = sort(new_motors);
    all_chunks = reshape(s_motors,numel(new_motors)./chunks,[]);

    [tVar] = atTouch_sorter(U{selected_unit},touch_window);
    touch_motors = normalize_var(tVar.allTouches.S_ctk(:,end),-1,1);

    chunked_responses = cell(1,size(all_chunks,2));
    figure(20);
    subplot(length(units_to_plot),1,g)
    for b = 1:size(all_chunks,2)
        [rsmall,rbig] = bounds(all_chunks(:,b));
        chunked_idx = intersect(find(touch_motors<rbig),find(touch_motors>rsmall));
        chunked_responses{b} = tVar.allTouches.R_ntk(chunked_idx,:);
        hold on; plot(touch_window,smooth(nanmean(chunked_responses{b}).*1000,10),'color',ones(1,3).*colors(b))
    end
    
    if g == length(units_to_plot)
        legend('far', 'middle','close')
    end
    
    xlabel('time from touch onset (ms)')
    ylabel('firing rate (Hz)')
    set(gca,'xtick',-25:25:50,'xlim',[-25 50])
end