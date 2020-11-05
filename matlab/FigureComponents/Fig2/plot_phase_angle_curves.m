function plot_phase_angle_curves(whisk_struct)

angle_tuned = find(cellfun(@(x) x.is_tuned==1,whisk_struct.angle));
tmp = cellfun(@(x) range(x.stim_response.bars_fit.mean),whisk_struct.phase(angle_tuned));
[~,idx] = sort(tmp);

figure(480);clf
check_cells = idx(end-10:end);

for b = 1:numel(check_cells)
    cellNum = angle_tuned(check_cells(b));
    ix = whisk_struct.angle{cellNum}.stim_response.bars_stim;
    iy = whisk_struct.angle{cellNum}.stim_response.bars_fit.mean;
    i_err = whisk_struct.angle{cellNum}.stim_response.bars_fit.confBands(:,2) - iy;
    x = cellfun(@median,(whisk_struct.angle{cellNum}.stim_response.raw_stim));
    y = cellfun(@mean,(whisk_struct.angle{cellNum}.stim_response.raw_response));
    
    ix_p = whisk_struct.phase{cellNum}.stim_response.bars_stim;
    iy_p = whisk_struct.phase{cellNum}.stim_response.bars_fit.mean;
    i_err_p = whisk_struct.phase{cellNum}.stim_response.bars_fit.confBands(:,2) - iy_p;
    x_p = cellfun(@median,(whisk_struct.phase{cellNum}.stim_response.raw_stim));
    y_p = cellfun(@mean,(whisk_struct.phase{cellNum}.stim_response.raw_response));
    
    
    subplot(2,numel(check_cells),b)
    hold on; bar(x,y,'k')
    shadedErrorBar(ix,iy,i_err,'c')

    subplot(2,numel(check_cells),numel(check_cells)+b)
    hold on;bar(x_p,y_p,'k')
    shadedErrorBar(ix_p,iy_p,i_err_p,'r')
    set(gca,'xlim',[-pi pi],'xtick',[])
end

suptitle('top = angle; bottom = phase')

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = ['matched_tuning_curves.eps'];
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])