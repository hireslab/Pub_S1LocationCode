function plot_cell_depth(U, touch_struct)

location_cells = cellfun(@(x) x.is_tuned==1,touch_struct.pole);

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

% fn = 'scatter_depth_firingrate.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
