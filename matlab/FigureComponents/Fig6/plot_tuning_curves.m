function plot_tuning_curves(touch_struct,whisk_struct,hilbert_var) 

tUnits = cellfun(@(x) x.is_tuned==1,touch_struct.(hilbert_var));
wUnits = cellfun(@(x) x.is_tuned==1,whisk_struct.(hilbert_var));

co_tuned = find(tUnits.*wUnits);
touch_only = setdiff(find(tUnits),co_tuned);
whisk_only = setdiff(find(wUnits),co_tuned);
tuned_mat = {co_tuned,touch_only,whisk_only};

num_cells_max = max(cellfun(@numel, tuned_mat));
t_tick = {'','','--'}; %tick labels for touch cells co-tune, touch only, and whisk only
w_tick = {'','--',''};
save_labels = {'co_tune','touch_only','whisk_only'};

for g = 1:length(save_labels)
    %find units that are built
    built_tx = find(cellfun(@(x) isfield(x,'stim_response'),touch_struct.(hilbert_var)));
    built_wx = find(cellfun(@(x) isfield(x,'stim_response'),whisk_struct.(hilbert_var)));
    tx_filt2 = find(cellfun(@(x) isfield(x.stim_response,'bars_fit'),touch_struct.(hilbert_var)(built_tx)));
    wx_filt2 = find(cellfun(@(x) isfield(x.stim_response,'bars_fit'),whisk_struct.(hilbert_var)(built_wx)));

    tx_ix = intersect(built_tx(tx_filt2),tuned_mat{g});
    wx_ix = intersect(built_wx(wx_filt2),tuned_mat{g});
    
    %build xyerr values based on bars fit for touch
    tx = cellfun(@(x) x.stim_response.bars_fit.x,touch_struct.(hilbert_var)(tx_ix),'uniformoutput',0);
    ty = cellfun(@(x) x.stim_response.bars_fit.mean,touch_struct.(hilbert_var)(tx_ix),'uniformoutput',0);
    nty = cellfun(@(x) normalize_var(x,0,1),ty,'uniformoutput',0); %normalized values to [0 1]
    terr = cellfun(@(x) x.stim_response.bars_fit.confBands,touch_struct.(hilbert_var)(tx_ix),'uniformoutput',0);
    ty_vals = {ty,nty};
    %build xyerr values based on bars fit for whisking
    wx = cellfun(@(x) x.stim_response.bars_fit.x,whisk_struct.(hilbert_var)(wx_ix),'uniformoutput',0);
    wy = cellfun(@(x) x.stim_response.bars_fit.mean,whisk_struct.(hilbert_var)(wx_ix),'uniformoutput',0);
    nwy = cellfun(@(x) normalize_var(x,0,1),wy,'uniformoutput',0);
    werr = cellfun(@(x) x.stim_response.bars_fit.confBands,whisk_struct.(hilbert_var)(wx_ix),'uniformoutput',0);
    wy_vals = {wy,nwy};
    
    %plotting in one giant elongated window
    figure(58);clf
    for k = 1:numel(wy_vals)
        for b = 1:numel(tuned_mat{g})
            if k ==1
                subplot(2,num_cells_max,b)
                try
                    shadedErrorBar(tx{b},ty_vals{k}{b},terr{b}(:,2) - ty_vals{k}{b},['r' t_tick{g}])
                end
                try
                    hold on; shadedErrorBar(wx{b},wy_vals{k}{b},werr{b}(:,2) - wy_vals{k}{b},['c' w_tick{g}])
                end
            elseif k ==2
                subplot(2,num_cells_max,b+num_cells_max)
                try
                    plot(tx{b},ty_vals{k}{b},['r' t_tick{g}])
                end
                try
                    hold on; plot(wx{b},wy_vals{k}{b},['c' w_tick{g}])
                end
            end
            if strcmp(hilbert_var,'phase')
                set(gca,'xlim',[-pi pi],'xtick',[])
                if k == 2
                    set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'})
                end
            end
        end   
    end
    
    suptitle(save_labels{g})
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
    fn = [hilbert_var save_labels{g} '_curves.eps'];
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
end