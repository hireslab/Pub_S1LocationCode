%How does each stimulus vary w/ pole position?
%dependent function to id all touches and pre/post decision touches
preDecisionTouches = preDecisionTouchMat(U);
viewWindow = [-25:50];

gc_vals = cell(length(U),1);
%stimulus and response variables definitions
for rec = 1:27
    array = U{rec};
    [tVar] = atTouch_sorter(array,viewWindow,preDecisionTouches{rec});
    pole = tVar.allTouches.S_ctk(:,end);
    it = tVar.allTouches.S_ctk(:,1:end-1);
    dt = tVar.allTouches.dtS_ctk;
    
    %heat of just each feature correlated to pole
    corr_values = corr([pole it dt]);
    gc_vals{rec} = corr_values;
    
    %plotting all stimulus features for one single neuron in a heatmap
    if rec == 27
        [~,idx] = sort(pole);
        
        mod_it = normalize_var(it(idx,:),0,1) + repmat(1:size(it,2),length(it),1);
        mod_dt = normalize_var(dt(idx,:),0,1) + repmat(size(it,2)+1:size(it,2)+size(dt,2),length(it),1);
        mod_mat = [mod_it mod_dt];
        figure(280);clf
        scatter(1:length(pole),normalize_var(pole(idx),1,0),'b')
        for g = 1:size(it,2)+size(dt,2)
            hold on; plot(1:length(pole),mod_mat(:,g),'k')
        end
        
        set(gca,'ytick',.5:size([pole it dt],2)+.5,'yticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames],...
            'xlim',[1 length(pole)],'ydir','reverse')
        
        
        stretch_resolution = 100 ;
        stretch_map = redbluecmap; 
        new_map = nan(stretch_resolution,3); 
        for j = 1:3
        new_map(:,j) = interp1(linspace(1,stretch_resolution,length(stretch_map)),stretch_map(:,j),1:stretch_resolution);
        end
        
        figure(238);clf
        imagesc([normalize_var(pole(idx),1,0)' ; normalize_var([it(idx,:) dt(idx,:)],0,1)'])
        set(gca,'yticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames])
        colorbar
        colormap(new_map)
        xlabel('touches sorted by pole position')
    end
    
end

%% Population feature correlation
vect = abs(cell2mat(cellfun(@(x) x(:)',gc_vals,'uniformoutput',0)));
figure(67);clf
imagesc(reshape(nanmedian(vect),size(gc_vals{1})))

caxis([0 1])
colorbar
colormap bone
set(gca,'xticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames]...
    ,'yticklabel',['pole' tVar.allTouches.itNames(1:end-1) tVar.allTouches.dtNames])
axis square