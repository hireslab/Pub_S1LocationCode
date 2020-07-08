function directional_modulation_struct = touch_directional_selectivity(uberarray,selected_cells,alpha_value,displayOpt)


if (nargin < 4), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    figure(9);clf
end

for rec = 1:length(selected_cells)
    array = uberarray{(selected_cells(rec))};
    
    spks = squeeze(array.R_ntk(:,:,:));
    response_window = array.meta.touchProperties.responseWindow(1):array.meta.touchProperties.responseWindow(2);
    touch_times = [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    phase = squeeze(array.S_ctk(5,:,:));
    velocity = squeeze(array.S_ctk(2,:,:));
    pt_vel = mean(velocity(touch_times - [-5:-1]),2);
    touch_phase = phase(touch_times);
    
    touch_response = mean(spks(touch_times + response_window),2);
    pro_touches = intersect(find(touch_phase<0),find(pt_vel>0));
    ret_touches = intersect(find(touch_phase>0),find(pt_vel<0));
    
    pro_resp = max(touch_response(pro_touches));
    ret_resp = max(touch_response(ret_touches));
    directional_modulation_struct{(selected_cells(rec))}.mod_idx_relative = (pro_resp-ret_resp) ./  (pro_resp+ret_resp);
    
    pro_responses = touch_response(pro_touches)*1000;
    ret_responses = touch_response(ret_touches)*1000;
    [~,p] = ttest2(pro_responses,ret_responses);
    ret_error = std(ret_responses)./sqrt(numel(ret_touches));
    pro_error = std(pro_responses)./ sqrt(numel(pro_touches));
    hold on; errorbar(mean(ret_responses),mean(pro_responses),pro_error,pro_error,ret_error,ret_error,'k','CapSize',0)
    if willdisplay
        if p < alpha_value
            if mean(pro_responses)>mean(ret_responses)
                hold on; scatter(mean(ret_responses),mean(pro_responses),'filled','markerfacecolor','r')
            else
                hold on; scatter(mean(ret_responses),mean(pro_responses),'filled','markerfacecolor','b')
            end
        else
            hold on; scatter(mean(ret_responses),mean(pro_responses),'filled','markerfacecolor',[.8 .8 .8])
        end
    end
    
end

if willdisplay
    hold on; plot([0 100],[0 100],'--k')
    set(gca,'xlim',[0 100],'ylim',[0 100],'ytick',0:25:100,'xtick',0:25:100)
    axis square
    xlabel('retraction responses');ylabel('protraction responses')
end