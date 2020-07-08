function wt_modulationIndex(hilbertTouch,hilbertWhisking,selectedCells)

fieldsToCompare = fields(hilbertTouch.R_ntk.allTouches);
figure(23);clf

for g = 1:4
    figure(23);subplot(2,2,g)
    currTArray = hilbertTouch.R_ntk.allTouches.(fieldsToCompare{g});
    currWArray = hilbertWhisking.R_ntk.(fieldsToCompare{g});

    chosenTArrays = currTArray(selectedCells);
    TouchR = cellfun(@(x) nanmean(x,2),chosenTArrays,'uniformoutput',0);
    chosenWArrays = currWArray(selectedCells);
    WhiskR = cellfun(@(x) nanmean(x,1),chosenWArrays,'uniformoutput',0);

    touchModIdx = cellfun(@(x) (max(x) - min(x))./nanmean(x),TouchR,'uniformoutput',1);
    whiskModIdx = cellfun(@(x) (max(x) - min(x))./nanmean(x),WhiskR,'uniformoutput',1);
    
    
    %modulation index
    scatter(whiskModIdx,touchModIdx)
    set(gca,'xlim',round([min([touchModIdx whiskModIdx]) max([touchModIdx whiskModIdx])]),'ylim',round([min([touchModIdx whiskModIdx]) max([touchModIdx whiskModIdx])]),'ytick',0:2:50,'xtick',0:2:50)
    hold on; plot(round([min([touchModIdx whiskModIdx]) max([touchModIdx whiskModIdx])]),round([min([touchModIdx whiskModIdx]) max([touchModIdx whiskModIdx])]),'-.k')
    axis square
    xlabel('whisking mod idx');ylabel('touch mod idx')
    title([fieldsToCompare{g} ' modulation idx'])
    
    
end