function touchFeatureBinned_plotter(cellStruct,popV,selectedCells,variableFields,touchOrderFields,window,gaussFilt)

[rc] = numSubplots(length(selectedCells));

plotrow = rc(1);
plotcolumn = rc(2);

for g = 1:numel(touchOrderFields)
    figure(322+g);clf
    for d = 1:length(selectedCells)
       
        currCell = selectedCells(d);
        array = cellStruct{currCell};
        
        %toss touches out that are >10 degree difference from nearest bin.
        toKeep = diff(popV{currCell}.allTouches.(variableFields{1}).range)<10;
        if ~isempty(find(toKeep==0))
            toKeep(find(toKeep==0):end) = 0;
        end
        
        toKeep = logical(ones(length(popV{currCell}.(touchOrderFields{g}).theta.range),1));
        
        rawHeat = popV{currCell}.(touchOrderFields{g}).(variableFields{1}).spikes(toKeep,:);
        smoothedHeat = imgaussfilt(popV{currCell}.(touchOrderFields{g}).(variableFields{1}).spikes(toKeep,:),gaussFilt,'padding','replicate');
        heatToPlot = smoothedHeat;
        
        subplot(plotrow,plotcolumn,d);
        imagesc(heatToPlot);
        colormap(gca,parula);
        set(gca,'Ydir','normal','ytick',1:sum(toKeep),'yticklabel',[popV{currCell}.(touchOrderFields{g}).(variableFields{1}).range(toKeep)],...
            'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
        for k=1:size(popV{currCell}.(touchOrderFields{g}).(variableFields{1}).counts(toKeep),1)
            text(20,k,num2str(popV{currCell}.(touchOrderFields{g}).(variableFields{1}).counts(k)),'FontSize',8,'Color','white')
        end
        hold on;plot([sum(window<0) sum(window<0)],[length(popV{currCell}.(touchOrderFields{g}).(variableFields{1}).range) 0],'w:')
        axis square
        caxis([0 prctile(heatToPlot(:),99)]) %99th percentile cut off so that outliers don't alias signal
        sampledRange(d) = max([popV{currCell}.(touchOrderFields{g}).(variableFields{1}).range(toKeep)]) - min([popV{currCell}.(touchOrderFields{g}).(variableFields{1}).range(toKeep)]);
        
        selIdx = [find(window==array.meta.responseWindow(1)):find(window==array.meta.responseWindow(2))];
        hold on; plot([selIdx(1) selIdx(1)],[0 30],'-.w')
        hold on; plot([selIdx(end) selIdx(end)],[0 30],'-.w')
    end
    suptitle(touchOrderFields{g})
end