function viewRaster(array)

%optional raster of FRs for tuned cells. 
figure(8);clf
allSpks = squeeze(array.R_ntk);
[~,idx] = sort(array.meta.motorPosition);
allSpks = allSpks(:,idx);
for k = 1:size(allSpks,2)
    st = find(allSpks(:,k)==1);
    if ~isempty(st)
        figure(8);hold on
        scatter(st,ones(length(st),1).*k,[],'.k')
    end
end
set(gca,'ylim',[1 array.k])
ylabel('sorted trials from near to far')
