function gof = decoderPerformance(mdl)

confusionMatrix = confusionmat(cell2mat(mdl.io.trueY),cell2mat(mdl.io.predY));
chance = 1/length(confusionMatrix);

accuracy = nan(1,length(mdl.io.trueY)); 
for g = 1:length(mdl.io.trueY)
    confusionMatrix_accuracy = confusionmat(mdl.io.trueY{g},mdl.io.predY{g});
    accuracy(g) = sum(diag(confusionMatrix_accuracy)) ./ sum(confusionMatrix_accuracy(:));
end

figure(10);clf;
subplot(1,2,1)
predProb = confusionMatrix ./ sum(confusionMatrix,2); %sum(cmat,1) for precision and sum(Cmat,2) for recall/sensitivity
imagesc(predProb)
caxis([0 prctile(predProb(:),99)])
caxis([0 .40])
set(gca,'xtick',[],'ytick',[])
xlabel('predicted');ylabel('true')
colorbar
colormap turbo
axis square


distance_from_true = cellfun(@(x,y) abs(x-y),mdl.io.trueY,mdl.io.predY,'uniformoutput',0);

distances = 0:max(cell2mat(distance_from_true));
probs = nan(numel(distance_from_true),numel(distances));
for i = 1:length(distances)
    check_dist = distances(i);
    probs(:,i) = cellfun(@(x) mean(x<=check_dist),distance_from_true);
end

meanProb = mean(probs);
SEM = nanstd(probs)./ sqrt(size(probs,1));
ts = tinv(.95,size(probs,1));
CI = SEM.*ts;

subplot(1,2,2);
shadedErrorBar(distances,meanProb,std(probs),'k')
set(gca,'ylim',[0 1],'xlim',[0 8],'xtick',0:2:10,'xticklabel',0:.5:5) %hard coded xticklabels for single touch prediction of pole position using population of OL tuned cells
xlabel('distance from prediction (mm)');ylabel('p(prediction)')
axis square

gof.accuracy = accuracy;
gof.cmat_raw = confusionMatrix;
gof.cmat = predProb;
gof.resolutionNames = {'distances','mean probability','std probability'};
gof.resolution = [distances' meanProb' std(probs)'];

% hold on;
% binedges = max(mdl.io.Y.normal)-.5;
% plotedges = max(mdl.io.Y.normal)-1;
% n = histcounts(mdl.io.trueXpredicted(:,2)-mdl.io.trueXpredicted(:,1),[-binedges:1:binedges]);
% x = -plotedges:1:plotedges;
% y = n/sum(n);
% plot(x,y,'k');
% hold on; plot(x,ones(1,length(x)).*chance,'-.k')
% axis square
% set(gca,'ylim',[0 .4])
% xlabel('distance of prediction from true');
% ylabel('proportion of trials')
% xdata = y(find(x==0):end);
% ydata = x(find(x==0):end);
% resolution = spline(xdata(1:9),ydata(1:9),chance);
%
%



% DEPRECATED 190524 since no longer doing angles
% sampledSpanPerCell = cellfun(@(x) (max(x.allTouches.theta.range) - min(x.allTouches.theta.range)),selectedArray);
% binResolutionMean = (mean(sampledSpanPerCell)./size(predProb,1));
% binResolutionStDev = (std(sampledSpanPerCell)./size(predProb,1));
%
% resolution = nan(size(predProb,1),1);
% for i = 1:size(predProb,1)
%     if ~isempty(find(predProb(i,:)<chance))
%     resolution(i) = min(abs(i - find(predProb(i,:)<chance)));
%     end
% end
% statsResolution = [mean(resolution.*binResolutionMean) std(resolution.*binResolutionMean)];
% subplot(1,2,1)
% title(['chance = ' num2str(chance) '. Resolution = ' num2str(statsResolution(1)) '+/-' num2str(statsResolution(2))])
%
%
