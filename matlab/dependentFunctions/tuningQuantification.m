function tuningStruct = tuningQuantification(cellStruct,popV,selectedCells,variableFields,touchOrderFields,viewWindow,displayOpt)

%% ZSCORE CONVERSION and tuning width

for t = 1:length(variableFields)
    bounds = popV{1}.allTouches.(variableFields{t}).bounds;
    
    if (nargin < 7), displayOpt = 'on'; end
    willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
        | strcmp(displayOpt,'off'));
    
    if (willdisplay)
        [rc] = numSubplots(length(selectedCells));
        plotrow = rc(1);
        plotcolumn = rc(2);
    end
    
    tunedCellMat =  nan(length(touchOrderFields),length(popV));
    
    for g = 1:numel(touchOrderFields)
        
        fr= nan(length(selectedCells),2);
        pw=fr;
        pwzs = nan(length(selectedCells),50);
        pwtheta = nan(length(selectedCells),2);
        
        for d = 1:length(selectedCells)
            
            %ZSCORREEEEE
            currCell = selectedCells(d);
            
            counts = popV{currCell}.(touchOrderFields{g}).(variableFields{t}).counts;
            ranges = popV{currCell}.(touchOrderFields{g}).(variableFields{t}).range ;
            thetavect = [];
            for k = 1:length(ranges)
                thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
            end
            
            allspikes = cell2mat(popV{currCell}.(touchOrderFields{g}).(variableFields{t}).raw);
            meanbase = mean(mean(allspikes(:,1:find(viewWindow==0)),2));
            stdbase = std(mean(allspikes(:,1:find(viewWindow==0)),2));
            %     postresponses = allspikes(:,find(window==0)+5:find(window==0)+35);
            postresponses = allspikes(:,find(viewWindow==0)+ cellStruct{currCell}.meta.responseWindow(1) : find(viewWindow==0)+cellStruct{currCell}.meta.responseWindow(2));
            
            xresp = mean(postresponses,2);
            zscore = (xresp - meanbase) ./ stdbase;
            
            if strcmp(variableFields{t},'phase')
                [sortedZscores,~,~] = binslin(thetavect,zscore,'equalE',numel(bounds)+1,bounds(1),bounds(end));
                [sortedRawResponse] = binslin(thetavect,xresp,'equalE',numel(bounds)+1,bounds(1),bounds(end));
            else
                [sortedZscores,~,~] = binslin(thetavect,zscore,'equalE',numel(bounds),bounds(1),bounds(end));
                [sortedRawResponse] = binslin(thetavect,xresp,'equalE',numel(bounds),bounds(1),bounds(end));
            end
            
            
            samps = cellfun(@numel,sortedZscores);
            selBins = samps>(sum(samps)./100);
            
            zraw = transpose(cell2nanmat(sortedZscores));
            rrnanmat = transpose(cell2nanmat(sortedRawResponse));
            
            selBinMat = nan(size(zraw));
            selBinMat(selBins,:) = 1;
            
            zraw = zraw .* selBinMat;
            rraw = rrnanmat .* selBinMat;
            
            rawmeanresp =  nanmean(rraw,2);
            rawmeanresp = rawmeanresp(~isnan(rawmeanresp));
            rawSEM = nanstd(rraw,[],2) ./ sqrt(sum(~isnan(rraw),2));
            rawSEM = rawSEM(~isnan(rawSEM));
            
            %sig calculation
            [pZ,~,stats]=anova1(zraw',[],'off');
   
            vals =multcompare(stats,[],'off');
            x=cellfun(@str2num,stats.gnames);
            
            %CI BINS for zscores
            cibins = nan(size(x,1),1);
            SEMg = nanstd(zraw,[],2) ./ sqrt(sum(~isnan(zraw),2));
            for i = 1:length(x)
                rawx = zraw(x(i),:);
                SEM = SEMg(x(i));
                ts = tinv(.95,sum(~isnan(rawx),2)-1);      % T-Score
                cibins(i,:) = ts.*SEM;   %confidence intervals
            end
            
            
            y = stats.means;
            x= bounds(x)';

            xy = [x y'];

            if willdisplay
                figure(390+g)
                subplot(plotrow,plotcolumn,d);
                shadedErrorBar(xy(:,1),xy(:,2),cibins,'-k');
                set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
                    'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
                axis square
            end
            
            zs = nanmean(zraw,2);
            zs(isnan(zs))=[];
            
            sigs =  vals(vals(:,end)<.01,:);
            
            clear peaksig
            [~,maxidx]=max(zs);
            midxsel = find(sum(sigs(:,[1 2])==maxidx,2)==1);
            
            if sum(midxsel>0)
                peaksig = sigs(midxsel,1:2);
                
                [~,midx] = min(abs(diff(peaksig')));
                pw(d,:) = peaksig(midx,:);
                if willdisplay
                    figure(390+g);
                    hold on; scatter(xy(peaksig(midx,:),1),xy(peaksig(midx,:),2),[],'filled','g')
                end
                
                raws=zs(pw(d,1):pw(d,2));
                pwzs(d,1:length(raws)) = raws;
                pwtheta(d,:) = x(pw(d,:));
            end
            
%             tuningXYerr{selectedCells(d)} = [xy cibins rawmeanresp rawSEM];
%             tuningZraw{selectedCells(d)} = zraw;
            tuningRraw{selectedCells(d)} = rraw;

        end  
        ocellidx = find(~isnan(pw(:,1)));
        tunedCellMat(g,selectedCells(ocellidx)) = 1;
        
        tuningStruct.R_ntk.(touchOrderFields{g}).(variableFields{t}) = tuningRraw;
        tuningStruct.S_ctk.(variableFields{t})  = mean([bounds(1:end-1);bounds(2:end)]); 
    end
    tuningStruct.populationQuant.(variableFields{t}).rowFeatNames = touchOrderFields;
    tuningStruct.populationQuant.(variableFields{t}).matrix = tunedCellMat;
    
    
end

