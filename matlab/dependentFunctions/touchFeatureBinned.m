%this function will be used to build a struct necessary for feature
%tuning quantification

function popV = touchFeatureBinned(U,viewWindow)

popV = cell(length(U),1);
preDecisionTouches = preDecisionTouchMat(U);
    
for rec=1:length(U)
    countThresh = 5; %min touch in each bin to be considered

    [tVar] = atTouch_sorter(U{rec},viewWindow,preDecisionTouches{rec});
    
    % Run this line below if you want to look at ret touches.
    %varspikes(varspikes(:,5)<0,:)=[];
    
    touchOrder = fields(tVar);
%     touchOrder = {'allTouches'};
    
    for g = 1:numel(touchOrder)
        %%
        avail_fields = tVar.(touchOrder{g}).varNames;
        bounds = [{[-100:2:100]} {[-9750:125:9750]} {[-100:2:100]}  {[-100:2:100]}  {linspace(pi*-1,pi,12)} {[-.95:.05:.95]}];
        for d = [1 2 3 4 5 6] %for variables 1:6
            if d == 2
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).S_ctk(:,d),tVar.allTouches.R_ntk,'equalE',numel(bounds{d})+1,-10000,10000);
            elseif d == 5
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).S_ctk(:,d),tVar.allTouches.R_ntk,'equalX',numel(bounds{d})+1);
            elseif d == 6
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).S_ctk(:,d),tVar.allTouches.R_ntk,'equalE',numel(bounds{d})+1,-1,1);
            else
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).S_ctk(:,d),tVar.allTouches.R_ntk,'equalE',numel(bounds{d})+1,-100,100);
            end
            
            binrange = bounds{d};
            
            % Trimming unused ranges
            trims=[binrange' cell2mat(cellfun(@mean,sortedBy,'Uniformoutput',0))];
            indexHasTheta = ~isnan(trims(:,2));
            trims = trims(indexHasTheta, :);
            counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
            
            % Populating fields with values
            V.(avail_fields{d}).bounds = bounds{d}; 
            V.(avail_fields{d}).counts =counttmp(indexHasTheta==1,1);
            V.(avail_fields{d}).range = trims(:,1);
            V.(avail_fields{d}).spikes=cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
            V.(avail_fields{d}).spikes=V.(avail_fields{d}).spikes(indexHasTheta==1,:);
            V.(avail_fields{d}).raw = sorted(indexHasTheta==1);
            
            %Trimming bins with touchcounts below thresholds
            mintouches=find(V.(avail_fields{d}).counts<countThresh);
            V.(avail_fields{d}).counts(mintouches,:)=[];
            V.(avail_fields{d}).spikes(mintouches,:)=[];
            V.(avail_fields{d}).range(mintouches,:)=[];
            V.(avail_fields{d}).raw(mintouches,:) = [];
            
            popV{rec}.(touchOrder{g})=V;
        end    
        
    end
end
