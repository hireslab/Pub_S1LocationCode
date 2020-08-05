%%
clear U
layer = 'BVL5b';
cellNum = [1 29:38 41:59];
trialCutoffs = repmat([1 250],numel(cellNum),1);
trialCutoffs(10,:)=[1 95];
%%
clear U
layer = 'NL5b';
cellNum = 1:30;
cellNum = [1:7 9:11 13:29 31]; %some cells selected out because they are interneurons
% trialCutoffs = repmat([1 500],numel(cellNum),1);
%%
clear U
layer = 'L5bOut';
cellNum = [1:5 7];
trialCutoffs = repmat([1 200],numel(cellNum),1);
%%
clear U
layer = 'L5bInt';
cellNum = [1:6];
trialCutoffs = repmat([1 500],numel(cellNum),1);
%% 
clear U
layer = 'Phil';
cd(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\Projects\Characterization\' layer '\TArrays'])
files = dir('*.mat');
for y = 1:length(files)
    fns_trials(y) = str2num(files(y).name(13:16));
end

cd(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\Projects\Characterization\' layer '\Contacts'])
files = dir('*.mat');
for y = 1:length(files)
    fns_contacts(y) = str2num(files(y).name(7:10));
end

cellNum = intersect(fns_contacts,fns_trials);
cellNum = cellNum(cellNum ~= 2075); %tmp removal because double cell and contact/T array lengths are different
trialCutoffs = repmat([1 999],numel(cellNum),1);

%%

timeIdx_threshold = 4000; %this is the number of time points that must be present to use trial; 

for cellStep = 1:length(cellNum)
    
    cd(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\Projects\Characterization\' layer '\TArrays'])
    load(['trial_array_' num2str(cellNum(cellStep)) '.mat'])
    cd(['C:\Users\jacheung\Dropbox\HLabBackup\Jon\Projects\Characterization\' layer '\Contacts'])
    load(['ConTA_' num2str(cellNum(cellStep)) '.mat'])
    
    d.varNames = {'thetaAtBase', 'velocity', 'amplitude', 'setpoint', 'phase', ...
        'deltaKappa','M0Adj','FaxialAdj', 'firstTouchOnset', 'firstTouchOffset', ...
        'firstTouchAll', 'lateTouchOnset','lateTouchOffset','lateTouchAll','PoleAvailable','beamBreakTimes',...
        'kappa'};
    
    d.cellNum = cellNum(cellStep);
    
    d.t = max(cellfun(@(x)round(x.whiskerTrial.time{1}(end)*1000)+1,T.trials(T.whiskerTrialInds)));
    
    [~,useTrials] = intersect(T.trialNums,T.whiskerTrialNums);
    useTrials = useTrials(useTrials >= trialCutoffs(cellStep,1) & useTrials <= trialCutoffs(cellStep,2));
    
    %Trial checking to ensure both T array and contacts have same number of
    %timepoints JC190711
    timeStamps = zeros(1,length(useTrials));
    timeStampsContacts = zeros(1,length(useTrials));
    for i = 1:length(useTrials)
        timeStamps(i) = length(round(T.trials{useTrials(i)}.whiskerTrial.time{1}*1000)+1);
        timeStampsContacts(i) = length(contacts{useTrials(i)}.M0combo{1});
    end
    useTrials = useTrials(timeStamps == timeStampsContacts);
    
    d.k = length(useTrials);
    d.u = 1;
    d.c = 16;
    d.S_ctk = nan(d.c, d.t, d.k);
    d.R_ntk = zeros(1, d.t, d.k);
    
    
    useTrials = useTrials';
    traj = 1;
    
    
    for i = 1:length(useTrials)
        display(i)
        timeIdx = round(T.trials{useTrials(i)}.whiskerTrial.time{traj}*1000)+1;
        
        
        theta = nan(1,d.t);
        theta(timeIdx) = T.trials{useTrials(i)}.whiskerTrial.thetaAtBase{traj};
        
        nanidx = find(isnan(theta));
        nanidx = nanidx(nanidx > 2 & nanidx < length(timeIdx)-2);
        for j = nanidx
            if j>=3999
                theta(j)=nanmean(theta(j+[-2:0]));
            else
                theta(j) = nanmean(theta(j+[-2:2]));
            end
        end
        
        firstTouchOn    = [];
        firstTouchOff   = [];
        firstTouchAll   = [];
        
        lateTouchOn     = [];
        lateTouchOff    = [];
        lateTouchAll    = [];
        
        %FROM PHIL TO IMPLEMENT AUTOCURATOR CONTACTS
        numContactsTest = contacts{useTrials(i)}.contactInds{1};
        if isstr(numContactsTest) ||  isempty(numContactsTest) %if auto curator marked it as a non touch in preliminary sweep or if no touches
            segmentsTMP = [];
        else
            diffDetect = numContactsTest(:) - (1:length(numContactsTest))';
            [startNums , uniqueInds,~] = unique(diffDetect);
            startNums = numContactsTest(uniqueInds);
            lengthArray = circshift(uniqueInds, -1) - uniqueInds;
            lengthArray(end) = 1+length(numContactsTest)-uniqueInds(end);
            segmentsTMP = [startNums(:) , startNums(:)+lengthArray(:)-1];
            
        end
        
        contacts{useTrials(i)}.segmentInds{traj} =segmentsTMP;
        
        
        if isfield(contacts{useTrials(i)},'segmentInds')
            if ~isempty(contacts{useTrials(i)}.segmentInds{traj})
                firstTouchOn  = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(1,1)));
                firstTouchOff = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(1,2)));
                firstTouchAll = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(1,1):contacts{useTrials(i)}.segmentInds{traj}(1,2)));
                
                if size(contacts{useTrials(i)}.segmentInds{traj},1)>1
                    lateTouchOn  = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(2:end,1)));
                    lateTouchOff = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(2:end,2)));
                    
                    for j = 2:size(contacts{useTrials(i)}.segmentInds{traj},1)
                        lateTouchAll = cat(2,lateTouchAll,round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(j,1):contacts{useTrials(i)}.segmentInds{traj}(j,2))));
                    end
                    
                end
            end
        end
        
        [hh amplitude  filteredSignal setpoint amplitudeS setpointS phase phaseS] =  SAHWhiskerDecomposition(theta);
        
        travelIn = 0;
        travelOut = 250;
        
        vel = diff([0 theta])/.001;
        
        pinIn = round(1000*T.trials{useTrials(i)}.pinDescentOnsetTime + travelIn);
        pinOut = min([d.t round(1000*T.trials{useTrials(i)}.pinAscentOnsetTime + travelOut)]);
        
        d.S_ctk(1,:,i) = theta;
        d.S_ctk(2,:,i) = vel;
        d.S_ctk(3,:,i) = amplitude;
        d.S_ctk(4,:,i) = setpoint;
        d.S_ctk(5,:,i) = phase;
        d.S_ctk(6,timeIdx,i) = T.trials{useTrials(i)}.whiskerTrial.deltaKappa{traj};
        d.S_ctk(7,timeIdx,i) = contacts{useTrials(i)}.M0comboAdj{traj};
        d.S_ctk(8,timeIdx,i) = contacts{useTrials(i)}.FaxialAdj{traj};
        d.S_ctk(9,firstTouchOn,i) = 1;
        d.S_ctk(10,firstTouchOff,i) = 1;
        d.S_ctk(11,firstTouchAll,i) = 1;
        d.S_ctk(12,lateTouchOn,i) = 1;
        d.S_ctk(13,lateTouchOff,i) = 1;
        d.S_ctk(14,lateTouchAll,i) = 1;
        d.S_ctk(15,:,i) = 0;
        d.S_ctk(15,pinIn:pinOut,i) = 1;
        d.S_ctk(16, ceil(1000*T.trials{useTrials(i)}.beamBreakTimes(T.trials{useTrials(i)}.beamBreakTimes > 0 & T.trials{useTrials(i)}.beamBreakTimes < d.t/1000)),i) = 1;
        d.S_ctk(17,timeIdx,i) = T.trials{useTrials(i)}.whiskerTrial.meanKappa{traj};
        
        d.S_ctk(6,setdiff(1:d.t,contacts{useTrials(i)}.contactInds{traj}),i) = 0;
        d.S_ctk(7,setdiff(1:d.t,contacts{useTrials(i)}.contactInds{traj}),i) = 0;
        d.S_ctk(8,setdiff(1:d.t,contacts{useTrials(i)}.contactInds{traj}),i) = 0;
        
        
        
        
        
    end
    
    %end
    
    
    d.S_ctk(15,1:min(find(nansum(squeeze(d.S_ctk(9,:,:))')))-1,:)= 0;  % define pole availablitity onset
    for i =1:length(useTrials);
        
        spiketimes = round(T.trials{useTrials(i)}.spikesTrial.spikeTimes/10-T.whiskerTrialTimeOffset*1000);
        spiketimes = spiketimes(spiketimes>0 & spiketimes <=d.t);
        d.R_ntk(1,spiketimes,i) = 1;
    end
    
    U{cellStep} = d;
    
    U{cellStep}.meta.mouseName         = T.mouseName;
    U{cellStep}.meta.sessionName       = T.sessionName;
    U{cellStep}.meta.cellName          = T.cellNum;
    U{cellStep}.meta.cellCode          = T.cellCode;
    U{cellStep}.meta.usedTrialNums        = useTrials;
    U{cellStep}.meta.depth             = T.depth;
    U{cellStep}.meta.layer             = layer;
    
    U{cellStep}.meta.motorPosition     = cellfun(@(x)x.behavTrial.motorPosition,T.trials(useTrials));
    U{cellStep}.meta.goPosition        = cellfun(@(x)x.behavTrial.goPosition,T.trials(useTrials));
    U{cellStep}.meta.nogoPosition      = cellfun(@(x)x.behavTrial.nogoPosition,T.trials(useTrials));
    U{cellStep}.meta.trialType         = cellfun(@(x)x.behavTrial.trialType,T.trials(useTrials));
    U{cellStep}.meta.trialCorrect      = cellfun(@(x)x.behavTrial.trialCorrect,T.trials(useTrials));
    U{cellStep}.meta.poleOnset         = cellfun(@(x)x.behavTrial.pinDescentOnsetTime,T.trials(useTrials)) - T.whiskerTrialTimeOffset;
    U{cellStep}.meta.poleOffset        = cellfun(@(x)x.behavTrial.pinAscentOnsetTime,T.trials(useTrials)) - T.whiskerTrialTimeOffset;


    U{cellStep}.meta.ranges            = unique([cellfun(@(x)x.behavTrial.nogoPosition,T.trials(useTrials)),cellfun(@(x)x.behavTrial.goPosition,T.trials(useTrials))]);
    % U{cellStep}.meta.C2distance      = sqrt(SU.recordingLocation{cellNum(cellStep)}(1)^2+SU.recordingLocation{cellNum(cellStep)}(2)^2)
    %  U{cellStep}.meta.isC2           = SU.distance{cellNum(cellStep)} < .16;

    U{cellStep}.meta.stimTrials        = {};
    
    U{cellStep}.whisker.follicleX      = cellfun(@(x)x.whiskerTrial.follicleCoordsX,T.trials(useTrials));
    U{cellStep}.whisker.follicleY      = cellfun(@(x)x.whiskerTrial.follicleCoordsY,T.trials(useTrials));
    U{cellStep}.whisker.barPos         = cellfun(@(x)x.whiskerTrial.barPos,T.trials(useTrials),'UniformOutput',false);
end

cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')

