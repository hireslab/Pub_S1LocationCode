%this function will output the times of touches, time of touches
%for go, time of touches nogo, time in between touches predecision, and the
%length of each of those touches, #mean number of touches. All outputs are in a cell array with each
%cell being 1 trial.


function [Gprelixtimes,GprelixtimesGO,GprelixtimesNOGO,Gprelixpcth,Gprelixduration,touchespredecision] = predecision_variables(array)

Gprelixpcth=cell(1,array.k);
Gprelixtimes=cell(1,array.k);
Gprelixduration=cell(1,array.k);

licks = squeeze(array.S_ctk(16,:,:))';

fl = [];
for i = 1:size(licks,1)
    fl(i) = min([find(licks(i,:)==1) array.t]);
end
fl = fl(fl<3000);
fls = sort(fl);
flninety= fls(ceil(numel(fls)*.9));

for k =1:array.k %for each trial
    postsamp=array.meta.poleOnset(k)+.95; %pole onset time + time to be touchable (.2s) + sampling period (.75s)
    lixtmp=find(array.S_ctk(16,:,k)==1);
    lix=lixtmp(find(lixtmp>=postsamp*1000,1,'first')); %find first lick after sampling period
    alltouches=[find(array.S_ctk(9,:,k)==1), find(array.S_ctk(12,:,k)==1)]; %all touches in trial
    
    if lix>0 %all touches before first decision
        prelix=alltouches(alltouches<lix); 
    else
        prelix=alltouches(alltouches<flninety); %setting touch cut off to be 90 percentile of all lix
    end
    
    duration = [find(array.S_ctk(10,:,k)==1) , find(array.S_ctk(13,:,k)==1)] - alltouches;
    
    pcth=zeros(1,length(alltouches));
    for i=2:length(alltouches)+1
        tmp=[0 alltouches];
        pcth(i-1)=tmp(i)-tmp(i-1);
    end
    
    if numel(lix)>0 %looking at the period before lick. If CR, all touches in that trial are counted here.
        Gprelixpcth{k}=pcth(alltouches<lix);
        Gprelixduration{k}=duration(alltouches<lix);
    else
        Gprelixpcth{k}=pcth;
        Gprelixduration{k}=duration;
    end
    
    if isempty(Gprelixduration{k}) %padding trials with no prelick touches with 0 duration
        Gprelixduration{k}=0;
    end
        
    if array.meta.trialCorrect(k)==1 %adding a row below all durations for trial correct or not
        %Gprelixduration{k}=[Gprelixduration{k}(2:end);ones(1,length(Gprelixduration{k}(2:end)))];
        Gprelixduration{k}=[Gprelixduration{k};ones(1,length(Gprelixduration{k}))];
    else
        %Gprelixduration{k}=[Gprelixduration{k}(2:end);zeros(1,length(Gprelixduration{k}(2:end)))];
        Gprelixduration{k}=[Gprelixduration{k};zeros(1,length(Gprelixduration{k}))];
    end
    
    if numel(Gprelixpcth{k})>1 %selecting only those trials with touches > 1 in the prelick area
        Gprelixpcth{k}=Gprelixpcth{k}(2:end);
        %Gprelixduration{k}=[Gprelixduration{k}(2:end)] old version where
        %we just output duration 
    else
        Gprelixpcth{k}=[];
        %Gprelixduration{k}=[];
    end
    
    Gprelixtimes{k}=prelix;
    
end

%find avg number of touches predecision
touchespredecision=mean(cellfun(@numel,Gprelixtimes));

%sort pre lick times by GO and NOGO trials
GprelixtimesGO = Gprelixtimes(logical(array.meta.trialType));
GprelixtimesNOGO = Gprelixtimes(logical(1-array.meta.trialType));
Gprelixpcth = Gprelixpcth; %cell array w/ time between touches
Gprelixduration = Gprelixduration; %cell array w/ time for each touch
