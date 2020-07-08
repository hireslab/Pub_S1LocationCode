% mouseName = 'AH0716'; %main example
% sessionName ='170902';
% videoloc = 'JON';


mouseName = 'AH0761'; 
sessionName ='171102';
videoloc = 'JON';
d= (['E:\' videoloc filesep mouseName filesep sessionName filesep])

cd(d)

%%
%farthest trial, example trial, closest trial
% trialnums ={'14','94','110'}; main example
trialnums = {'25','56','124'}

rc = numSubplots(length(trialnums));
figure(4502);clf
for i = 1:length(trialnums)
    trialnum = trialnums{i};
WSTName = [mouseName 'x' sessionName '-' trialnum '_' 'WST' '.mat'];
% WSTName = [mouseName 'x170901-' trialnum '_' 'WST' '.mat'];
load(WSTName)%load file based on trial above to test mask 
tp = [.5  1.25 5]; %last value is downsample rate
%%%%%% plot any mask you want use trial number above  
subplot(rc(1),rc(2),i);
ws.plot_fitted_whisker_time_projection(0,'k',tp)
hold on; 
ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
ws.plot_mask(0,'g',tp);
ws.plot_follicle_position_time_projection(0,'b.',tp)
 axis square 
 tmp = load([mouseName 'x' sessionName '-' trialnum '.bar']);
% tmp = load([mouseName 'x170901-' trialnum '.bar']);
s=scatter(tmp(1,2),tmp(1,3),'filled','c');
s.SizeData = 200;

set(gca,'visible','off','xlim',[0 350],'ylim',[0 400])

end
% 
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
%     fn = 'pole_distance.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])

%% 73,103