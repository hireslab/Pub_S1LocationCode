%Load whisking and neural time series struct
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
U = struct_editor(U);

hilbertVar = 'angle';
% tStruct = object_location_quantification(U,touchCells,hilbertVar,'off');
tStruct = object_location_quantification(U,1:length(U),hilbertVar,'off');

loc_units = cellfun(@(x) x.is_tuned==1,tStruct);
%%
dk_corr = zeros(1,numel(U));
dt_corr = zeros(1,numel(U));
ccoeff = cell(1,numel(U));
fr_diff = cell(1,numel(U)); 
dk_peak_diff = nan(numel(U),2);

num_curvature_bins = 2; 
num_angle_bins = 8;

for i = 86
    
    % raw variables
    touchOn = [find(U{i}.S_ctk(9,:,:)==1)  ;find(U{i}.S_ctk(12,:,:)==1)];
    touchOff = [find(U{i}.S_ctk(10,:,:)==1)  ;find(U{i}.S_ctk(13,:,:)==1)];
    spikes = squeeze(U{i}.R_ntk(:,:,:));
    angle = squeeze(U{i}.S_ctk(1,:,:));
    touch_dkap = squeeze(U{i}.S_ctk(19,:,:));
    touch_dtheta = squeeze(U{i}.S_ctk(18,:,:));
    touch_tnum = ceil(touchOn./U{i}.t);
    poles = normalize_var(U{i}.meta.motorPosition,10,0);
    touch_poles = poles(touch_tnum);
    
    % during/at touch features (dkappa, dtheta, theta)
    dk = zeros(length(touchOn),1);
    dt = zeros(length(touchOn),1);
    for g = 1:length(touchOn)
        dk_tps = touch_dkap(touchOn(g):touchOff(g));
        [~,dk_idx] = max(abs(dk_tps));
        dk(g) = dk_tps(dk_idx);
        
        dt_tps = touch_dtheta(touchOn(g):touchOff(g));
        [~,dt_idx] = max(abs(dt_tps));
        dt(g) = dt_tps(dt_idx);
    end
    touchTheta = angle(touchOn);
    
    % touch response
    try
        % response calculation
        response_window = U{i}.meta.touchProperties.responseWindow(1):U{i}.meta.touchProperties.responseWindow(2);
        tresponse = mean(spikes(touchOn + response_window),2) * 1000;
        cmap_response = spikes(touchOn + [-25:50]) * 1000; 
        
        
        % protraction dk/dt correlation w/ pole location
        pt = dk<0;
        dk_corr(i) = corr(touch_poles(pt)',dk(pt));
        dt_corr(i) = corr(touch_poles(pt)',dt(pt));
        
        % filtering out retraction touches
        pt_spikes = tresponse(pt);
        pt_dk = dk(pt);
        pt_theta = touchTheta(pt);
        pt_cmap = cmap_response(pt,:);

        % plot 1 - 
        norm_pt = normalize_var(pt_spikes,1,0.1);
        [sort_colors,sort_idx] = sort(norm_pt,'descend');
        figure(03);clf
        subplot(3,2,[1 2])
        scatter(pt_theta(sort_idx),pt_dk(sort_idx),20,repmat(sort_colors,1,3),'filled');
        xlabel('whisker angle at touch')
        ylabel('max curvature')
        title('protraction touch analysis')
        
        %% plot 2 - stratification analyses (low and high curvature      
%         [dk_sorted, dksby] = binslin(pt_dk,[pt_theta pt_spikes pt_cmap],'equalN',num_curvature_bins);
%         
%         sorted = cell(1,numel(dksby));
%         subplot(3,2,[3 4])
%         for p = 1:numel(dksby)
%             [sorted{p}, sortedby] = binslin(dk_sorted{p}(:,1),dk_sorted{p}(:,2),'equalN',num_angle_bins,min(pt_theta),max(pt_theta));     
%             hold on; plot(cellfun(@median, sortedby),cellfun(@mean,sorted{p}))
%             
%             vals = cellfun(@mean,sorted{p});
%             [peak_fr, max_idx] = max(vals);
%             avg_not_max = mean(vals(setdiff(1:length(vals),max_idx))); % mean of not max
%             dk_peak_diff(i,p) = peak_fr - avg_not_max;
%         end
%         
%         avg_curves = cell2mat(cellfun(@(x) cellfun(@mean, x),sorted,'uniformoutput',0));
%         normed = normalize_var(avg_curves(:),0,1);
%         norm_avg_curves = reshape(normed,num_angle_bins,2);
%         ccoeff{i} = corr(avg_curves);
%         fr_diff{i} = norm_avg_curves(:,1)-norm_avg_curves(:,2);
%         
%         xlabel('whisker angle at touch')
%         ylabel('firing rate (spks/s)')
%         legend('high dk','low dk')
%         
%         peak_val = zeros(1,2); 
%         for p = 1:numel(dksby)
%             [sorted{p}, sortedby] = binslin(dk_sorted{p}(:,1),dk_sorted{p}(:,3:end),'equalN',num_angle_bins,min(pt_theta),max(pt_theta));     
%             subplot(3,2,4+p)
%             strat_cmap = cell2mat(cellfun(@mean,sorted{p},'uniformoutput',0));
%             peak_val(p) = max(strat_cmap(:)); 
%             hold on; imagesc(strat_cmap)
%             set(gca,'ytick',1:length(sortedby),'yticklabel',round(cellfun(@mean, sortedby)),...
%             'ylim',[1 length(sortedby)], 'xtick',[0:25:75],'xticklabel',[-25:25:50])
%             colorbar;
%         end
%         
%         subplot(3,2,5); caxis([0 max(peak_val)])
%         subplot(3,2,6); caxis([0 max(peak_val)])
%         
%         suptitle(['cell num ' num2str(i)])
%% plot 3 - stratification analyses - location first then curvature change    

        [dk_sorted, thetasby] = binslin(pt_theta,[pt_dk pt_spikes pt_cmap],'equalN',num_angle_bins,min(pt_theta),max(pt_theta));
        
        sorted = cell(1,numel(thetasby));
        sortedmap = cell(1,numel(thetasby));
       
        for p = 1:numel(thetasby)
            [sorted{p}, sortedby] = binslin(dk_sorted{p}(:,1),dk_sorted{p}(:,2),'equalN',num_curvature_bins);
            [sortedmap{p}, sortedby] = binslin(dk_sorted{p}(:,1),dk_sorted{p}(:,3:end),'equalN',num_curvature_bins);
        end
        % plot of tuning curves
        responses = cell2mat(cellfun(@(x) cellfun(@mean, x) ,sorted,'uniformoutput',0))';
        stim =  cellfun(@median,thetasby);
        subplot(3,2,[3 4])
        plot(stim,responses)
        xlabel('whisker angle at touch')
        ylabel('firing rate (spks/s)')
        legend('high dk','low dk')
        
        % plot heatmaps for different dk
        maps = cell(1,num_curvature_bins);
        for b = 1:numel(maps)
            maps{b} = cell2mat(cellfun(@(x) mean(x{b}), sortedmap,'uniformoutput',0)');
            subplot(3,2,4+b)
            imagesc(imgaussfilt(maps{b},.5))
            set(gca,'ytick',1:length(stim),'yticklabel',round(stim),...
            'ylim',[1 length(stim)], 'xtick',[0:25:75],'xticklabel',[-25:25:50]...
            ,'ydir','normal')
            colorbar;
            
        end
        suptitle(['cell num ' num2str(i)])
        if mean( cellfun(@numel,sorted{1})) < 20
            disp('warning: number of samples per angle bin < 20')
        end
        
        normed = normalize_var(responses(:),0,1);
        norm_avg_curves = reshape(normed,num_angle_bins,2);
        ccoeff{i} = corr(responses);
        fr_diff{i} = norm_avg_curves(:,1)-norm_avg_curves(:,2);
        
        
    catch
        touch_cells = cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),U);
        rw = mean(cell2mat(cellfun(@(x) x.meta.touchProperties.responseWindow,U(touch_cells),'uniformoutput',0)'));
        response_window = round(rw(1):rw(2));
        
        disp(['cell ' num2str(i) ' is not a touch cell. Skipping']) 
    end
    
    pause 
%   if i ==29  
      saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
      fn = ['angle_X_curvature_ex2.eps'];
      export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
      fix_eps_fonts([saveDir, fn])
%   end
end
%% analyses
% peak vs non peak 
loc_diff = dk_peak_diff(loc_units,:);
c_diff = loc_diff(~isnan(loc_diff(:,1)),:);
norm_c_diff = c_diff./max(c_diff,[],2);
figure(23040);clf
plot(norm_c_diff','o-', 'color',[.8 .8 .8])
hold on; errorbar([1 2],mean(norm_c_diff),std(norm_c_diff),'bo-')
set(gca,'xlim',[0.5 2.5])


%analyzing average correlation between high and low dk
nonempty = ccoeff(~cellfun(@isempty, ccoeff));
c_accum = zeros(size(nonempty{1}));
for i = 1:numel(nonempty)
    c_accum = c_accum + nonempty{i};
end

%analyzing diff in normed firing rate between high and low dk.  
figure(230);clf
avg_diff = mean(cell2mat(fr_diff),2);
avg_sd = std(cell2mat(fr_diff),[],2);
avg_sem = avg_sd./sqrt(length(cell2mat(fr_diff)));

shadedErrorBar(linspace(-1,1,numel(avg_diff)),avg_diff*100,avg_sem*100,'k')
hold on; plot([-1 1],[0 0],'--k')
set(gca,'ylim',[-10 30],'xlim',[-1.1 1.1],'xtick',[-1 0 1])
xlabel('normalized pole location')
ylabel('% change in firing rate (high dk - low dk)')
title('on average, higher dk leads to increased firing rates')

% testing to see if rate of change from one position is sig diff from another
anova1(cell2mat(fr_diff)')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
fn = ['angle_X_curvature_diff.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% simple glm for touch cells only
fitCoeffs = cell(1,numel(U));
fitDevExplained = zeros(1,numel(U));
rsq = zeros(1,numel(U));

    % as of 200522, model is garbage at predicting spikes using poisson and
    % dk and theta at touch for each touch. Response probabilites are just
    % way too low. Model using mean predicts FR suggesting that features
    % arent' very useful in driving neural activity over touch response
    % window. 
%     if isfield(U{i}.meta.touchProperties,'responseWindow')
%         
%         glmnetOpt = glmnetSet;
%         glmnetOpt.standardize = 0;
%         glmnetOpt.alpha = 0.95;
%         glmnetOpt.xfoldCV = 5;
%         glmnetOpt.numIterations = 25;
%         
%         DmatX = [pt_dk pt_theta];
%         norm_DmatX = (DmatX - mean(DmatX)) ./ std(DmatX);
%         DmatY = pt_spikes;
%         
%         % test and train split
%         exampleIdx = 1:length(pt_spikes);
%         shuffIdx = exampleIdx(randperm(length(exampleIdx)));
%         trainIdx = shuffIdx(1:round(length(pt_spikes)*.8));
%         testIdx = setdiff(shuffIdx,trainIdx);
%         
%         trainDmatX = norm_DmatX(trainIdx',:);
%         trainDmatY = DmatY(trainIdx',:);
%         testDmatX = norm_DmatX(testIdx',:);
%         testDmatY = DmatY(testIdx',:);
%         
%         % model fitting
%         try
%             cv = cvglmnet(trainDmatX,trainDmatY,'poisson',glmnetOpt,[],glmnetOpt.xfoldCV);
%             
%             % model evaluation
%             fitLambda = cv.lambda_1se;
%             iLambda = find(cv.lambda == fitLambda);
%             fitCoeffs{i} = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
%             
%             model = exp([ones(length(testDmatX),1),testDmatX]*fitCoeffs{i});
%             mu = mean(testDmatY); % null poisson parameter
%             nullLogLikelihood = sum(log(poisspdf(testDmatY,mu)));
%             saturatedLogLikelihood = sum(log(poisspdf(testDmatY,testDmatY)));
%             fullLogLikelihood = sum(log(poisspdf(testDmatY,model)));
%             fitDevExplained(i) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
%             MSE = mean((testDmatY - model).^2);
%             MSE_bl = mean((testDmatY-mu).^2);
%             rsq(i) = 1 - (MSE/MSE_bl);
%         catch
%             fitDevExplained(i) = nan;
%             fitCoeffs{i} = nan(size(DmatX,2)+1,1);
%             rsq(i) = nan;
%         end
%     else
%         fitDevExplained(i) = nan;
%         fitCoeffs{i} = nan(size(DmatX,2)+1,1);
%         rsq(i) = nan;
%     end




