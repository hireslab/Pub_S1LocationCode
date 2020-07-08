%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
% U = defTouchResponse(U,.95,'off');
% selectedCells = find(cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),U)~=0);
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));

pole_tuned = object_location_quantification(U,selectedCells,'pole','off'); %for old see object_location_v1.0

tuned_cells = find(cellfun(@(x) x.is_tuned,pole_tuned)==1);
% untuned_cells = intersect(find(~(cellfun(@(x) x.is_tuned,pole_tuned)==1)),selectedCells);
% builtUnits = find(cellfun(@(x) isfield(x,'stim_response'),pole_tuned));
sampledSpace = cellfun(@(x) range(x.stim_response.values(:,1)),pole_tuned(tuned_cells)) ./ 2;
units_2_use = tuned_cells(sampledSpace > .8);

%% Justification for decoder - fano factor check
tuned_cells = find(cellfun(@(x) x.is_tuned,pole_tuned)==1);
[ff_bin] = poisson_sampling_justification(U,pole_tuned); %only plots for location tuned units

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
% fn = 'poisson_fano_factor.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn]);


mod_idx_relative = cellfun(@(x) x.calculations.mod_idx_relative,pole_tuned(tuned_cells));
regressed_ff = cellfun(@(x) nanmean(x.fano_factor),ff_bin);

figure(8);clf
scatter(mod_idx_relative,regressed_ff,'filled','k')
axis square
xlabel('location mod. depth');ylabel('"regressed" FF')
lm = fitlm(mod_idx_relative,regressed_ff);
x = (0:.1:1)';
y = lm.predict(x);
hold on; plot(x,y,'r')
title(num2str(lm.Rsquared.ordinary));
axis square
xlabel('location mod. depth');ylabel('"regressed" FF')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
fn = 'poisson_fano_factor_x_mod_depth.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn]);

%% population at touch pole decoding
% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 10;
glmnetOpt.pctSamplingThreshold = .80; %what percent of total pole positions must be sampled before using that unit
glmnetOpt.interpResolution = 40; %10mm / numBins (e.g. 40 = .25mm resolution)
glmnetOpt.samplingOption = 'poisson';
glmnetOpt.numberResamples = 50;

fileName = 'glm_location_decoder';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] =designMatrixBlocks_poleDecoder(glmModel,pole_tuned,glmnetOpt); %current model build is using all touch cells
end

mdlResults = {};

DmatXraw = cell2mat(cellfun(@(x) x.io.DmatX,glmModel,'uniformoutput',0));
mdlResults.io.Xnorm = (DmatXraw - mean(DmatXraw)) ./ std(DmatXraw);
mdlResults.io.Y.normal = glmModel{1}.io.DmatY;

mdlResults = multinomialModel(mdlResults,mdlResults.io.Xnorm,mdlResults.io.Y.normal,glmnetOpt);
%% decoding resolution and probability of guesses
gof = decoderPerformance(mdlResults);

suptitle([ 'Per touch location decoding using ' num2str(size(DmatXraw,2)) ' tuned units'])

%% real psychometric curve 
publication_data_location = 'C:\Users\jacheung\Dropbox\LocalizationBehavior\DataStructs\publication';
dataStructLocation = publication_data_location;
load([dataStructLocation filesep 'behavioral_structure.mat'])
psychometric_curve = rawPsychometricCurves(BV);
interp_psycho = cellfun(@(x) interp1(linspace(1,40,10),x(:,3),1:40),psychometric_curve,'uniformoutput',0);
real_psycho_mean = fliplr(mean(cell2mat(interp_psycho'))); 
real_psycho_sem = fliplr(std(cell2mat(interp_psycho')) ./ sqrt(numel(BV))); 

%% number of neurons for resolution
usedUnits = cellfun(@(x) x.params.cellNum,glmModel);

numNeurons = [1:25];
numIterations = 500;

mdl_mean = median(cell2mat(cellfun(@(x) x(:),mdlResults.fitCoeffs,'uniformoutput',0)),2);
reshaped_coeffs = reshape(mdl_mean,size(mdlResults.fitCoeffs{1}));

nframe = numNeurons;
neurometric_curve = cell(1,length(numNeurons));
% v = VideoWriter('resolution_heatmap.avi');
% v.FrameRate = 1;
% open(v)
resamp_mdl = [];
gofmetrics = []; 
for g = 1:length(numNeurons)
    
    pred = cell(1,numIterations);
    true = cell(1,numIterations);
    for d = 1:numIterations
%           neuron_to_use = datasample(1:numel(usedUnits),numNeurons(g),'Replace',false);
        neuron_to_use = datasample(1:numel(usedUnits),numNeurons(g));
        
        coeffs_to_use = reshaped_coeffs([1 neuron_to_use+1],:);
        dmatX = [ones(length(mdlResults.io.Y.normal),1) mdlResults.io.Xnorm(:,neuron_to_use)];
        dmatY = mdlResults.io.Y.normal;
        
        predicts = dmatX * coeffs_to_use;
        probability =  1 ./ (1+exp(predicts*-1)); %convert to probability by using mean function (for binomial, it is the sigmoid f(x) 1/1+exp(-predicts))
        
        %highest probability = predict class
        [~,pred{d}] = max(probability,[],2);
        true{d} = dmatY;
    end
    
    resamp_mdl{g}.numNeurons = numNeurons(g);
    resamp_mdl{g}.io.trueY = true';
    resamp_mdl{g}.io.predY = pred';
    gofmetrics{g} = decoderPerformance(resamp_mdl{g});
    suptitle(['Number of neurons = ' num2str(numNeurons(g))])
    
    %calculate neurometric curve from matrix
    %how do we do this? create first a confusion matrix
    %
    % pred(Go)    pred(nogo)
    % aka lick    aka no lick
    %%%%%%%%%%%%%%%%%%%%%%%%
    %           |           |
    %    HIT    |    MISS   |  true(Go)
    %           |           |
    %%%%%%%%%%%%|%%%%%%%%%%%|
    %           |           |
    %    FA     |    CR     |  true(Nogo)
    %           |           |
    %%%%%%%%%%%%|%%%%%%%%%%%|
    %
    for b = 1:numIterations
        raw_mat = confusionmat(true{b},pred{b});
        prob_mat = raw_mat./ sum(raw_mat,2);
        mat_shape = size(prob_mat,1);
        lix_pred = prob_mat(:,1:(mat_shape/2));
        neurometric_curve{g}(b,:) = sum(lix_pred,2);
    end
    
    %     frame=getframe(gcf);
    %     writeVideo(v,frame);
end

%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
%     fn = 'decoding_heat.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])
% close(v)



%%
figure(80);clf
subplot(2,2,1);
imagesc(gofmetrics{g}.cmat)
caxis([0 .4])
set(gca,'xtick',[],'ytick',[]);
colormap turbo
colorbar
axis square

numNeurons_toPlot = [1 5 7 15 20 25]; 
[~,idx] = intersect(numNeurons,numNeurons_toPlot);

my_Map = flipud(jet(length(idx)));
subplot(2,2,2)
for d = 1:length(idx)
    hold on;
    %     shadedErrorBar(gofmetrics{d}.resolution(:,1),gofmetrics{d}.resolution(:,2),gofmetrics{d}.resolution(:,3),'linecolor','b');
    plot(gofmetrics{idx(d)}.resolution(:,1),gofmetrics{idx(d)}.resolution(:,2),'color',my_Map(d,:));
end
set(gca,'ylim',[0 1],'xlim',[0 8],'xtick',0:2:10,'xticklabel',0:.5:5,'ytick',0:.25:1) %hard coded xticklabels for single touch prediction of pole position using population of OL tuned cells
xlabel('mm w/in prediction');ylabel('p (prediction)')
axis square
title('resolution')
legend([num2str(numNeurons(idx)')])

%plot neurometric curve
neuro_mean = cellfun(@nanmean ,neurometric_curve,'uniformoutput',0);
neuro_sem = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),neurometric_curve,'uniformoutput',0);
neuro_std = cellfun(@(x) nanstd(x),neurometric_curve,'uniformoutput',0);
MAE = cellfun(@(x) abs(x-real_psycho_mean),neurometric_curve,'uniformoutput',0);
mean_MAE = cellfun(@(x) mean(mean(x,2)),MAE);
sem_MAE = cellfun(@(x) std(mean(x,2))./sqrt(numIterations),MAE);
std_MAE = cellfun(@(x) std(mean(x,2)),MAE);
[minMAE,minidx] = min(mean_MAE);

MAE_for_anova = cellfun(@(x) median(x,2),MAE,'uniformoutput',0);
[p,~,stats] = anova1(cell2mat(MAE_for_anova),[],'off');
alpha_value = 0.01;
if p<alpha_value
    compare_table = multcompare(stats,'Display','off');
    max_compares = compare_table(any(compare_table(:,[1 2]) == minidx,2),:);
    sig_max_compares = max_compares(max_compares(:,end) < alpha_value,:);
    compare_idx = sig_max_compares(:,[1 2]);
    other_idx = compare_idx(compare_idx ~= minidx);
    
    left_tuning_idx = other_idx(find(other_idx<minidx,1,'last'));
    right_tuning_idx = other_idx(find(other_idx>minidx,1,'first'));
end

subplot(2,2,3)
% %average neurometric distance from average pyschometric
% MAE = cellfun(@(x) mean(abs((x-real_psycho_mean))),neuro_mean);
% [minMAE,minidx] = min(MAE);
% scatter(numNeurons,MAE,'ko')


raw_accuracy = cellfun(@(x) mean(x.meta.trialCorrect),BV);
mean_accuracy = mean(raw_accuracy);
sem_accuracy = std(raw_accuracy)./sqrt(numel(raw_accuracy));
std_accuracy = std(raw_accuracy); 
lix_neuro = cellfun(@(x) x(:,1:20)>.5,neurometric_curve,'uniformoutput',0);
nolix_neuro = cellfun(@(x) x(:,21:40)<=.5,neurometric_curve,'uniformoutput',0);
pcorrect_neuro = cellfun(@(x,y) (sum(x,2) + sum(y,2)) ./ (size(x,2) + size(y,2)),lix_neuro,nolix_neuro,'uniformoutput',0);
[~,p] = cellfun(@(x) ttest2(x,raw_accuracy),pcorrect_neuro);
numNeurons(strfind(p<0.05,[0 1]) + 1)

% shadedErrorBar(numNeurons_toPlot,cellfun(@mean, pcorrect_neuro),cellfun(@(x) std(x)./sqrt(numel(x)),pcorrect_neuro))
errorbar(numNeurons,cellfun(@mean, pcorrect_neuro),cellfun(@(x) std(x),pcorrect_neuro),'ko-')
hold on; errorbar(30,mean_accuracy,std_accuracy,'ko')
set(gca,'xlim',[1 31],'ylim',[0.5 1],'ytick',0:.25:1)
axis square
% shadedErrorBar(numNeurons,mean_MAE,sem_MAE,'k-')
% hold on;scatter(numNeurons(minidx),mean_MAE(minidx),'filled','r')
% try
%     hold on;scatter(numNeurons(left_tuning_idx),mean_MAE(left_tuning_idx),'filled','g')
%     hold on;scatter(numNeurons(right_tuning_idx),mean_MAE(right_tuning_idx),'filled','g')
% end

% set(gca,'ylim',[0 .4],'xtick',0:5:25,'xlim',[0 25])
% xlabel('number of neurons')
% ylabel('mean absolute error')
% axis square

subplot(2,2,4)
for d = 1:length(idx)
    %     subplot(rc(1),rc(2),d)
    hold on;
    h=shadedErrorBar(linspace(-1,1,numel(neuro_mean{idx(d)})),neuro_mean{idx(d)},neuro_sem{idx(d)});
    h.patch.FaceColor = my_Map(d,:);
    %     h.patch.FaceAlpha = .5;
    h.mainLine.Color = [0 0 0];
    %     plot(linspace(-1,1,numel(neuro_mean{d})),neuro_mean{d},'color',boneMap(d,:));
    set(gca,'xlim',[-1 1],'xtick',-1:1:1,'ylim',[0 1],'ytick',0:.25:1)
end

hold on; shadedErrorBar(linspace(-1,1,numel(neuro_mean{d})),real_psycho_mean,real_psycho_sem,'k')
title(['neurometric curve (opt = ' num2str(numNeurons(strfind(p<0.05,[0 1]) + 1)) ')'])
ylabel('lick probability')
xlabel('normalized pole location')
axis square
suptitle('number of neurons affecting prediction of:')

    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'resolution_neurometric_replace.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])




