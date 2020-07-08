function [glmModel] = designMatrixBlocks_session(selectedArray,glmnetOpt,glmModel)


for i = 1:length(selectedArray)
    currentCell = selectedArray{i};
     
    %Defining features
    raw.touch = ~isnan(squeeze(currentCell.S_ctk(9,:,:))) + ~isnan(squeeze(currentCell.S_ctk(12,:,:)));
    raw.spikes = squeeze(currentCell.R_ntk);
    raw.phase = squeeze(currentCell.S_ctk(5,:,:)); %these are at touch features we care about
    raw.amplitude = squeeze(currentCell.S_ctk(3,:,:));  %these are at touch features we care about
    raw.midpoint = squeeze(currentCell.S_ctk(4,:,:)) ; %these are at touch features we care about
    raw.angle = squeeze(currentCell.S_ctk(1,:,:)) ; %these are at touch features we care about
    raw.kappa = squeeze(currentCell.S_ctk(17,:,:)) ; %kappa at touch
    raw.DKappa = squeeze(currentCell.S_ctk(6,:,:));
    raw.DTheta = squeeze(currentCell.S_ctk(18,:,:));
    
    all_fields = fields(raw);
    
    %toss trials where more than 20% of samples are missing 
    toss_trials = cell(1,length(all_fields));
    for d = 1:length(all_fields)
        toss_trials{d} = find((sum(isnan(raw.(all_fields{d}))) ./ currentCell.t) > .2); 
    end
    toss_trials = unique(cell2mat(toss_trials)); 

    ds = {};
   %Downsampling 
   if mod(numel(raw.touch),glmnetOpt.downsampling_rate) == 0
       for g = 1:length(all_fields)
           raw.(all_fields{g})(:,toss_trials) = []; %Tossing out empty trials 
           if strcmp(all_fields{g},'touchMat') || strcmp(all_fields{g},'spikes')
               ds_tmp = nansum (reshape(raw.(all_fields{g})(:), glmnetOpt.downsampling_rate,...
                   numel(raw.(all_fields{g}))./glmnetOpt.downsampling_rate));
           elseif strcmp(all_fields{g},'DKappa') || strcmp(all_fields{g},'DTheta')
               reshaped_tmp = reshape(raw.(all_fields{g})(:), glmnetOpt.downsampling_rate,...
                   numel(raw.(all_fields{g}))./glmnetOpt.downsampling_rate);
               [~,max_idx] = max(abs(reshaped_tmp));
               ds_tmp = reshaped_tmp(sub2ind(size(reshaped_tmp),max_idx,1:size(reshaped_tmp,2)));
           else
               ds_tmp = nanmean( reshape(raw.(all_fields{g})(:), glmnetOpt.downsampling_rate,...
                   numel(raw.(all_fields{g}))./glmnetOpt.downsampling_rate));
           end
           
           if strcmp(glmnetOpt.interpOption,'yes') && any(isnan(ds_tmp))
               disp(['interpolating for ' num2str(sum(isnan(ds_tmp))) ' timepoint(s) in ' all_fields{g}])
               nanx = isnan(ds_tmp);
               t    = 1:numel(ds_tmp);
               ds_tmp(nanx) = interp1(t(~nanx), ds_tmp(~nanx), t(nanx));
               if ~isempty(find(isnan(ds_tmp)))
                   disp(['interpolating ' num2str(numel(find(isnan(ds_tmp)))) ' edge timepoint w/ n-1'])
                  ds_tmp(find(isnan(ds_tmp))) = ds_tmp(find(isnan(ds_tmp))-1);
               end
               
           end
            
           ds.(all_fields{g}) = reshape(ds_tmp ,currentCell.t ./ glmnetOpt.downsampling_rate, currentCell.k - numel(toss_trials));
       end
   else
       error('downsampling rate must be divisible by trial length')
   end
    

    %lag shifting and then populating output matrix 
    all_fields_ds = fields(ds);
    shift_values = glmnetOpt.shift;
    for b = 1:length(all_fields_ds)
        cs_nan_pad = [];
        nan_pad = [nan(abs(min(shift_values)),currentCell.k - numel(toss_trials)) ; ds.(all_fields_ds{b}) ; nan(abs(min(shift_values)),currentCell.k - numel(toss_trials))];
        for d = 1:length(shift_values)
            cs_nan_pad(:,d) = circshift(nan_pad(:),shift_values(d));
        end
        if strcmp(all_fields_ds{b},'spikes')
            glmModel{i}.io.DmatY = cs_nan_pad(~any(isnan(cs_nan_pad),2), shift_values==0);
        else
            glmModel{i}.io.components.(all_fields_ds{b}) = cs_nan_pad(~any(isnan(cs_nan_pad),2),:);
        end
        
    end
    
    glmModel{i}.io.tossed_trials = toss_trials; 
    glmModel{i}.raw.pole = currentCell.meta.motorPosition;
    glmModel{i}.raw.pole(toss_trials) = []; 
    
end