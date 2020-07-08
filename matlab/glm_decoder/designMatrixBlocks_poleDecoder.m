function mdl = designMatrixBlocks_poleDecoder(mdl,tuning_structure,glmnetOpt)

tuned_cells = find(cellfun(@(x) x.is_tuned,tuning_structure)==1);
sampledSpace = cellfun(@(x) range(x.stim_response.values(:,1)),tuning_structure(tuned_cells)) ./ 2;
units_2_use = tuned_cells(sampledSpace > glmnetOpt.pctSamplingThreshold);

for g = 1:length(units_2_use)
    current = tuning_structure{units_2_use(g)};
    stim = normalize_var(current.stim_response.values(:,1),-1,1);
    stretchResolution = linspace(-1,1,glmnetOpt.interpResolution)'; 
    interpResp = interp1(stim,current.stim_response.values(:,2),stretchResolution);
    interpStd = interp1(stim,current.stim_response.values(:,3),stretchResolution);
            
    
    resampX = [];
    resampY = [];
    if strcmp(glmnetOpt.samplingOption,'poisson')
        for d = 1:glmnetOpt.numberResamples
            resampX = [resampX ; poissrnd(interpResp)];
            resampY = [resampY ; (1:length(interpResp))'];
        end
    elseif strcmp(glmnetOpt.samplingOption,'normal')
        for d = 1:glmnetOpt.numberResamples
            resampX = [resampX ; normrnd(interpResp,interpStd)];
            resampY = [resampY ; (1:length(interpResp))'];
        end
    else
        error('input glmnetOpt.samplingOption as "normal" or "poisson"')
    end
    
    
    mdl{g}.buildBlocks.valueNames = current.stim_response.varNames;
    mdl{g}.buildBlocks.values.raw = current.stim_response.values;
    mdl{g}.buildBlocks.values.interpolated = [stretchResolution interpResp interpStd]; 
    
    mdl{g}.io.DmatX = resampX;
    mdl{g}.io.DmatY = resampY;
    
    mdl{g}.params = glmnetOpt;
    mdl{g}.params.cellNum = units_2_use(g); 
end