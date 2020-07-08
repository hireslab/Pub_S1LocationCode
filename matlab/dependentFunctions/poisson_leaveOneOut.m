function glmModel = poisson_leaveOneOut(glmModel,glmnetOpt,selectedUnits)

for i=1:length(selectedUnits)
    sel_cell = selectedUnits(i); 
    
    coefficients = glmModel{sel_cell}.coeffs.raw;
    numVariables = 1:size(coefficients,1);
    
    fitDevExplained = nan(size(coefficients,1),glmnetOpt.numIterations); 
   
    for g = 2:size(coefficients,1) %for each coefficient 
        sel_coefficients = coefficients(setdiff(numVariables,g),:); %for leave one out 
%         sel_coefficients = coefficients; %for shuffle
        
        for d = 1:glmnetOpt.numIterations
            coeff_check = sel_coefficients(:,d);
            DmatX = glmModel{sel_cell}.predicted.inputX{d};
            full_DmatX = [ones(size(DmatX,1),1) DmatX];
            DmatX_check = full_DmatX(:,setdiff(numVariables,g)); %for leave one out 
%             DmatX_check = full_DmatX;%for shuffle
%             DmatX_check(:,g) = DmatX_check(randperm(size(DmatX,1)),g); %for shuffle
            
            DmatY = glmModel{sel_cell}.predicted.spikeTestRaw{d}; 
            model = exp(DmatX_check*coeff_check);
            mu = mean(DmatY); % null poisson parameter
            nullLogLikelihood = sum(log(poisspdf(DmatY,mu)));
            saturatedLogLikelihood = sum(log(poisspdf(DmatY,DmatY)));
            fullLogLikelihood = sum(log(poisspdf(DmatY,model)));
            fitDevExplained(g,d) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            %devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
        end
    end
    
    glmModel{sel_cell}.gof.LOO_names = {'row = feature','column = iteration #'};
    glmModel{sel_cell}.gof.LOO_de = fitDevExplained;
    glmModel{sel_cell}.gof.LOO_importance = 1 - (glmModel{sel_cell}.gof.LOO_de ./ glmModel{sel_cell}.gof.devExplained); %proportion de explained. 
   
%     fitDevExplained = nan(1,size(coefficients,1)); 
%     meanCoeff = mean(coefficients,2);
%     for g = 2:size(coefficients,1) %for each coefficient
%         sel_coefficients = meanCoeff(setdiff(numVariables,g),:);
%         
%         DmatX = glmModel{sel_cell}.io.DmatXNormalized;
%         full_DmatX = [ones(size(DmatX,1),1) DmatX];
%         DmatX_check = full_DmatX(:,setdiff(numVariables,g));
%         DmatY = glmModel{sel_cell}.io.DmatY;
%         
%         model = exp(DmatX_check*sel_coefficients);
%         mu = mean(DmatY); % null poisson parameter
%         nullLogLikelihood = sum(log(poisspdf(DmatY,mu)));
%         saturatedLogLikelihood = sum(log(poisspdf(DmatY,DmatY)));
%         fullLogLikelihood = sum(log(poisspdf(DmatY,model)));
%         fitDevExplained(g) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
%         %devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
%         
%     end
%     
%     mdl_full = exp(full_DmatX * meanCoeff);
%     mu = mean(DmatY); % null poisson parameter
%     nullLogLikelihood = sum(log(poisspdf(DmatY,mu)));
%     saturatedLogLikelihood = sum(log(poisspdf(DmatY,DmatY)));
%     fullLogLikelihood = sum(log(poisspdf(DmatY,mdl_full)));
%     fitDevExplained_full = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
%     
%         glmModel{sel_cell}.gof.LOO_names = {'row = feature','column = iteration #'};
%     glmModel{sel_cell}.gof.LOO_de = (1-(fitDevExplained ./ fitDevExplained_full)) ./ (nansum(1-(fitDevExplained ./ fitDevExplained_full)));

end

