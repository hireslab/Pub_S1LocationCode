function [glmModel] = designMatrixBlocks_simplified(selectedArray,glmnetOpt,glmModel)

preDecisionTouches = preDecisionTouchMat(selectedArray);

for i = 1:length(selectedArray)
    array = selectedArray{i};
    [tVar] = atTouch_sorter(array,glmnetOpt.buildIndices,preDecisionTouches{i});
    it = tVar.allTouches.S_ctk(:,1:end-1);
    dt = tVar.allTouches.dtS_ctk;
    
    response_window = find(glmnetOpt.buildIndices == array.meta.touchProperties.responseWindow(1)) : find(glmnetOpt.buildIndices == array.meta.touchProperties.responseWindow(2)); 
    response =  sum(tVar.allTouches.R_ntk(:,response_window),2);
    
    %Feature blocks for design matrix construction 
    glmModel{i}.io.components.angle = it(:,1);
    glmModel{i}.io.components.phase = it(:,5);
    glmModel{i}.io.components.amplitude = it(:,3);
    glmModel{i}.io.components.midpoint = it(:,4);
    glmModel{i}.io.components.kappa = it(:,6);
    glmModel{i}.io.components.pt_velocity = it(:,2);
    glmModel{i}.io.components.DKappa = dt(:,2);
    glmModel{i}.io.components.DTheta = dt(:,3); 
    
    %Responses
    glmModel{i}.io.DmatY = response;
    
    %raw 
    glmModel{i}.raw.pole = tVar.allTouches.S_ctk(:,end); 
end
