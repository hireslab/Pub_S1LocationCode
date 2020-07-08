%adding parameters on top of uber_array; 
function array = struct_editor(array)

for k = 1:length(array)
    curr_array = array{k};
    %find touch indices
    touchOnIdx = [find(curr_array.S_ctk(9,:,:)==1); find(curr_array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(curr_array.S_ctk(10,:,:)==1); find(curr_array.S_ctk(13,:,:)==1)];
    
    dTheta_mat = zeros(curr_array.t,curr_array.k);
    dKappa_mat = zeros(curr_array.t,curr_array.k);
    for i = 1:length(touchOnIdx)
        dTheta_mat(touchOnIdx(i):touchOffIdx(i)) = curr_array.S_ctk(1,touchOnIdx(i):touchOffIdx(i)) - curr_array.S_ctk(1,touchOnIdx(i));
        dKappa_mat(touchOnIdx(i):touchOffIdx(i)) = curr_array.S_ctk(17,touchOnIdx(i):touchOffIdx(i)) - curr_array.S_ctk(17,touchOnIdx(i));
    end
    
   
    array{k}.varNames{18} = 'deltaTheta';
    array{k}.varNames{19} = 'deltaKappa_corrected';
    array{k}.S_ctk(18,:,:) = dTheta_mat;
    array{k}.S_ctk(19,:,:)  = dKappa_mat;
end

%gender addition
for k = 1:length(array)
    if any(strcmp(array{k}.meta.mouseName,{'AH0286','AH0641','AH0669','AH0716','AH0761','AH0762','AH0976','AH1014','AH1015','AH1016','AH1017'}))
        array{k}.meta.gender = 'f';
    elseif any(strcmp(array{k}.meta.mouseName,{'AH0638','AH0657','AH0927','AH0991','AH0992'}))
        array{k}.meta.gender = 'm';
    else
        array{k}.meta.gender = 'unknown';
    end
end
