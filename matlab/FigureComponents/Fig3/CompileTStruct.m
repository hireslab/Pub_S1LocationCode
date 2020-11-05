function touch_struct = CompileTStruct(data_directory, feature_list)

for b = 1:numel(feature_list)
    fileName = ['touch_' feature_list{b}];
    
    if exist([data_directory 'Touch\' fileName '.mat'], 'file')
        % Load existing touch structure
        disp(['Loading touch stucture for ' feature_list{b}])
        load([data_directory 'Touch\' fileName '.mat']);
        touch_struct.(feature_list{b}) = tStruct;
    
    else
        % Build touch structure from scratch
        disp('No touch structure found. Building from scratch...')
        
        if ~exist('U')
            disp('Loading all neural recordings')
            load([data_directory 'Raw\excitatory_clean.mat']); %L5b excitatory cells recorded by Jon and Phil
        end
        
        disp(['Building touch structure for ' feature_list{b}])
        disp('This may take some time so consider loading pre-built structures')
%         selected_cells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U)); % get touch cells
        selected_cells = 1:length(U);
        tStruct = object_location_quantification(U,selected_cells,feature_list{b},'off');      
        save(['touch_' feature_list{b}],'tStruct') 
        touch_struct.(feature_list{b}) = tStruct;
    end
end