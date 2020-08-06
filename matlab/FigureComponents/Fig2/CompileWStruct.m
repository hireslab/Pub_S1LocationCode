function whisk_struct = CompileWStruct(data_directory)
hilbert_feature = {'angle','phase','midpoint','amplitude','velocity'};

for b = 1:numel(hilbert_feature)
    fileName = ['whisk_' hilbert_feature{b} '_window'];
    
    if exist([data_directory 'Whisking\' fileName '.mat'], 'file')
        disp(['Loading whisking stucture for ' hilbert_feature{b}])
        load([data_directory 'Whisking\' fileName '.mat']);
        whisk_struct.(hilbert_feature{b}) = wStruct;
    else
        disp('No whisking structure found. Building from scratch...')
        
        if ~exist('U')
            disp('Loading complete all neural recordings')
            load([data_directory 'Raw\excitatory_all.mat']); %L5b excitatory cells recorded by Jon and Phil
        
        end
        disp(['Building whisking structure for ' hilbert_feature{b}])
        disp('This may take some time so consider loading pre-built structures')
        whisk_struct.(hilbert_feature{b}) = whisking_location_quantification(U,1:length(U),hilbert_feature{b},'off');
        
    end
end
