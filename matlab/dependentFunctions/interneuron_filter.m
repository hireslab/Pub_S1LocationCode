load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Raw\excitatory_all.mat')
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Raw\';
T_directory = [data_directory 'Interneurons\TArrays\'];
T_array_list = dir([T_directory '/*.mat']);

interneuron_list = cell(1,length(T_array_list));

for g = 1:length(T_array_list)
    load([T_directory T_array_list(g).name])
    interneuron_list{g}=[T.cellNum T.cellCode];
end

U_list = cellfun(@(x) [x.meta.cellName x.meta.cellCode],U,'uniformoutput',0);

U_interneuron_filter = ones(1,length(U));
for g = 1:length(U_list)
    for k = 1:length(interneuron_list)
        if U_list{g} == interneuron_list{k}
            U_interneuron_filter(g) = 0;
        end
    end
end

U = U(logical(U_interneuron_filter)); 