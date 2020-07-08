function [trainIdx,testIdx] = stratifiedIndices(DmatY,proportionTrain)

classes = unique(DmatY);

trainIdx = cell(1,numel(classes)); 
testIdx = cell(1,numel(classes)); 

for d = 1:numel(classes)
    classidx = find(DmatY==classes(d));
    shuffIdx = classidx(randperm(length(classidx)));
    
    trainIdx{d} = shuffIdx(1:round(length(classidx)*proportionTrain));
    testIdx{d} = setdiff(shuffIdx,trainIdx{d});
end


