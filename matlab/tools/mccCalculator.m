function mcc = mccCalculator(true,predicted)

cmat=confusionmat(true,predicted);
TP = cmat(1);FP = cmat(3);
TN = cmat(4);FN = cmat(2);

%MCC FULL
top  = TP*TN - FP*FN;
bottom = sqrt((TP+FP) * (TP+FN) * (TN+FP)* (TN+FN));
mcc = top./bottom;
if isnan(mcc)
    mcc = 0 ;
end
