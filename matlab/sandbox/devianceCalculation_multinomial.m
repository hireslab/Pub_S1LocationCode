p = [1/2 1/3 1/6];
n = 10;
x1 = 0:n;
x2 = 0:n;
[X1,X2] = meshgrid(x1,x2);
X3 = n-(X1+X2);

Y = mnpdf([X1(:),X2(:),X3(:)],repmat(p,(n+1)^2,1))

mnpdf(trainDmatY',ones(1,length(trainDmatY))*.1)

%Built to be inserted into multinomial_classifier function 
% X and PROB are m-by-k matrices or 1-by-k vectors

% where k is the number of multinomial bins or categories (e.g. 10 here)

% Each row of PROB must sum to one 
% and the sample sizes for each observation (rows of X) are given by the row sums sum(X,2)

trainpredValues = sigmoid([ones(length(trainDmatX),1) trainDmatX] * mdl.fitCoeffs{h})
sumToOne = trainpredValues./sum(trainpredValues,2);
largeY = repmat(1:10,size(trainpredValues,1),1);
trainnullLL = sum(log(mnpdf(1:10,ones(size(sumToOne)).*.1)))
trainfullLL = sum(log(mnpdf(1:10,sumToOne)))

 
actualNull = cv.glmnet_fit.nulldev
actualFull = cv.glmnet_fit.dev(iLambda)



trainnullLL = sum(log(mnpdf(trainDmatY,1,(ones(length(trainDmatY),1)*mean(trainDmatY)))));
trainpredProb = sigmoid([ones(length(trainDmatX),1) trainDmatX] *fitCoeffs);
trainfullLL = sum(log(binopdf(trainDmatY,1, trainpredProb)));
trainDevianceExplained = 1-(-trainfullLL)/(-trainnullLL);
devExplained = cv.glmnet_fit.dev(iLambda);
