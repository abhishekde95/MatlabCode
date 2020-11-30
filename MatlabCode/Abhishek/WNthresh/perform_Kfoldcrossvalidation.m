function [h1,p1,MSEmodel1,MSEmodel2] = perform_Kfoldcrossvalidation(K,responseS1,rawrgb_onbasisvec1)
class0resp = responseS1(responseS1==0);
class1resp = responseS1(responseS1==1);
c0 = cvpartition(class0resp,'KFold',K);
c1 = cvpartition(class1resp,'KFold',K);
class0stim = rawrgb_onbasisvec1(responseS1==0);
class1stim = rawrgb_onbasisvec1(responseS1==1);
MSEmodel1 = []; MSEmodel2 = [];
for ii = 1:K
    stim = [class0stim(:,c0.training(ii)) class1stim(:,c1.training(ii))];
    response = [class0resp(c0.training(ii)) class1resp(c1.training(ii))];
    teststim = [class0stim(:,c0.test(ii)) class1stim(:,c1.test(ii))];
    testresponse = [class0resp(c0.test(ii)) class1resp(c1.test(ii))];
    model1 =  fitglm(stim',response','linear','Distribution','binomial','Link','logit');
    model2 =  fitglm(stim',response','quadratic','Distribution','binomial','Link','logit');
    modelprediction1 = predict(model1,teststim'); % linear model prediction
    modelprediction2 = predict(model2,teststim'); % quadratic model prediction
    MSEmodel1 = [MSEmodel1; mean((modelprediction1-testresponse').^2)]; % mean squared error
    MSEmodel2 = [MSEmodel2; mean((modelprediction2-testresponse').^2)];
end
[h1,p1] = ttest(MSEmodel1,MSEmodel2);

end

