function [h1,p1,MSE1,MSE2] = performfitting_resampling(K,responseS1,rawrgb_onbasisvec1)
class0resp = responseS1(responseS1==0);
class1resp = responseS1(responseS1==1);
class0stim = rawrgb_onbasisvec1(:,responseS1==0);
class1stim = rawrgb_onbasisvec1(:,responseS1==1);
MSE1 = []; MSE2 = [];
for ii = 1:K
    ind = datasample(1:numel(class0resp),numel(class1resp),'Replace',false);
    stim = [class0stim(:,ind) class1stim];
    response = [class0resp(ind) class1resp];
    model1 =  fitglm(stim',response','linear','Distribution','binomial','Link','logit');
    model2 =  fitglm(stim',response','quadratic','Distribution','binomial','Link','logit');
    probval1 = predict(model1,rawrgb_onbasisvec1');
    probval2 = predict(model2,rawrgb_onbasisvec1');
    MSE1 = [MSE1; sum((probval1-responseS1').^2)];
    MSE2 = [MSE2; sum((probval2-responseS1').^2)];
end
[h1,p1] = ttest(MSE1,MSE2);
end

