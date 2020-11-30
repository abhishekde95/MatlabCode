function [mu1,mu2,sigma2] = calcsignalnoisedist(Hits,Miss,CR,FA,decisioncriterion)
global prophits propFA
error_val = 10000;
options.MaxIter = 1000000;
options.MaxFunEvals = 1000000;
options.TolFun = 1e-4;
options.TolX = 1e-4;
prophits = Hits/(Hits+Miss); 
propFA = FA/(FA+CR); 

[model, fval, success,~,~,~,hess] = fmincon(@finddist, [1 -1 1],[0 0 -1],[eps],[],[],[],[],[],options);
mu1 = model(1);
mu2 = model(2);
sigma2 = model(3);

    function out = finddist(input)
        predictedhits = 1-normcdf(decisioncriterion,input(1),1);
        predictedFA = 1-normcdf(decisioncriterion,input(2),input(3));
        out = (prophits- predictedhits)^2 + (propFA- predictedFA)^2;
    end

end

