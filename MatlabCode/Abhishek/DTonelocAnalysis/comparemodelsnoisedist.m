function LLratio = comparemodelsnoisedist(CRL,FAL,CRNL,FANL,decisioncriterion)
global propFAL propFANL sigmaNL
error_val = 10000;
options.MaxIter = 1000000;
options.MaxFunEvals = 1000000;
options.TolFun = 1e-4;
options.TolX = 1e-4; 
propFAL = FAL/(FAL+CRL);
propFANL = FANL/(FANL+CRNL);

[modelL1, fval1, success1,~,~,~,hess1] = fmincon(@findnoisedist1, [-1 -1 1],[0 0 -1],eps,[],[],[],[],[],options); % allowed to change mu and sigma
[modelL2, fval2, success2,~,~,~,hess2] = fmincon(@findnoisedist2, [-1 -1 1 1],[0 0 -1 0;0 0 0 -1],[eps eps],[],[],[],[],[],options); % only sigma changes
% keyboard;
LLratio = (fval2/fval1); % Degree of freedom
out = chi2inv(0.95,1)<LLratio;

    function out = findnoisedist1(input)
        predictedFANL = 1-normcdf(decisioncriterion,input(1),input(3));
        predictedFAL = 1-normcdf(decisioncriterion,input(2),input(3));
        out = (propFANL- predictedFANL)^2 + (propFAL- predictedFAL)^2;
    end

    function out = findnoisedist2(input)
        predictedFANL = 1-normcdf(decisioncriterion,input(1),input(3));
        predictedFAL = 1-normcdf(decisioncriterion,input(2),input(4));
        out = (propFANL- predictedFANL)^2 + (propFAL- predictedFAL)^2;
    end


end

