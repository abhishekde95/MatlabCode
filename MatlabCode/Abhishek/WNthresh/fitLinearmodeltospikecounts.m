function [model,fval,FRsurface,FRpredictedvalues] = fitLinearmodeltospikecounts(WTS,SPC,linemodelparams)
% Fitting a linear model to the spike counts
% options.MaxIter = 1000000000;
% options.MaxFunEvals = 1000000000;
% options.TolFun = 1e-4;
% options.TolX = 1e-4;
options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
initial_guess = [linemodelparams max(SPC) 1 0.5 0];

A1 = [-1 0 0 0 0 0;...
    0 -1 0 0 0 0;...
    0 0 -1 0 0 0;...
    0 0 1 0 0 0;...
    0 0 0 -1 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 -1 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 -1];
b = [-eps -eps -max(SPC)/2 5*max(SPC) -eps 5 -0.003 5 -eps];
[model,fval] = fmincon(@(x)lse(x,WTS,SPC),initial_guess,A1,b,[],[],[],[],[],options);
[X,Y] = meshgrid(linspace(-2,2,51));
FRsurface = (model(3) * (max(model(1)*X+model(2)*Y,0).^model(4))./((max(model(1)*X+model(2)*Y,0).^model(4)) + model(5)^model(4))) + model(6);
FRpredictedvalues = (model(3) * (max(model(1)*WTS(:,1)+model(2)*WTS(:,2),0).^model(4))./((max(model(1)*WTS(:,1)+model(2)*WTS(:,2),0).^model(4)) + model(5)^model(4))) + model(6);
    
% least square error
    function err = lse(input,wts,spikecounts)
        A = input(1);
        B = input(2);
        Rmax = input(3);
        N = input(4);
        Chalf = input(5);
        bl = input(6);
        generator_signal = A*wts(:,1)+B*wts(:,2);
        mu = (Rmax * (max(generator_signal,0).^N)./((max(generator_signal,0).^N) + Chalf^N)) + bl;
        mu(mu==0) = eps;
        logprob = spikecounts.*log(mu) - mu - log(factorial(spikecounts)); % Poisson error model
        err = -1*sum(logprob); % negative log-likelihood
    end
end

