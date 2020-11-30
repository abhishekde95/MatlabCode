function [model] = sigmoidalfit_AD(x,y)
% Naka-Rushton fit to be precise
options.MaxIter = 1000000;
options.MaxFunEvals = 1000000;
options.TolFun = 1e-4;
options.TolX = 1e-4;
guesses = [1 mean(x) 2];
model = fminsearch(@sigmoidalerr,guesses,options);

    function err = sigmoidalerr(params)
        A = params(1);
        c50 = params(2);
        n = params(3);
        output = A*(x.^n)./(x.^n+c50.^n);
        err = sum((output - y).^2,1);
    end

end

