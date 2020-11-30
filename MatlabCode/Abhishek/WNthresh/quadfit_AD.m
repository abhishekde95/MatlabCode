function [final_model,fval] = quadfit_AD(RHO, THETA, outofgamut,initial_guess)
    % Author - Abhishek De 3/2017 for fitting a curve in RHO - THETA space
    % such that the fitted curve becomes a line in cartesian plane. As of
    % now just exclude the out of gamut points
    
    options.MaxIter = 100000000;
    options.MaxFunEvals = 100000000;
    options.TolFun = 1e-6;
    options.TolX = 1e-6;
    num_trials = 1;
    fval = 1000;
    final_model = [];
    for ii = 1:num_trials
        [model,val] = fminsearch(@lse, initial_guess, options);
        initial_guess = initial_guess + rand(1,3);
%         if val < fval
            final_model = model;
            fval = val;
%         end
    end
    
    % least square error
    function err = lse(input)
%         keyboard;
        pred = 1./(input * [cos(THETA).^2 sin(THETA).^2 cos(THETA).*sin(THETA)]');
        L = pred<0;
        pred(L) = 10000000;
        resid = log(pred)-log(RHO');
        resid(outofgamut)=0;
        err = sum(resid.^2)+ 2*sum(L);
    end
end
