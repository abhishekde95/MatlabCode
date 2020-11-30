function [final_model,fval] = linefit_AD(RHO, THETA, outofgamut,initial_guess)
    % Author - Abhishek De 3/2017 for fitting a curve in RHO - THETA space
    % such that the fitted curve becomes a line in cartesian plane. As of
    % now just exclude the out of gamut points
    
    options.MaxIter = 100000000;
    options.MaxFunEvals = 100000000;
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    num_trials = 200;
    fval = 1000;
    final_model = [];
    for ii = 1:num_trials
        [model,val] = fminsearch(@lse, initial_guess, options);
        initial_guess = randn(3,1);
        if val < fval
            final_model = model;
            fval = val;
        end
    end
%     final_model = model;
%     fval = val;
    
    % least square error
    function err = lse(input)
        a = input(1); 
        b = input(2);
        c = input(3);
        pred = 1*(b*sin(THETA*pi/180)+a*cos(THETA*pi/180))/c;
        resid = log(abs(pred))-log(1./RHO);
        resid(outofgamut & (abs(pred)<1./RHO))=0;
        err = sum(resid.^2);
    end
end