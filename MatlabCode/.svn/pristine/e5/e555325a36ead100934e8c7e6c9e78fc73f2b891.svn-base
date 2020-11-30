function [final_model,fvalR,fvalOLS] = tmp_linefit(RHO, THETA,not_oog_idx,outofgamut,initial_guess)
    % Author - Abhishek De 3/2017 for fitting a curve in RHO - THETA space
    % such that the fitted curve becomes a line in cartesian plane. As of
    % now just exclude the out of gamut points
    
    options.MaxIter = 1000000000;
    options.MaxFunEvals = 1000000000;
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    [model,fvalOLS] = fminsearch(@(x)lse(x,1), initial_guess, options); % first fitting using LSE
    [final_model,fvalR] = fminsearch(@(x)lse(x,2), model, options); % Then fitting using Robust regression
    % least square error
    function err = lse(input,mode)
        pred = 1./(input * [cos(THETA) sin(THETA)]');
        err = tmp_calcerror(pred,RHO,not_oog_idx,outofgamut,mode);
    end
end
