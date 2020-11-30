function [final_model, fval] = ellipsefit_AD(RHO, THETA, outofgamut)
    % Author - Abhishek De 3/2017 for fitting a curve in RHO - THETA space
    % such that the fitted curve becomes an "ELLIPSE" in cartesian plane.
    % Have incorporated the outgamut penalty
    options.MaxIter = 100000000;
    options.MaxFunEvals = 100000000;
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    num_trials = 200;
    initial_guess = [1; 1; 1];
    fval = 100;
    final_model =[];
    for ii = 1:num_trials
        [model,val] = fminsearch(@lse, initial_guess, options);
        initial_guess = randn(3,1);
        if val<fval
            fval = val;
            final_model = model;
        end 
    end
    
    % least square error
    function err = lse(input)
        a = input(1); % Major/minor axis length
        b = input(2);  % Major/minor axis length
        c = input(3);
        A = cos(THETA*pi/180);
        B = sin(THETA*pi/180);
        predr2 = 1./(((A+c*B)/a).^2 + ((c*A-B)/b).^2);
        pred = sqrt(predr2);
        resid = log(abs(pred))-log(RHO);
        resid(outofgamut & (abs(pred)>=RHO))=0;
        err = sum(resid.^2);
    end

end

