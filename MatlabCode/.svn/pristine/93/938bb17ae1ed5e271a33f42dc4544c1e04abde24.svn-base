function [final_model,fvalR] = tmp_quadfit(RHO,THETA,not_oog_idx,outofgamut,initial_guess)
    % Author - Abhishek De 1/18 for fitting a curve in RHO - THETA space
    % such that the fitted curve becomes a line in cartesian plane. As of
    % now just exclude the out of gamut points
    
    options.MaxIter = 10000000;
    options.MaxFunEvals = 1000000;
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    final_model = [];
    fval = 10000;
    valOLS = fval;
    count = 30;
    while count>=0
        [modeltmp,valOLStmp] = fminunc(@(x)lse(x,1),initial_guess,options); % first fitting using LSE
        if valOLStmp < valOLS
            valOLS = valOLStmp;
            model = modeltmp;
        end
        initial_guess = randn(1,5);
        count = count - 1;
    end
    [model,valR] = fminsearch(@(x)lse(x,2), model, options); % Then fitting using Robust regression
    if valR < fval && isreal(valR)
        final_model = model;
        fvalR = valR;
        fvalOLS = valOLS;
    end
   
    clear RHO THETA not_oog_idx
    % least square error
    function err = lse(input,mode)
        p = [input(1)*(cos(THETA).^2)+input(2)*(sin(THETA).^2)+input(3)*(cos(THETA).*sin(THETA)) input(4)*cos(THETA)+input(5)*sin(THETA) -1*ones(numel(THETA),1)];
        pred = [];
        for ii = 1:size(p,1)
            rts = roots(p(ii,:));
            if all(imag(rts)==0) & ~isempty(rts)
                if all(rts>0)
                    r = 100000; %min(rts);
                elseif all(rts<0)
                    r = 100000;
                else
                    r = max(rts);
                end
            else
                r = 100000;
            end
            pred = [pred r];
        end

        err = tmp_calcerror(pred,RHO,not_oog_idx,outofgamut,mode);
    end
end

