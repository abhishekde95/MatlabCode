function err = tmp_calcerror_wo_oog(pred,RHO,mode)
% This function is same as calcerror.m but does not take into account the
% out of gamut points.
if all(imag(pred)~=0)
    err = 100000000;
else
    %     keyboard;
    L = pred<0;    
    if mode == 1 % LSE
        if any(L)
            tmp = log(pred)-log(RHO');
            resid = tmp(~L);
            err = sum(resid.^2) + 1000*sum(L);
        else
            resid = log(pred)-log(RHO');
            err = sum(resid.^2);
        end
    elseif mode == 2 % Tukey bisquare regression
        
        if any(L)
            tmp = log(pred)-log(RHO');
            resid = tmp(~L);
            resid1 = resid(abs(resid)>=4.685);
            resid2 = resid(abs(resid)<4.685);
            if isempty(resid1)
                resid1 = 0;
                error1 = 0;
            else
                error1 = numel(resid1)*(4.685^2)/6;
            end
            error1 = numel(resid1)*(4.685^2)/6;
            error2 = sum(((4.685^2)/6)*(1-(1-(resid2.^2/4.685^2)).^3));
            error3 = sum(L)*(4.685^2)/6;
            err = error1 + error2 + error3;
        else
            resid = log(pred)-log(RHO');
            resid1 = resid(abs(resid)>=4.685);
            resid2 = resid(abs(resid)<4.685);
            if isempty(resid1)
                resid1 = 0;
                error1 = 0;
            else
                error1 = numel(resid1)*(4.685^2)/6;
            end
            error2 = sum(((4.685^2)/6)*(1-(1-(resid2.^2/4.685^2)).^3));
            err = error1 + error2;
        end
    end
end
end

