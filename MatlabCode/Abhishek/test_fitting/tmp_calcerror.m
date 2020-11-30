function err = tmp_calcerror(pred,RHO,not_oog_idx,outofgamut,mode)
% Takes into account the out of gamut points for fitting the data

if all(imag(pred(not_oog_idx))~=0)
    err = 100000000;
else
    L = pred(not_oog_idx)<0 | imag(pred(not_oog_idx))~=0;
    L1 = pred>0 & pred<RHO' & outofgamut' & imag(pred)==0;
    try
        oogresid = (log(pred(L1))-log(RHO(L1)')); % OOG errors if the prediction is less than the out of gamut contrast
    catch
        keyboard;
    end
    if mode == 1 % LSE
        if any(L)
            tmp = log(pred(not_oog_idx))-log(RHO(not_oog_idx)');
            resid = [tmp(~L) oogresid];
            err = sum(resid.^2);% + 1000*sum(L);
        else
            resid = log(pred(not_oog_idx))-log(RHO(not_oog_idx)');
            resid = [resid oogresid];
            err = sum(resid.^2);
        end
    elseif mode == 2 % Tukey bisquare regression
        if any(L)
            tmp = log(pred(not_oog_idx))-log(RHO(not_oog_idx)');
            resid = [tmp(~L) oogresid];
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
            resid = log(pred(not_oog_idx))-log(RHO(not_oog_idx)');
            resid = [resid oogresid];
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

