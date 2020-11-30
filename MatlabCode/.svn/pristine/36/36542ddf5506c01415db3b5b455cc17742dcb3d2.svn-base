function err = calcSSE(pred,RHO,mode)
% This function is same as calcerror.m but does not take into account the
% out of gamut points.
if ~isreal(pred)
   if mode == 1 % LSE
        err = 100000000;
    elseif mode == 3 % Huber's regression 
        K = 1.3456;
        err = numel(pred)*(K^2);
    elseif mode == 2 % Tukey-Bisquare regression
        K = 4.685;
        err = numel(pred)*(K^6);
    end
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
    elseif mode == 3 % Huber regression
        
        tmp = log(pred)-log(RHO');
        K = 1.3456;
        %         keyboard;
        if any(L)
            resid = tmp(~L);
            resid1 = resid(abs(resid)>=K);
            resid2 = resid(abs(resid)<K);
            if isempty(resid1)
                resid1 = 0;
                error1 = 0;
            else
                error1 = sum((K./abs(resid1)).*(K*abs(resid1)-0.5*K^2));
            end
            error2 = sum((resid2.^2)/2);
            error3 = sum(L)*(K^2);% assuming infinite error
            
            err = error1 + error2 + error3;
        else
            resid = tmp;
            resid1 = resid(abs(resid)>=K);
            resid2 = resid(abs(resid)<K);
            if isempty(resid1)
                resid1 = 0;
                error1 = 0;
            else
                error1 = sum((K./abs(resid1)).*(K*abs(resid1)-0.5*K^2));
            end
            error2 = sum((resid2.^2)/2);
            err = error1 + error2;
        end
    elseif mode == 2 % Tukey bisquare regression
        K = 4.685;
        if any(L)
            tmp = log(pred)-log(RHO');
            resid = tmp(~L);
            resid1 = resid(abs(resid)>=K);
            resid2 = resid(abs(resid)<K);
            if isempty(resid1)
                resid1 = 0;
                error1 = 0;
            else
                error1 = numel(resid1)*(K^2)/6;
            end
            weight1 = 1;
            weight2 = 1;
            weight3 = 1;
            error2 = (((K^2)/6)*(1-(1-(resid2.^2/K^2)).^3));
            error3 = sum(L)*(K^2)/6;
            err = weight1*error1 + sum(weight2.*error2) + weight3*error3;
        else
            resid = log(pred)-log(RHO');
            resid1 = resid(abs(resid)>=K);
            resid2 = resid(abs(resid)<K);
            if isempty(resid1)
                resid1 = 0;
                error1 = 0;
            else
                error1 = numel(resid1)*(K^2)/6;
            end
            weight1 = 1;
            weight2 = 1;
            error2 = (((K^2)/6)*(1-(1-(resid2.^2/K^2)).^3));
            err = weight1*error1 + sum(weight2.*error2);
        end
    end
end
end

