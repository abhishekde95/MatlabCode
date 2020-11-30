% called to estimate the parameters of a line fit to data

function [estimates, model] = fitStraight(x, y)

options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'FunValCheck', 'on');
start_point = rand(1, 2);  
model = @expfun;
estimates = fminsearch(model, start_point, options); 
     function [sse, FittedCurve, p, ErrorVector] = expfun(p)
        FittedCurve = p(1) .*x + p(2); % p(1) = slope % p(2) = y-intercept
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2); 
        p = [p(1), p(2)];
    end
end