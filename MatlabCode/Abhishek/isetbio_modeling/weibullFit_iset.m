function [final_model, final_fval, final_success] = weibullFit_iset(contrasts, data)
% Author - Abhishek De, 2016, CSHL project
error_val = 1000;
options.MaxIter = 1000000;
options.MaxFunEvals = 1000000;
options.TolFun = 1e-50;
options.TolX = 1e-50;
final_model = [];
final_fval = [];
final_success = [];

alphaGuess = .01;
betaGuess = 2; %the apparent slope if fit with line
guesses = [alphaGuess, betaGuess];

lb = [0.00001 0.5];
ub = [1.0 5];

%[model, fval, success] = fminsearch(@sse, guesses, options);
[model, fval, success] = fmincon(@(x)sser(x,contrasts,data), guesses,[-1 0 ;0 -1],[0 0],[],[],lb,ub,[],options);

if fval < error_val
    error_val  = fval;
    final_fval = fval;
    final_model = model;
    final_success = success;
end
end

% least square error function
function sse = sser(input,contrasts,data)

a = input(1);
b = input(2);
fit = 1 - 0.5.*exp(-((contrasts./a).^b));
sse = sum((data-fit).^2);

end