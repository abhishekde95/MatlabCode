function model = fitsigmoidMLE(X,data)
% Writing a function for fitting a saturating exponential: sigmoid function
% Fitting it using a MLE approach, assuming a binomial distribution
% data is [K N-K]

options.MaxIter = 1000000;
options.MaxFunEvals = 1000000;
options.TolFun = 1e-6;
options.TolX = 1e-6;
guesses = [-3; 3; 1.0];

% [model, fval, success] = fmincon(@mlerr, guesses,[],[],[],[],[-10 -5 eps],[-eps 5 1],[],options);
[model, fval, success] = fminsearch(@mlerr, guesses,options);

% likelihood error function
    function lik = mlerr(input)
        fit = input(3)./(1+exp(input(1)*(X-input(2)))); % Objective Function
        fit(fit>1-eps) = 1-eps;  % Hack to avoid -Inf errors.
        logP = log10([fit(:) 1-fit(:)]);
        lik = sum(sum(data.*logP));
        lik = -lik;
%         keyboard;
    end
end

