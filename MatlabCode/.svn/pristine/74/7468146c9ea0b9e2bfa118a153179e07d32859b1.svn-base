function [final_model, final_fval, final_success, hess] = weibullFit_AD(contrasts, data, I, fit_mode, method, guesses)
    % Originally written by Charlie, modified by Abhishek De, 2016
    % Performs MLE assuming a mean using Weibull cdf and variance around
    % each point to be Binomial distributed
    
    % Has the option of using two different kind of search methods -
    % fminsearch and Simulated Annealing
    error_val = 10000;
    options.MaxIter = 1000000;
    options.MaxFunEvals = 1000000;
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    final_model = [];
    final_fval = [];
    final_success = [];
    hess = 1000;
    
    %determine guesses for alpha and beta
    %     keyboard
    if nargin<6
        if strcmp(fit_mode,'multiple')
            alphaGuess = .2;
            betaGuess = 1; %the apparent slope if fit with line
            alpha2Guess = .2;
            guesses = [alphaGuess, betaGuess, alpha2Guess]; 
            
        elseif strcmp(fit_mode,'single')
            alphaGuess = .1;
            betaGuess = 1; %the apparent slope if fit with line
            guesses = [alphaGuess, betaGuess]; 
            
        end
    end
    
    if strcmp(fit_mode,'multiple')
        lb = [0 0 0];
        ub = [1 10 1];
        iter = 10;
    elseif strcmp(fit_mode,'single')
        lb = [0 0];
        ub = [1 10];
        iter = 2;
    end
    
    %clean up the inputs. We won't consider data that comes from zero
    %contrast trials, so chuck these data.
    if numel(data) == length(data)
        data = data(:); %make SSE data a column so that the following deletions work for SSE & MLE
    end
    l_zeros = softEq(0, contrasts, 5); %anything less than 0.00001 gets considered a zero
    contrasts(l_zeros) = [];
    data(l_zeros,:) = [];
    %call the appropriate fitting routine and perform the fit
    for i = 1:iter
        if i > 1 
            if strcmp(fit_mode,'multiple')
                alphaGuess = rand();
                betaGuess = randi(2);
                alpha2Guess = rand();
                guesses = [alphaGuess, betaGuess, alpha2Guess];
            elseif strcmp(fit_mode,'single')
                alphaGuess = rand();
                betaGuess = randi(2);
                guesses = [alphaGuess, betaGuess];
            end
        end
                
        if strcmp(method,'fminsearch')
            [model, fval, success] = fminsearch(@mlerr, guesses, options);
            try
                [hess,~] = hessian(@mlerr,model);
            catch
                keyboard;
            end

        elseif strcmp(method, 'sim_anneal')
            % Doesn't work well at all. I would rather use fminsearch
            [model, fval, success] = simulannealbnd(@mlerr, guesses);
        elseif strcmp(method,'fmincon')
%             keyboard
            [model, fval, success,~,~,~,hess] = fmincon(@mlerr, guesses,[-1 0 0; 0 -1 0; 0 0 -1],[eps eps eps],[],[],lb,ub,[],options);
 
        end

        if strcmp(fit_mode,'single') && fval < error_val
            error_val  = fval;
            final_fval = fval;
            final_model = model;
            final_success = success;
        elseif strcmp(fit_mode,'multiple') && fval < error_val
            if model(3)>0
                error_val  = fval;
                final_fval = fval;
                final_model = model;
                final_success = success;
            else
                i = i-1;
            end
        end
    end
       
    % likelihood error function
    function lik = mlerr(input)
        if strcmp(fit_mode,'single')
            a = input(1);
            b = input(2);
            fit = 1 - 0.5.*exp(-((contrasts./a).^b));
        elseif strcmp(fit_mode,'multiple')
            a = input(1);
            b = input(2);
            c = input(3);
            idx  = (I == 1);
            fit = [1 - 0.5*exp(-((contrasts(idx)/a).^b)); 1 - 0.5*exp(-((contrasts(~idx)/c).^b))];    
        end
        if any(fit > 1)
            lik = inf;
            return
        end
        if a < 0
            lik = inf;
            return
        end
        fit(fit>1-eps) = 1-eps;  % Ugly hack to avoid log10(0) = -Inf errors.
        logP = log10([fit(:), 1-fit(:)]);
        lik = sum(sum(data.*logP));
        lik = -lik;
%         if lik < 0
%             keyboard;
%         end
    end
end