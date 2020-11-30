function [alpha, beta, gamma, success, modErrs,fval] = weibullFitforDToneloc(contrasts, data, errFun, guesses,I)
    
    options.MaxIter = 10000000;
    options.MaxFunEvals = 10000000;
    
    %determine guesses for alpha and beta
    useModGamma = 1;

    alphaGuess = .1;
    betaGuess = 1; %the apparent slope if fit with line
    gammaGuess = 0.5;

    if ~exist('guesses')
        if ~exist('I')
            guesses = [alphaGuess, betaGuess, gammaGuess];
        else
            guesses = [alphaGuess, alphaGuess+0.1, betaGuess/2, gammaGuess];
        end
    end
    if nargin == 4 & ~isempty(guesses)
        if length(guesses) < 3
            useModGamma = 0;
        end
    end
      
    %call the appropriate fitting routine and perform the fit
    if strcmpi(errFun, 'mle')
        if ~exist('I')
            A = [1 0 0;...
                -1 0 0;...
                0 1 0;...
                0 -1 0;...
                0 0 1;...
                0 0 -1];
            b = [10; -eps; 10; -eps; 1; eps-1];
        else
            A = [1 0 0 0;...
                -1 0 0 0;...
                0 1 0 0;...
                0 -1 0 0;...
                0 0 1 0;...
                0 0 -1 0;...
                0 0 0 1;...
                0 0 0 -1];
            b = [10; -eps; 10; -eps; 10; -eps; 1; eps-1];
        end
       [model, fval, success] = fmincon(@mlerr,guesses,A,b,[],[],[],[],[],options);
       if exist('I')
           fval_init = fval;
           for ii = 1:5 
               guesses = [rand(1), rand(1), 5*rand(1), rand(1)];
               try
                   [model1, tmp, success1] = fmincon(@mlerr,guesses,A,b,[],[],[],[],[],options);
               catch
                   keyboard;
               end
               if tmp < fval_init
                   model = model1;
                   fval_init = tmp;
                   success = success1;
               end
           end
           fval= fval_init;
       end
           
           
           
       if nargout > 4 %if the SE of the model estimates are desired
          [hess, ~] = hessian(@mlerr, model);
          if rcond(hess) > (2*eps) %i.e., the matrix is poorly invertible
              invHess = inv(hess);
              modErrs = sqrt(abs(diag(invHess)));
          else
              modErrs = nan(2,1);
          end
       end
    end
    
    
    %determine the outputs of the function
    if ~exist('I')
        alpha = model(1);
        beta = model(2);
        gamma = model(3);
    else
        alpha = [model(1) model(2)];
        beta = model(3);
        gamma = model(4);
    end
    
    % likelihood error function
    function lik = mlerr(input) 
        if ~exist('I')
            a = input(1);
            b1 = input(2);
            g = input(3);
            fit = g*(1-exp(-((contrasts./a).^b1)));
            if a < 0
                lik = inf;
                return
            end
        else
            a1 = input(1);
            a2 = input(2);
            b0 = input(3);
            g = input(4);
            fit = g*([1 - exp(-((contrasts(I)./a1).^b0)); 1 - exp(-((contrasts(~I)./a2).^b0))]);
            if a1 < 0 | a2 < 0
                lik = inf;
                return
            end
        end
        if any(fit > 1)
           lik = inf; 
           return
        end
       fit(fit>1-eps) = 1-eps;  % Ugly hack to avoid log10(0) = -Inf errors.
       logP = log10([fit(:), 1-fit(:)]);
       lik = sum(sum(data.*logP));
       lik = -lik;
    end
end