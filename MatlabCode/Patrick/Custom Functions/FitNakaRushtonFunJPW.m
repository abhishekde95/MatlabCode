 function [f] = FitNakaRushtonFunJPW(params,contrast,response,FITSTR,error,showplot)
% [f] = FitNakaRushtonFun(params,contrast,response,FITSTR)
% 
% Evaluate model fit and return measure of goodness of fit (f).
%
% 8/1/05    dhb, pr     Wrote it.
% 8/2/07    dhb         Get rid of silly call to ComputeNakaRushtonError.
% 11/21/12  JPW         Made personalized function.
% 8/14      JPW         Added input variable to specify type of error.
% 2015      JPW         Added Negative Binomial Error

if nargin < 5
    disp('Must specify type of error to be used...')
    return
elseif (nargin < 6)
    showplot = 0;
end


prediction = ComputeNakaRushtonJPW(params,contrast,FITSTR);

if strcmp(error,'Poisson') || strcmp(error,'poisson')
    f = sum(prediction)-(log(prediction')*response);  % -1 * log-likelihood (Poisson)
elseif strcmp(error,'Gaussian') || strcmp(error,'gaussian')
    f = sum((prediction-response).^2); % Sum of squarred error (Gaussian)
elseif strcmp(error,'Bernoulli') || strcmp(error,'bernoulli')
    f = -1 * (log(prediction')*response + (log(1-prediction)'*(1-response))); % Bernoulli error
elseif strcmp(error,'NegativeBinomial') || strcmp(error,'negativebinomial')
    %f = -sum(gammaln(r+x)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r)); % from nbinfit
    %f = -sum(gammaln(r+x)) + sum(gammaln(r)) - (log(p)'*response) - (log(1-p)'*r); % from wikipedia
    
    kappa = params(end);
    mu = prediction;
    sigsq = mu + kappa * mu.^2;
    p = (sigsq - mu) ./ sigsq;
    r = mu.^2 ./ (sigsq - mu);
    if kappa == 0
        f = sum(prediction)-(log(prediction')*response);  % -1 * log-likelihood (Poisson)
    elseif any(r <= 0)
        disp('r <= 0')
        f = 5000000000;
        keyboard
    elseif any(p>1) || any(p<0)
        disp('p > 1 or p < 0')
        f = 5000000000;
        keyboard
    elseif any(~isreal(mu))
        disp('mu has immaginary component')
        f = 5000000000;
        keyboard
    else
        f = -sum(gammaln(r+response)) + sum(gammaln(r))...
            - (log(p)'*response) - (log(1-p)'*r); % from wikipedia
        %f = gammaln(r+response) + gammaln(r) - (log(p).*response) - (log(1-p).*r)
    end
end

if showplot == 1
    if strcmp(FITSTR,'asymmetric') 
        figure(500); clf; hold on; grid on;
        plot(contrast,response,'ko')
        plot(contrast,prediction,'m*');
        set(gca,'ylim',[0 max(prediction)],...
            'xlim',[min(contrast(:)) max(contrast(:))],...
            'ylim',[min(contrast(:)) max(contrast(:))]);
        drawnow
    end
    if strcmp(FITSTR,'surface7') || strcmp(FITSTR,'surface8') || strcmp(FITSTR,'conicsection')
        figure(500); 
        campos = get(gca,'cameraposition');
        clf; hold on; grid on;
        plot3(contrast(:,1),contrast(:,2),response,'ko')
        plot3(contrast(:,1),contrast(:,2),prediction,'m*')
        set(gca,'ylim',[0 max(prediction)],...
            'xlim',[min(contrast(:)) max(contrast(:))],...
            'ylim',[min(contrast(:)) max(contrast(:))],...
            'cameraposition',campos);
        drawnow
    end
end

% Handle bizarre parameter values.
if (isnan(f))
    f = 500000000;
    disp('LL is nan...');
    %keyboard
end

if (all(prediction == response))
    f = 500000000;
    disp('EXACT MATCH!');
    %keyboard
end

if ~isreal(f)
    f = 500000000;
    disp('Irrational -LL...')
    %keyboard
end

if any(~isreal(prediction))
    f = 500000000;
    disp('Irrational prediction...')
    %keyboard
end

end


