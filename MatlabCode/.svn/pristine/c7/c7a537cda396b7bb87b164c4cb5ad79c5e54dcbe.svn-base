function [beta, fitInfo] = WN_LARS(X, y, lassoon, preproc, nsteps)
%%
% This applys LASSO regularization using least angle regression (LARS)
%  
% Inputs:
%   'x' [nstim, ncovariates] - 2D array of stimuli data
%   'y' [nstim] - Vector of response

%   'preproc' [bool] - Center and scale input data
%   'nsteps' [int] - number of steps to compute (max iterations)
%
% Outputs:
%   'beta' [nweights, nsteps] - 2D array of weights for each step
%   'fitInfo' [struct] - Contains information about model goodness of fit
%       '.df' [nstep] - Model degrees of freedom for each step
%       '.cp' [nstep] - Model Cp confidence value. (Lower is better)
%% 

% If unspecified, set to run LASSO, preprocess data, and return all steps
if nargin < 3
    lassoon = true;
end
if nargin < 4
    preproc = true;
end
if nargin < 5
    nsteps = 512*size(X,2);
end

%% Preprocess data if necessary
% Ensure y is a column vector. Take out when working with higher dimensional responses.
if size(y, 1) == 1 && size(y,2) > 1; y = y(:); end

% Center and normalize response and input if requested
if preproc
    X = X-mean(X,1);
    Xc = zeros(size(X));
    for nn = 1:size(X,2)
    	Xc(:,nn) = X(:,nn)/norm(X(:,nn));
    end
    yc = y - mean(y);
else
    yc = y;
    Xc = X;
end

% Initialize weights vector as zeros
lassomod = false;

% Initialize beta
beta = zeros(size(X, 2), 1);

% Initialize mu (current position)
mu = zeros(size(X, 1), 1);

% Calculate OLS solution to use in confidence parameter calculation
bOLS = Xc\yc;

% Define possible covariates
inact = 1:size(X,2); % Inactive indices
act = []; % Active indices
stopcondition = false; % If while loop should be stopped
ii = 1; % Step Counter

%% Loop through steps 
wtbar = waitbar(0, 'Computing LARS...');
while ~stopcondition
    
    % Define correlations between each covariate and current residual vector 
    c = Xc'*(yc-mu);    
    
    % Update active and inactive set
    % If lasso condition is met remove idx from active indices
    if lassomod
        % Update active and inactive set by removing 
        inact = [inact act(gammaidx)];
        act(gammaidx) = [];
        
        % Pick direction to perform initial search in
        [C, actidx] = max(abs(c(inact)));
        lassomod = false;
    else
        % Pick direction to perform initial search in
        [C, actidx] = max(abs(c(inact)));
        
        % Update active and inactive set
        act = [act inact(actidx)];
        inact(inact == act(end)) = [];
    end
    
    % Calculate grammian
    G_A = Xc(:,act)'*Xc(:,act);    
    
    % Get signs of active correlations
    s = sign(c);
    
    % Compute inverse of G_A
    G_Achol = chol(G_A);
    G_Ainv = (G_Achol\(G_Achol'\eye(size(G_A,1))))*s(act);
    
    % Compute direction of next step
    A_A = sum(G_Ainv.*s(act),'all').^(-0.5);
    w = sum(A_A*G_Ainv,2);
    u_A = Xc(:,act)*w;
    
    % Calculate step size
    a = Xc'*u_A;    
    gamma_temp = [(C-c(inact))./(A_A-a(inact));(C+c(inact))./(A_A+a(inact)); C/A_A];
    gamma = min(gamma_temp(gamma_temp>0));
    
    % LASSO modification to LARS
    if lassoon
        gamma_lasso = -beta(act,ii)./w;
        [gamma_tilde] = min(gamma_lasso(gamma_lasso>0));
        if gamma_tilde < gamma
            lassomod = true;
            gamma = gamma_tilde;
            gammaidx = find(gamma_lasso == gamma_tilde);
        end
    end
    
    % Update mu
    mu = mu + gamma*u_A;
    
    % Update beta
    beta(act,ii+1) = beta(act, ii) +  gamma*w;
    
    % Calculate model goodness of fit parameters
    % Define degree of freedom for step
    fitInfo.df(ii) = length(nonzeros(beta(:,ii+1)));
    % Calculate cp - risk value (models with a higher value are more likely
    % to be poor)
    fitInfo.cp(ii) = sum((yc-mu).^2)/var(yc-Xc*bOLS) - length(mu) + 2*fitInfo.df(ii);
    % Calculate pe - percent explained (unsure if this is correct)
%     fitInfo.pe(ii) = 1 - sum((mu-yc).^2)/sum(yc.^2);
    
    % Check stop condition
    if size(beta, 2) == nsteps || length(nonzeros(beta(:,ii+1))) == size(X,2)
        stopcondition = true;
    end
    
    % Display progress 
    waitbar(length(nonzeros(beta(:,ii)))./(size(X,2)), wtbar);
    
    ii = ii + 1;
end
close(wtbar);
% Trim beta
% % beta = beta(:,2:end);

