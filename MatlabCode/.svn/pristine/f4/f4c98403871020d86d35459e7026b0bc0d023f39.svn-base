function [params,nLL,hessval] = GLMSGUI_AUC(sub,surftype,paramvals)
global GLMP

%Preallocate space and set up variables
errortype = 'Bernoulli';
xvals = GLMP.subunit{sub}.uniqueLcc;
yvals = GLMP.subunit{sub}.uniqueMcc;
zvals = GLMP.subunit{sub}.AUC;
Aguess = 1;
sigguess = [.04 0.4];
expguess = 2;
blguess = .5;
ang = paramvals(end-1);
nLL = Inf;
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','UseParallel',1,'Hessian','bfgs',...
    'display','off');

% Fit surface according to which surf best fit the neural data, preserving preferred direction
if strcmp(surftype,'surface7')
    
    %Set up some variables
    vlb = [.5  .001  .001  .001  .5   -pi];
    vub = [1     10    10   10   .5    pi];
    guessIdx = fullfact([numel(sigguess) numel(sigguess)]);
    
    % Roll through all guesses
    for n = 1:size(guessIdx,1)
        
        % Generating an initial guess
        params0 = [Aguess sigguess(guessIdx(n,1)) sigguess(guessIdx(n,2))...
            expguess blguess ang];
        
        % Fit all rotated points
        [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,surftype,errortype);
        if fval < nLL
            params = f1;
            nLL = fval;
            hessval = hess;
        end
    end    
    
elseif strcmp(surftype,'surface8')
    
    %Set up some variables
    vlb = [.5 0.001 0.001 0.001 0.001 .001  .5  -pi];
    vub = [1   10    10    10    10    10   .5   pi];
    guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(sigguess) numel(sigguess)]);
   
    % Roll through all guesses
    for n = 1:size(guessIdx,1)
        
        % Generating an initial guess
        params0 = [Aguess sigguess(guessIdx(n,1)) sigguess(guessIdx(n,2))...
            sigguess(guessIdx(n,3)) sigguess(guessIdx(n,4))...
            expguess blguess ang];
        
        % Fit all rotated points
        [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,surftype,errortype);
        if fval < nLL
            nLL = fval;
            params = f1;
            hessval = hess;
        end
    end    
    
end



