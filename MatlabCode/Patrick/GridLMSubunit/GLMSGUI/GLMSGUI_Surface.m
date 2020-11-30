function [] = GLMSGUI_Surface(~,~,varargin)
% This code fits surfaces to GLMP datasets.

SetUpFig()

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% For running the analysis
if ~isempty(varargin) 
    if numel(varargin) == 2
        conpanel.subselect = varargin{2};
        set(conpanel.subbuttons,'selectedobject',conpanel.(['sub' num2str(conpanel.subselect)]))
        
        %Save variables
        set(surffig.conpanel,'UserData',conpanel);
        set(surffig.fitspanel,'UserData',fitspanel);
        set(gcf,'UserData',surffig);
        
    end
    LoadPreviousFits(varargin{1})
end

end


%%% Interactive functions %%%

function startanalysis(~,~)
global GLMP
% This is the heart of the analysis. Here, we will route the analysis to
% different subfunctions.

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');

disp('')
disp(['Beginning Surface Analysis on ' GLMP.datafile...
    '(sub # ' num2str(conpanel.subselect) ')']);

% Run through fits and compare
OneDFit()
TwoDFit()
CheckOneDFit()
CheckTwoDFit()
LogLikelihoodRatioTest()
NormalizedLogLikelihoodTest()

disp('Surface Analysis Completed.')

end

function Subsel(~,eventdata)
global GLMP
% Set up subunit selection

% Load variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

subval = get(eventdata.NewValue,'string');
if strcmp(subval,'Subunit #1')
    conpanel.subselect = 1;
elseif strcmp(subval,'Subunit #2')
    conpanel.subselect = 2;
elseif strcmp(subval,'Subunit #3')
    conpanel.subselect = 3;
end

% Replot
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
axes(fitspanel.axes.oneDfit); cla; hold on; grid on;
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
axes(fitspanel.axes.twoDfit); cla; hold on; grid on;
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end

% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(gcf,'UserData',surffig);

end


%%% Analyses %%%

function OneDFit()
global GLMP
disp('Fitting a 1D function to data...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load data parameters
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Set bounds for 1D fit
vub = fitspanel.ub.oneD;
vlb = fitspanel.lb.oneD;

% Set initial guesses and options for 1D fit
angs = unique(GLMP.subunit{sub}.theta);
Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
sigguess = max(GLMP.subunit{sub}.rho)/2;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end


%%% Bidirectional Fits %%%

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
bidnLL = Inf;
for rot = 1:size(angs)
    
    % New initial guess
    paramsGuess = [Aguess 1/sigguess 1/sigguess 0 expguess blguess angs(rot) kappaguess];
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
        
    % If better than previous best fit, replace parmeter values
    if fval < bidnLL
        bidparams = f1;
        bidnLL = fval;
        bidhess = hess;
    end
end


%%% Unidirectional Fits %%%

% Set sig2 bounds to Inf;
vlb0 = vlb;
vlb0(3) = 0;
vub0 = vub;
vub0(3) = 0;
unidnLL = Inf;
for rot = 1:size(angs)
    
    % New initial guess
    paramsGuess = [Aguess 1/sigguess 0 0 expguess blguess angs(rot) kappaguess];
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb0,vub0,[],fitspanel.options,...
        [Lcc Mcc],nsp,fitspanel.surftype,fitspanel.errortype);
        
    % If better than previous best fit, replace parmeter values
    if fval < unidnLL
        unidparams = f1;
        unidnLL = fval;
        unidhess = hess;
    end
end

% Test special cases
% Unidirectional surface (responsive in only 1 direction) is the base case. (6 params)
% Suppression and excitation surfaces (7 params) along the second axis are 
% each tested against this case. P-vals are compared against a threshold.

%%% 1D version of best 2D fit %%%
if ~isempty(fitspanel.params.twoD)
    paramsGuess = fitspanel.params.twoD;
    paramsGuess(4) = 0;
    vlb0 = vlb;
    vub0 = vub;
      
    % If sig2 is 0, guess is unidirectional
    if paramsGuess(3) == 0
        vlb0(3) = 0;
        vub0(3) = 0;
    end

    % Test 2D fit
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb0,vub0,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
    if f1(3) == 0
        if fval < unidnLL
            unidnLL = fval;
            unidparams = f1;
            unidhess = hess;
        end
    else
        if fval < bidnLL
            bidnLL = fval;
            bidparams = f1;
            bidhess = hess;
        end
    end
end

% Testing for suppression
% paramsGuess = bidparams;
% paramsGuess(3) = -.5;
% [supparams,supnLL,~,~,~,~,suphess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
%     [],[],[],[],vlb0,vub0,[],fitspanel.options,[Lcc Mcc],nsp,...
%     fitspanel.surftype,fitspanel.errortype);

% Testing unidirectional case
% paramsGuess = bidparams;
% paramsGuess(3) = 0;
% vlb0 = vlb;
% vub0 = vub;
% vlb0(3) = 0;
% vub0(3) = 0;
% [unidparams,unidnLL,~,~,~,~,unidhess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
%     a,b,[],[],vlb0,vub0,[],fitspanel.options,[Lcc Mcc],nsp,...
%     fitspanel.surftype,fitspanel.errortype);


%%% Custom guess %%%
if any(str2double(get(conpanel.paramvals.oneD.A,'string')))
    %disp('Making a 1D custom guess...')
    paramsGuess(1) = str2double(get(conpanel.paramvals.oneD.A,'string'));
    paramsGuess(2) = 1/str2double(get(conpanel.paramvals.oneD.sig1,'string'));
    paramsGuess(3) = 1/str2double(get(conpanel.paramvals.oneD.sig2,'string'));
    paramsGuess(4) = 0;
    paramsGuess(5) = str2double(get(conpanel.paramvals.oneD.exp,'string'));
    paramsGuess(6) = str2double(get(conpanel.paramvals.oneD.bl,'string'));
    paramsGuess(7) = str2double(get(conpanel.paramvals.oneD.rot,'string'))/180*pi;
    paramsGuess(8) = str2double(get(conpanel.paramvals.oneD.kappa,'string'));
    vub0 = vub;
    vlb0 = vlb;
    
    % If second axis is 0, guess is unidirectional
    if paramsGuess(3) == 0
        %disp('Custom guess is unidirectional.')
        vlb0(3) = 0;
        vub0(3) = 0;
    else
        %disp('Custom guess is bidirectional.')
    end
    
    % Fit 1D function
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb0,vub0,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
    if f1(3) == 0
        if fval < unidnLL
            unidnLL = fval;
            unidparams = f1;
            unidhess = hess;
        end
    else
        if fval < bidnLL
            bidnLL = fval;
            bidparams = f1;
            bidhess = hess;
        end
    end
end


%%% Compare unidirectional with bidirectional fits %%%
diffLL = 2 * (unidnLL - bidnLL); % larger model nLL - smaller model LL
bidpval = 1 - chi2cdf(diffLL,1);
% if bidpval < fitspanel.thresh
%     disp('1D Surface is bidirectional.')
% else
%     disp('1D Surface is unidirectional.')
% end


%%% Record and display results %%%

% Record LL, pvals, params, hess for each fit
nLLs = cat(1,unidnLL,bidnLL);
pvals = cat(1,0,bidpval);
parvals = cat(1,unidparams,bidparams);
hessvals{1} = unidhess;
hessvals{2} = bidhess;
names{1} = 'unidirectional';
names{2} = 'bidirectional';
fitspanel.fitcomps.oneD.names = names;
fitspanel.fitcomps.oneD.nLLs = nLLs;
fitspanel.fitcomps.oneD.pvals = pvals;
fitspanel.fitcomps.oneD.parvals = parvals;
fitspanel.fitcomps.oneD.hess = hessvals;

% Find surface with greatest significance
%sigidx = find(pvals < fitspanel.thresh);
%[~,idx] = min(nLLs(sigidx));
%whichsurf = sigidx(idx);
whichsurf = 2; % just chosing bichromatic surface for now
params = parvals(whichsurf,:);
nLL = nLLs(whichsurf);
hessval = hessvals{whichsurf};

% Record mean R-squared 
predresp = ComputeNakaRushtonJPW(params,[Lcc Mcc],fitspanel.surftype);
sstot = sum((nsp - mean(nsp)).^2);
ssres = sum((nsp - predresp).^2);
fitspanel.rsq.oneD = 1 - (ssres/sstot);

% Record Hessian std
std = sqrt(1/hessval(end-1,end-1))/pi*180;
set(conpanel.paramvals.oneD.std,'string',std);

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.oneDfit); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
fitspanel.surf.oneD = surfc(xx,yy,surface);
set(fitspanel.surf.oneD(1),'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);

% Fill in parmeter values in figure
params1 = round(params,2);
set(conpanel.paramvals.oneD.sig1,'string',round(1/params1(2),3));
set(conpanel.paramvals.oneD.sig2,'string',round(1/params1(3),3));
set(conpanel.paramvals.oneD.A,'string',params1(1));
set(conpanel.paramvals.oneD.exp,'string',params1(5));
set(conpanel.paramvals.oneD.bl,'string',params1(6));
set(conpanel.paramvals.oneD.rot,'string',round(params1(7)/pi*180,0));
set(conpanel.paramvals.oneD.kappa,'string',params1(8));
set(conpanel.paramvals.oneD.LL,'string',-nLL);

% Save results
fitspanel.LL.oneD = -nLL;
fitspanel.params.oneD = params;
fitspanel.Hessian.oneD = hessval;
fitspanel.std.oneD = std;

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',surffig)

end

function TwoDFit()
global GLMP
disp('Fitting a 2D function to data...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load data variables
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Upper and lower bounds on fit parameters
vub = fitspanel.ub.twoD;
vlb = fitspanel.lb.twoD;

% Parameter guesses
Aguess = max(GLMP.nspikes)*.8;
sigguess = max(GLMP.subunit{sub}.rho)/2;
orthsigguess = [-sigguess Inf sigguess];
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
angs = linspace(vlb(end-1),vub(end-1),7)';
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < 0
    kappaguess = 0;
end
guessIdx = fullfact([numel(Aguess) 1 numel(sigguess) numel(orthsigguess) numel(expguess) numel(angs)]);

%%% Bidirectional fits %%%
bidnLL = Inf;
warning('off','MATLAB:singularMatrix'); 
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    paramsGuess = [Aguess(guessIdx(rot,1)) 1/sigguess 1/sigguess(guessIdx(rot,3))...
        1/orthsigguess(guessIdx(rot,4)) expguess(guessIdx(rot,5))... 
        blguess angs(guessIdx(rot,6)) kappaguess];
    
    % Using paramsGuess as initial guess into fit
    warning off MATLAB:nearlySingularMatrix
    warning off MATLAB:illConditionedMatrix
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);

    % If -LL is less than current best fit, replace it
    if fval < bidnLL
        bidparams = f1;
        bidnLL = fval;
        bidhess = hess;
    end
end

%%% Unidirectional fits %%%
unidnLL = Inf;
vub0 = vub;
vlb0 = vlb;
vub0(3) = 0;
vlb0(3) = 0;
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    paramsGuess = [Aguess(guessIdx(rot,1)) 1/sigguess(guessIdx(rot,2)) 0 ...
        1/orthsigguess(guessIdx(rot,4)) expguess(guessIdx(rot,5))... 
        blguess angs(guessIdx(rot,6)) kappaguess];
    
    % Using paramsGuess as initial guess into fit
    warning off MATLAB:nearlySingularMatrix
    warning off MATLAB:illConditionedMatrix
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb0,vub0,[],fitspanel.options,...
        [Lcc Mcc],nsp,fitspanel.surftype,fitspanel.errortype);

    % If -LL is less than current best fit, replace it
    if fval < unidnLL
        unidparams = f1;
        unidnLL = fval;
        unidhess = hess;
    end
end


%%% 2D version of best 1D fit %%%
if ~isempty(fitspanel.params.oneD)
    paramsGuess = fitspanel.params.oneD;
    vlb0 = vlb;
    vub0 = vub;
      
    % If sig2 is 0, guess is unidirectional
    if paramsGuess(3) == 0
        vlb0(3) = 0;
        vub0(3) = 0;
    end

    % Test 1D fit
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb0,vub0,[],fitspanel.options,...
        [Lcc Mcc],nsp,fitspanel.surftype,fitspanel.errortype);
    if f1(3) == 0
        if fval < unidnLL
            unidnLL = fval;
            unidparams = f1;
            unidhess = hess;
        end
    else
        if fval < bidnLL
            bidnLL = fval;
            bidparams = f1;
            bidhess = hess;
        end
    end
end


%%% Custom guess %%%
if any(str2double(get(conpanel.paramvals.twoD.A,'string')))
    %disp('Making a 2D custom guess...')
    paramsGuess(1) = str2double(get(conpanel.paramvals.twoD.A,'string'));
    paramsGuess(2) = 1/str2double(get(conpanel.paramvals.twoD.sig1,'string'));
    paramsGuess(3) = 1/str2double(get(conpanel.paramvals.twoD.sig2,'string'));
    paramsGuess(4) = 1/str2double(get(conpanel.paramvals.twoD.orthsig,'string'));
    paramsGuess(5) = str2double(get(conpanel.paramvals.twoD.exp,'string'));
    paramsGuess(6) = str2double(get(conpanel.paramvals.twoD.bl,'string'));
    paramsGuess(7) = str2double(get(conpanel.paramvals.twoD.rot,'string'))/180*pi;
    paramsGuess(8) = str2double(get(conpanel.paramvals.twoD.kappa,'string'));
    vub0 = vub;
    vlb0 = vlb;
    
    % If second axis is 0, guess is unidirectional
    if paramsGuess(3) == 0
        vlb0(3) = 0;
        vub0(3) = 0;
    end
    
    % Fit 1D function
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        fitspanel.a,fitspanel.b,[],[],vlb0,vub0,[],fitspanel.options,[Lcc Mcc],nsp,...
        fitspanel.surftype,fitspanel.errortype);
    if f1(3) == 0
        if fval < unidnLL
            unidnLL = fval;
            unidparams = f1;
            unidhess = hess;
        end
    else
        if fval < bidnLL
            bidnLL = fval;
            bidparams = f1;
            bidhess = hess;
        end
    end
end


%%% Compare unidirectional with bidirectional fits %%%
diffLL = unidnLL - bidnLL; % larger model nLL - smaller model LL
bidpval = 1 - chi2cdf(2*diffLL,1);
if bidpval < fitspanel.thresh
    disp('1D Surface is bidirectional.')
else
    disp('1D Surface is unidirectional.')
end

%%% Record and display results %%%

% Record LL, pvals, params, hess for each fit
nLLs = cat(1,unidnLL,bidnLL);
pvals = cat(1,0,bidpval);
parvals = cat(1,unidparams,bidparams);
hessvals{1} = unidhess;
hessvals{2} = bidhess;
names{1} = 'unidirectional';
names{2} = 'bidirectional';
fitspanel.fitcomps.twoD.names = names;
fitspanel.fitcomps.twoD.nLLs = nLLs;
fitspanel.fitcomps.twoD.pvals = pvals;
fitspanel.fitcomps.twoD.parvals = parvals;
fitspanel.fitcomps.twoD.hess = hessvals;

% Find surface with greatest significance
% sigidx = find(pvals < fitspanel.thresh);
% [~,idx] = min(nLLs(sigidx));
% whichsurf = sigidx(idx);
whichsurf = 2; % just chosing bichromatic surface for now
params = parvals(whichsurf,:);
nLL = nLLs(whichsurf);
hessval = hessvals{whichsurf};

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.twoDfit); cla; hold on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
fitspanel.surf.twoD = surfc(xx,yy,surface);
set(fitspanel.surf.twoD(1),'EdgeColor','none')
alpha(.5);
set(gca,'view',[0 90])

%Fill in parameter values
params1 = round(params,2);
set(conpanel.paramvals.twoD.A,'string',params1(1));
set(conpanel.paramvals.twoD.sig1,'string',round(1/params1(2),3));
set(conpanel.paramvals.twoD.sig2,'string',round(1/params1(3),3));
set(conpanel.paramvals.twoD.orthsig,'string',round(1/params1(4),3));
set(conpanel.paramvals.twoD.exp,'string',params1(5));
set(conpanel.paramvals.twoD.bl,'string',params1(6));
set(conpanel.paramvals.twoD.rot,'string',round(params1(7)/pi*180,0));
set(conpanel.paramvals.twoD.kappa,'string',params1(8));

% Save parameter values
fitspanel.params.twoD = params1;
fitspanel.LL.twoD = -nLL;
set(conpanel.paramvals.twoD.LL,'string',round(-nLL,1));

% Record mean R-squared 
predresp = ComputeNakaRushtonJPW(params1,[Lcc Mcc],fitspanel.surftype);
sstot = sum((nsp - mean(nsp)).^2);
ssres = sum((nsp - predresp).^2);
fitspanel.rsq.twoD = 1 - (ssres/sstot);

% Determine standard deviation
fitspanel.Hessian.twoD = hessval;
std = sqrt(1/hessval(end-1,end-1))/pi*180;
fitspanel.std.twoD = std;
set(conpanel.paramvals.twoD.std,'string',round(std,1));

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',surffig)

end

function CheckOneDFit()
global GLMP
disp('Rechecking 1D fit...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load data variables
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
unidnLL = fitspanel.fitcomps.oneD.nLLs(1);
bidnLL = fitspanel.fitcomps.oneD.nLLs(2);
unidparams = fitspanel.fitcomps.oneD.parvals(1,:);
bidparams = fitspanel.fitcomps.oneD.parvals(2,:);
unidhess = fitspanel.fitcomps.oneD.hess{1};
bidhess = fitspanel.fitcomps.oneD.hess{2};

% Upper and lower bounds on fit parameters
vub = fitspanel.ub.oneD;
vlb = fitspanel.lb.oneD;

% Use 2D fit as 1D guess
paramsGuess = fitspanel.params.twoD;
paramsGuess(4) = 0;

% Which surface type is being suggested? Interpret second sigma as indicator.
% If second axis is 0, guess is unidirectional
if paramsGuess(3) == 0
    vlb(3) = 0;
    vub(3) = 0;
end

% Fit 1D function
[f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
    fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,...
    [Lcc Mcc],nsp,fitspanel.surftype,fitspanel.errortype);
if f1(3) == 0
    if fval < unidnLL
        unidnLL = fval;
        unidparams = f1;
        unidhess = hess;
    end
else
    if fval < bidnLL
        bidnLL = fval;
        bidparams = f1;
        bidhess = hess;
    end
end

% Compare unidirectional with bidirectional fit
diffLL = 2 * (unidnLL - bidnLL); % larger model nLL - smaller model LL
bidpval = 1 - chi2cdf(diffLL,1);
if bidpval < fitspanel.thresh
    disp('1D Surface is bidirectional.')
else
    disp('1D Surface is unidirectional.')
end

% Collect all unid and bid params
nLLs = cat(1,unidnLL,bidnLL);
pvals = cat(1,0,bidpval);
parvals = cat(1,unidparams,bidparams);
hessvals{1} = unidhess;
hessvals{2} = bidhess;
fitspanel.fitcomps.oneD.nLLs = nLLs;
fitspanel.fitcomps.oneD.pvals = pvals;
fitspanel.fitcomps.oneD.parvals = parvals;
fitspanel.fitcomps.oneD.hess = hessvals;

% Find smallest significant surface
% sigidx = find(pvals < fitspanel.thresh);
% [~,idx] = min(nLLs(sigidx));
% whichsurf = sigidx(idx);
whichsurf = 2; % just chosing bichromatic surface for now
params = parvals(whichsurf,:);
nLL = nLLs(whichsurf);
hessval = hessvals{whichsurf};

% Record mean R-squared 
predresp = ComputeNakaRushtonJPW(params,[Lcc Mcc],fitspanel.surftype);
sstot = sum((nsp - mean(nsp)).^2);
ssres = sum((nsp - predresp).^2);
fitspanel.rsq.oneD = 1 - (ssres/sstot);

% Determine standard deviation
fitspanel.Hessian.oneD = hessval;
std = sqrt(1/hessval(end-1,end-1))/pi*180;
fitspanel.std.oneD = std;
set(conpanel.paramvals.oneD.std,'string',round(std,1));

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.oneDfit); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
fitspanel.surf.oneD = surfc(xx,yy,surface);
set(fitspanel.surf.oneD(1),'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);

% Fill in parmeter values in figure
params1 = round(params,2);
set(conpanel.paramvals.oneD.A,'string',params1(1));
set(conpanel.paramvals.oneD.sig1,'string',round(1/params1(2),3));
set(conpanel.paramvals.oneD.sig2,'string',round(1/params1(3),3));
set(conpanel.paramvals.oneD.exp,'string',params1(5));
set(conpanel.paramvals.oneD.bl,'string',params1(6));
set(conpanel.paramvals.oneD.rot,'string',round(params1(7)/pi*180,0));
set(conpanel.paramvals.oneD.kappa,'string',params1(8));
set(conpanel.paramvals.oneD.LL,'string',round(-nLL,1));

% Save results
fitspanel.LL.oneD = -nLL;
fitspanel.params.oneD = params;
    
% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',surffig)

end

function CheckTwoDFit()
global GLMP
disp('Rechecking 2D fit...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load data variables
sub = conpanel.subselect;
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
unidnLL = fitspanel.fitcomps.twoD.nLLs(1);
bidnLL = fitspanel.fitcomps.twoD.nLLs(2);
unidparams = fitspanel.fitcomps.twoD.parvals(1,:);
bidparams = fitspanel.fitcomps.twoD.parvals(2,:);
unidhess = fitspanel.fitcomps.twoD.hess{1};
bidhess = fitspanel.fitcomps.twoD.hess{2};

% Upper and lower bounds on fit parameters
vub = fitspanel.ub.twoD;
vlb = fitspanel.lb.twoD;

% Use 2D fit as 1D guess
paramsGuess = fitspanel.params.oneD;

% Which surface type is being suggested? Interpret second sigma as indicator.
% If second axis is 0, guess is unidirectional
if paramsGuess(3) == 0
    vlb(3) = 0;
    vub(3) = 0;
end

% Fit using 1D parameters
[f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
    fitspanel.a,fitspanel.b,[],[],vlb,vub,[],fitspanel.options,...
    [Lcc Mcc],nsp,fitspanel.surftype,fitspanel.errortype);
if f1(3) == 0
    if fval < unidnLL
        unidnLL = fval;
        unidparams = f1;
        unidhess = hess;
    end
else
    if fval < bidnLL
        bidnLL = fval;
        bidparams = f1;
        bidhess = hess;
    end
end

% Compare unidirectional with bidirectional fit
diffLL = 2 * (unidnLL - bidnLL); % larger model nLL - smaller model LL
bidpval = 1 - chi2cdf(diffLL,1);
% if bidpval < fitspanel.thresh
%     disp('2D Surface is bidirectional.')
% else
%     disp('2D Surface is unidirectional.')
% end


% Collect all unid and bid params
nLLs = cat(1,unidnLL,bidnLL);
pvals = cat(1,0,bidpval);
parvals = cat(1,unidparams,bidparams);
hessvals{1} = unidhess;
hessvals{2} = bidhess;
fitspanel.fitcomps.twoD.nLLs = nLLs;
fitspanel.fitcomps.twoD.pvals = pvals;
fitspanel.fitcomps.twoD.parvals = parvals;
fitspanel.fitcomps.twoD.hess = hessvals;

% Find smallest significant surface
%sigidx = find(pvals < fitspanel.thresh);
%[~,idx] = min(nLLs(sigidx));
%whichsurf = sigidx(idx);
whichsurf = 2; % Just chosing bichromatic fit for now
params = parvals(whichsurf,:);
nLL = nLLs(whichsurf);
hessval = hessvals{whichsurf};

% Record mean R-squared 
predresp = ComputeNakaRushtonJPW(params,[Lcc Mcc],fitspanel.surftype);
sstot = sum((nsp - mean(nsp)).^2);
ssres = sum((nsp - predresp).^2);
fitspanel.rsq.twoD = 1 - (ssres/sstot);

% Determine standard deviation
fitspanel.Hessian.twoD = hessval;
std = sqrt(1/hessval(end-1,end-1))/pi*180;
fitspanel.std.twoD = std;
set(conpanel.paramvals.twoD.std,'string',round(std,1));

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(fitspanel.axes.twoDfit); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
fitspanel.surf.twoD = surfc(xx,yy,surface);
set(fitspanel.surf.twoD(1),'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);

% Fill in parmeter values in figure
params1 = round(params,2);
set(conpanel.paramvals.twoD.A,'string',params1(1));
set(conpanel.paramvals.twoD.sig1,'string',round(1/params1(2),3));
set(conpanel.paramvals.twoD.sig2,'string',round(1/params1(3),3));
set(conpanel.paramvals.twoD.exp,'string',params1(5));
set(conpanel.paramvals.twoD.bl,'string',params1(6));
set(conpanel.paramvals.twoD.rot,'string',round(params1(7)/pi*180,0));
set(conpanel.paramvals.twoD.kappa,'string',params1(8));
set(conpanel.paramvals.twoD.LL,'string',round(-nLL,1));

% Save results
fitspanel.LL.twoD = -nLL;
fitspanel.params.twoD = params;
    
% Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',surffig)

end

function LogLikelihoodRatioTest()

disp('Performing Log Likelihood Ratio Test...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Prepare parameters
nfp2 = numel(fitspanel.params.twoD);
nfp1 = nfp2 - 1;
LL1 = fitspanel.LL.oneD;
LL2 = fitspanel.LL.twoD;

% Actual statistical test
if LL2 < LL1
    set(conpanel.paramvals.pval,'string','nan')
    fitspanel.LLRTpval = nan;
    disp('1D Likelihood greater than 2D likelihood!')
else
    diffLL = 2 * (LL2 - LL1);
    pval = 1 - chi2cdf(diffLL,nfp2-nfp1);
    set(conpanel.paramvals.pval,'string',pval);
    fitspanel.LLRTpval = pval;
end

% Save User Variables
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.conpanel,'UserData',conpanel);
set(gcf,'UserData',surffig);

end

function NormalizedLogLikelihoodTest()
global GLMP
% Normalized log likelhood between mean and exact fit

disp('Calculating Normalized Log Likelihood...')

% Load Figure Variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load data variables
sub = conpanel.subselect;
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
kappaguess = fitspanel.params.oneD(end);

% Find minimum LL (using Negative Binomial)
mu = repmat(mean(nsp),size(nsp));
x = nsp;
[kappa,f] = fmincon('minimizeKappa',kappaguess,[],[],[],[],0,[],[],fitspanel.options,mu,x);
fitspanel.normLL.lbLL = -f;
fitspanel.normLL.lbkappa = kappa;

% Find maximum LL (using negative binomial)
kappaguess = fitspanel.params.twoD(end);
nsp = [];
prednsp = [];
for t = 1:numel(GLMP.subunit{sub}.uniqueIdx)
    nsp = cat(1,nsp,GLMP.subunit{sub}.nspikes(GLMP.subunit{sub}.uniqueIdx{t}));
    prednsp = cat(1,prednsp,repmat(GLMP.subunit{sub}.meannspikes(t),numel(GLMP.subunit{sub}.uniqueIdx{t}),1));
end
meanbl = repmat(mean(GLMP.subunit{sub}.blnspikes),numel(GLMP.subunit{sub}.blnspikes),1);
prednsp = cat(1,prednsp,meanbl);
nsp = cat(1,nsp,GLMP.subunit{sub}.blnspikes);
mu = prednsp;
x = nsp;
[kappa,f] = fmincon('minimizeKappa',kappaguess,[],[],[],[],0,[],[],fitspanel.options,mu,x);
fitspanel.normLL.ubLL = -f;
fitspanel.normLL.ubkappa = kappa;

% Calculate normalized LL values
ub = fitspanel.normLL.ubLL;
lb = fitspanel.normLL.lbLL;
fitspanel.fitcomps.oneD.normLL = (-fitspanel.fitcomps.oneD.nLLs - lb)/(ub-lb);
fitspanel.fitcomps.twoD.normLL = (-fitspanel.fitcomps.twoD.nLLs - lb)/(ub-lb);
oneDLL = fitspanel.LL.oneD;
twoDLL = fitspanel.LL.twoD;
fitspanel.normLL.oneDnormLL = (oneDLL-lb)/(ub-lb);
fitspanel.normLL.twoDnormLL = (twoDLL-lb)/(ub-lb);

% Fill in fields
set(conpanel.paramvals.oneD.normLL,'string',round(fitspanel.normLL.oneDnormLL,2))
set(conpanel.paramvals.twoD.normLL,'string',round(fitspanel.normLL.twoDnormLL,2))
set(conpanel.paramvals.normLL,'string',round(fitspanel.normLL.twoDnormLL-fitspanel.normLL.oneDnormLL,4))

% Save User Variables
set(surffig.fitspanel,'UserData',fitspanel);
set(surffig.conpanel,'UserData',conpanel);
set(gcf,'UserData',surffig);

end


%%% Setup %%%

function SetUpFig()
global GLMP

% Set up figure
figure(60); clf;
set(gcf,'units','pixels','pos',[300 150 700 700],'NumberTitle','off',...
    'Name',['Surface Analyses (' GLMP.datafile ')']);

% Set up panels
surffig = get(gcf,'UserData');
surffig.fitspanel = uipanel('Pos',[.01 .5 .98 .49],'Parent',gcf,'title','Fits to Data');
surffig.conpanel = uipanel('pos',[.01 .01 .98 .48],'parent',gcf);

% Set up control panel
conpanel.startanalysis = uicontrol('style','pushbutton',...
    'parent',surffig.conpanel,'units','normalized','pos',[.8 .05 .19 .1],...
    'string','Start Analysis','fontsize',12,'callback',@startanalysis,...
    'backgroundcolor',[.2 1 0]);
conpanel.subbuttons = uibuttongroup('Parent',surffig.conpanel,...
    'units','normalized','pos',[.8 .65 .19 .31],...
    'SelectionChangeFcn',@Subsel,'title','Subunit Selection');
conpanel.sub1 = uicontrol('style','radiobutton','units','normalized',...
    'parent',conpanel.subbuttons,'pos',[.1 .675 .8 .25],...
    'string','Subunit #1','fontsize',12);
conpanel.sub2 = uicontrol('style','radiobutton','units','normalized',...
    'parent',conpanel.subbuttons,'pos',[.1 .375 .8 .25],...
    'string','Subunit #2','fontsize',12);
conpanel.sub3 = uicontrol('style','radiobutton','units','normalized',...
    'parent',conpanel.subbuttons,'pos',[.1 .075 .8 .25],...
    'string','Subunit #3','fontsize',12);
if numel(GLMP.subunit) == 1
    set(conpanel.sub2,'enable','off')
    set(conpanel.sub3,'enable','off')
end

% Set up stats panel
conpanel.params1pan = uipanel('parent',surffig.conpanel,'units','normalized',...
    'pos',[.01 .01 .35 .98],'title','1D Parameters');
conpanel.params2pan = uipanel('parent',surffig.conpanel,'units','normalized',...
    'pos',[.37 .01 .42 .98],'title','2D Parameters');
conpanel.stats = uipanel('parent',surffig.conpanel,'units','normalized',...
    'pos',[.8 .18 .19 .45],'title','Statistics');

% 1D panel
conpanel.paramslabels.oneD.A = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .9 .29 .08],'fontsize',12,'string','A =');
conpanel.paramslabels.oneD.sig = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .8 .29 .08],'fontsize',12,'string','Sigma =');
conpanel.paramslabels.oneD.exp = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .7 .29 .08],'fontsize',12,'string','Exponent =');
conpanel.paramslabels.oneD.bl = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .6 .29 .08],'fontsize',12,'string','Baseline =');
conpanel.paramslabels.oneD.rot = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .5 .29 .08],'fontsize',12,'string','Rotation =');
conpanel.paramslabels.oneD.kappa = uicontrol(conpanel.params1pan,'style','text',...
    'units','normalized','pos',[.01 .4 .29 .08],'fontsize',12,'string','Kappa =');

conpanel.paramvals.oneD.A = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.4 .9 .25 .08],'fontsize',12);
conpanel.paramvals.oneD.sig1 = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.4 .8 .25 .08],'fontsize',12);
conpanel.paramvals.oneD.sig2 = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.7 .8 .25 .08],'fontsize',12);
conpanel.paramvals.oneD.exp = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.4 .7 .25 .08],'fontsize',12);
conpanel.paramvals.oneD.bl = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.4 .6 .25 .08],'fontsize',12);
conpanel.paramvals.oneD.rot = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.4 .5 .25 .08],'fontsize',12);
conpanel.paramvals.oneD.kappa = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.4 .4 .25 .08],'fontsize',12);

% 2D panel
conpanel.paramslabels.twoD.A = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .9 .2 .08],'fontsize',12,'string','A =');
conpanel.paramslabels.twoD.sig = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .8 .2 .08],'fontsize',12,'string','Sigma =');
conpanel.paramslabels.twoD.orthsig = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.63 .8 .15 .08],'fontsize',12,'string','Orthsig =');
conpanel.paramslabels.twoD.exp = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .7 .2 .08],'fontsize',12,'string','Exp =');
conpanel.paramslabels.twoD.bl = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .6 .24 .08],'fontsize',12,'string','Baseline =');
conpanel.paramslabels.twoD.rot = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .5 .2 .08],'fontsize',12,'string','Rotation =');
conpanel.paramslabels.twoD.kappa = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .4 .2 .08],'fontsize',12,'string','Kappa =');

conpanel.paramvals.twoD.A = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.25 .9 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.sig1 = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.25 .8 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.sig2 = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.44 .8 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.orthsig = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.82 .8 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.exp = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.25 .7 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.bl = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.25 .6 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.rot = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.25 .5 .15 .08],'fontsize',12);
conpanel.paramvals.twoD.kappa = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.25 .4 .15 .08],'fontsize',12);

% Logliklihood values
conpanel.paramslabels.oneD.LL = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .25 .55 .08],'fontsize',12,'string','Log Likelihood =');
conpanel.paramvals.oneD.LL = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.6 .25 .3 .08],'fontsize',12);
conpanel.paramslabels.twoD.LL = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .25 .45 .08],'fontsize',12,'string','Log Likelihood =');
conpanel.paramvals.twoD.LL = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.5 .25 .25 .08],'fontsize',12);

% Normalized Loglikelihood Values
conpanel.paramlabels.oneD.normLL = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .15 .55 .08],'fontsize',12,'string','Normalized LL =');
conpanel.paramvals.oneD.normLL = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.6 .15 .3 .08],'fontsize',12);
conpanel.paramslabels.twoD.normLL = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .15 .45 .08],'fontsize',12,'string','Normalized LL =');
conpanel.paramvals.twoD.normLL = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.5 .15 .25 .08],'fontsize',12);

% Confidence Intervals
conpanel.paramslabels.oneD.std = uicontrol('style','text','parent',conpanel.params1pan,...
    'units','normalized','pos',[.01 .05 .55 .08],'fontsize',12,'string','Fit Std =');
conpanel.paramvals.oneD.std = uicontrol('style','edit','parent',conpanel.params1pan,...
    'units','normalized','pos',[.6 .05 .3 .08],'fontsize',12);
conpanel.paramslabels.twoD.std = uicontrol('style','text','parent',conpanel.params2pan,...
    'units','normalized','pos',[.01 .05 .45 .08],'fontsize',12,'string','Fit Std =');
conpanel.paramvals.twoD.std = uicontrol('style','edit','parent',conpanel.params2pan,...
    'units','normalized','pos',[.5 .05 .25 .08],'fontsize',12);

% Stats Panel
conpanel.paramlabels.pval = uicontrol('style','text','parent',conpanel.stats,...
    'units','normalized','pos',[.01 .75 .98 .15],'fontsize',12,'string','P-value =');
conpanel.paramvals.pval = uicontrol('style','edit','parent',conpanel.stats,...
    'units','normalized','pos',[.1 .55 .8 .2],'fontsize',12);
conpanel.paramlabels.normLL = uicontrol('style','text','parent',conpanel.stats,...
    'units','normalized','pos',[.01 .35 .98 .15],'fontsize',12,'string','normLL Diff =');
conpanel.paramvals.normLL = uicontrol('style','edit','parent',conpanel.stats,...
    'units','normalized','pos',[.1 .15 .8 .2],'fontsize',12);

% Set up fits panel (for displaying surfaces)
fitspanel.axes.oneDfit = axes('parent',surffig.fitspanel,'units','normalized',...
    'pos',[.1 .15 .35 .75],'box','on','cameraposition',[5 -5 4]);
title('1D Fit'); xlabel('Lcc'); ylabel('Mcc'); zlabel('Response'); grid on; hold on;
fitspanel.axes.twoDfit = axes('parent',surffig.fitspanel,'units','normalized',...
    'pos',[.6 .15 .35 .75],'box','on','cameraposition',[5 -5 4]);
title('2D Fit'); xlabel('Lcc'); ylabel('Mcc'); zlabel('Response'); grid on; hold on;

% Default subunit is 1
conpanel.subselect = 1;

% Set up some fitting specifics
fitspanel.options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
fitspanel.ub.twoD = [max(GLMP.nspikes)             300  300  300 10 max(GLMP.nspikes)  pi 5];
fitspanel.lb.twoD = [min(GLMP.nspikes) 1/max(GLMP.rho)    0 -300  1              .001 -pi 0];
fitspanel.ub.oneD = fitspanel.ub.twoD;
fitspanel.lb.oneD = fitspanel.lb.twoD;
fitspanel.ub.oneD(4) = 0;
fitspanel.lb.oneD(4) = 0;
fitspanel.surftype = 'conicsection_xy';
fitspanel.errortype = 'NegativeBinomial';
fitspanel.thresh = .05; % For comparing types of 1D fits
fitspanel.params.oneD = [];
fitspanel.params.twoD = [];
fitspanel.a(1,:) = [0 -1 1 0 0 0 0 0];
fitspanel.a(2,:) = [-1 0 0 0 0 1 0 0];
fitspanel.b = zeros(size(fitspanel.a,1),1);

%Save variables
set(surffig.conpanel,'UserData',conpanel);
set(surffig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',surffig);


end

function LoadPreviousFits(analyze)
global GLMP GLMSPopData

% Load figure variables
surffig = get(gcf,'UserData');
conpanel = get(surffig.conpanel,'UserData');
fitspanel = get(surffig.fitspanel,'UserData');

% Load previously saved fits from Pop Data
sub = conpanel.subselect;
dfsub = strcat(GLMP.datafile,num2str(sub));
datatypes = GLMSPopData(1,:);
popdfs = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
popsubs = [GLMSPopData{2:end,strcmp(datatypes,'Subunit')}]';
popdfsubs = strcat(popdfs,num2str(popsubs));
popidx = find(strcmp(dfsub,popdfsubs)) + 1;
Lcc = GLMP.subunit{sub}.Lcc;
Mcc = GLMP.subunit{sub}.Mcc;
nsp = GLMP.subunit{sub}.nspikes;

% Fill in parameter values
if ~isempty(popidx)
    
    surfinfo = GLMSPopData{popidx,strcmp(datatypes,'Surface Parameters')};
    
    % Statistics Panel
    set(conpanel.paramvals.normLL,'string',rndofferr(surfinfo.twoD.normLL-surfinfo.oneD.normLL,4));
    set(conpanel.paramvals.pval,'string',rndofferr(surfinfo.LLRTpval,4));
    fitspanel.LLRTpval = surfinfo.LLRTpval;
    fitspanel.fitcomps = surfinfo.fitcomps;
    
    % 1D surface
    parvals = surfinfo.oneD.parvals;
    fitspanel.params.oneD = parvals;
    fitspanel.LL.oneD = surfinfo.oneD.LL;
    fitspanel.normLL.oneDnormLL = surfinfo.oneD.normLL;
    fitspanel.Hessian.oneD = surfinfo.oneD.Hessian;
    fitspanel.rsq.oneD = surfinfo.oneD.rsq;
    params1 = round(parvals,2);
    set(conpanel.paramvals.oneD.A,'string',params1(1));
    set(conpanel.paramvals.oneD.sig1,'string',round(1/params1(2),3));
    set(conpanel.paramvals.oneD.sig2,'string',round(1/params1(3),3));
    set(conpanel.paramvals.oneD.exp,'string',params1(5));
    set(conpanel.paramvals.oneD.bl,'string',params1(6));
    set(conpanel.paramvals.oneD.rot,'string',round(params1(7)/pi*180,0));
    set(conpanel.paramvals.oneD.kappa,'string',params1(8));
    set(conpanel.paramvals.oneD.LL,'string',round(surfinfo.oneD.LL,1));
    set(conpanel.paramvals.oneD.std,'string',round(sqrt(1/surfinfo.oneD.Hessian(end-1,end-1))/pi*180,1));
    set(conpanel.paramvals.oneD.normLL,'string',rndofferr(surfinfo.oneD.normLL,2));
    
    % Display 1D surface
    x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
    [xx,yy] = meshgrid(x,x);
    surface = ComputeNakaRushtonJPW(surfinfo.oneD.parvals,[xx(:) yy(:)],surfinfo.oneD.surftype);
    surface = reshape(surface,size(xx));
    axes(fitspanel.axes.oneDfit); cla; hold on; grid on;
    fitspanel.surf.oneD = surfc(xx,yy,surface);
    set(fitspanel.surf.oneD(1),'edgecolor','none')
    alpha(.5);
    set(gca,'view',[0 90]);
    uniquestim = unique([Lcc Mcc],'rows');
    maxnsp = max(GLMP.subunit{sub}.meannspikes);
    for i = 1:size(uniquestim,1)
        L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
        mn = mean(nsp(L))/maxnsp*10;
        h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
        set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
    end

    % 2D surface
    parvals = surfinfo.twoD.parvals;
    fitspanel.params.twoD = parvals;
    fitspanel.LL.twoD = surfinfo.twoD.LL;
    fitspanel.normLL.twoDnormLL = surfinfo.twoD.normLL;
    fitspanel.Hessian.twoD = surfinfo.twoD.Hessian;
    fitspanel.rsq.twoD = surfinfo.twoD.rsq;
    params1 = round(parvals,2);
    set(conpanel.paramvals.twoD.A,'string',params1(1));
    set(conpanel.paramvals.twoD.sig1,'string',round(1/params1(2),3));
    set(conpanel.paramvals.twoD.sig2,'string',round(1/params1(3),3));
    set(conpanel.paramvals.twoD.orthsig,'string',round(1/params1(4),3));
    set(conpanel.paramvals.twoD.exp,'string',params1(5));
    set(conpanel.paramvals.twoD.bl,'string',params1(6));
    set(conpanel.paramvals.twoD.rot,'string',round(params1(7)/pi*180,0));
    set(conpanel.paramvals.twoD.kappa,'string',params1(8));
    set(conpanel.paramvals.twoD.LL,'string',round(surfinfo.twoD.LL,1));
    set(conpanel.paramvals.twoD.std,'string',round(sqrt(1/surfinfo.twoD.Hessian(end-1,end-1))/pi*180,1));
    set(conpanel.paramvals.twoD.normLL,'string',rndofferr(surfinfo.twoD.normLL,2));
    
    % Display 2D surface
    x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
    [xx,yy] = meshgrid(x,x);
    surface = ComputeNakaRushtonJPW(surfinfo.twoD.parvals,[xx(:) yy(:)],surfinfo.twoD.surftype);
    surface = reshape(surface,size(xx));
    axes(fitspanel.axes.twoDfit); cla; hold on;
    fitspanel.surf.twoD = surfc(xx,yy,surface);
    set(fitspanel.surf.twoD(1),'edgecolor','none')
    alpha(.5);
    set(gca,'view',[0 90]);
    uniquestim = unique([Lcc Mcc],'rows');
    maxnsp = max(GLMP.subunit{sub}.meannspikes);
    for i = 1:size(uniquestim,1)
        L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
        mn = mean(nsp(L))/maxnsp*10;
        h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
        set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
    end
    
    drawnow;
     
end

% Save figure variables
set(surffig.conpanel,'userdata',conpanel);
set(surffig.fitspanel,'userdata',fitspanel);
set(gcf,'userdata',surffig);

% Execute function if requested
if ~isempty(analyze)  
    startanalysis()
end

end
