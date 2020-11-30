% for building new surfaces

x = linspace(-1,1,51);
[X,Y] = meshgrid(x,x);

% Assign parameter values
A = 50;
bl = 0;
exp = 2;
sig1 = .5;
sig2 = 500;
sig3 = -.5;
sig4 = 1000;
rot = 0;

%% Surface 9
% 1D Surface. Allows for suppression.  Negative sigma indicates suppression negative. 
% Upper asymptote is bl, lower asymptote is 0.

% Assign inputs
params = [A sig1 sig3 exp bl rot];
surftype = 'surface9';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(1); clf; hold on; grid on; box on;
surface(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)


%% Surface 9a
% 2D version of Surface9. Allows for suppression.  Negative sigma indicates suppression negative. 
% Upper asymptote is bl, lower asymptote is 0. Both sigma and A change as an ellipse.

% Assign inputs
params = [A sig1 sig2 sig3 sig4 exp bl rot];
surftype = 'surface9a';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(2); clf; hold on; grid on; box on;
surface(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)

%% Surface 9b
% 2D version of Surface9. Allows for suppression.  Negative sigma indicates suppression negative. 
% Upper asymptote is bl, lower asymptote is 0. 
% Orthogonal 1D mechanisms are summed.

% Assign inputs
params = [A sig1 sig2 sig3 sig4 exp bl rot];
surftype = 'surface9b';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(3); clf; hold on; grid on; box on;
surface(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)

%% Surface 9c
% 2D version of Surface9. Allows for suppression.  Negative sigma indicates suppression negative. 
% Upper asymptote is bl, lower asymptote is 0. 
% Orthogonal 1D mechanisms are averaged.

% Assign inputs
params = [A sig1 sig2 sig3 sig4 exp bl rot];
surftype = 'surface9c';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(4); clf; hold on; grid on; box on;
surface(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)

%% Surface 9d
% 2D version of Surface9. Allows for suppression.  Negative sigma indicates suppression negative. 
% Upper asymptote is bl, lower asymptote is 0. 
% Orthogonal 1D mechanisms are maxed.

% Assign inputs
params = [A sig1 sig2 sig3 sig4 exp bl rot];
surftype = 'surface9d';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(5); clf; hold on; grid on; box on;
surface(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)

%% Surface 10
% 2D Surface. Allows for hyper tuning (suppression to orthogonal stimuli)

% Assign inputs
s = 2;
params = [A sig1 sig2 exp s bl rot];
surftype = 'surface10';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(6); clf; hold on; grid on; box on;
surface(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)

%% Surface 11
% Allows for boardly tuned fits

% Assign parameter values
s = 2;
params = [A sig1 sig2 exp s bl rot];
surftype = 'surface11';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(7); clf; hold on; grid on; box on;
%surface(X,Y,Z)
surfc(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)

%% Surface 12
% sandbox with greg

% Assign parameter values
orthocoef = 5;
anticoef = 0;
exp = 2;
params = [A sig1 orthocoef anticoef exp bl rot];
surftype = 'surface12';
Z = ComputeNakaRushtonJPW(params,[X(:) Y(:)],surftype);
Z = reshape(Z,size(X));

% Plot resulting surface
figure(7); clf; hold on; grid on; box on;
%surface(X,Y,Z)
surfc(X,Y,Z)
alpha(.5)
zlim([0 max(Z(:))])
set(gcf,'Name',surftype)


%% Conic Section
% Allows for sigmas that change as a conic section

% generate grid
x = linspace(-1,1,51);
[X,Y] = meshgrid(x,x);

% Generate fake neural data
A = 40;
sig1 = .5;
sig2 = Inf;
orthosig = .5;
exp = 3;
bl = 0;
rot = 0;
params0 = [A 1/sig1 1/sig2 1/orthosig exp bl rot];
surftype0 = 'conicsection';

figure(10); cla; hold on; grid on; box on;
% stim = 2*rand(500,2)-1;
% resp = poissrnd(ComputeNakaRushtonJPW(params0,stim,surftype0));
% plot3(stim(:,1),stim(:,2),resp,'ko');
stim = [X(:) Y(:)];
resp = ComputeNakaRushtonJPW(params0,stim,surftype0);
surf(X,Y,reshape(resp,size(X)));
xlim([min(x) max(x)])
ylim([min(x) max(x)])

%% Load in real datasets to test new surface
global GLMSPopData

% Load data from pop structure
popidx = 100;
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{popidx+1,strcmp(datatypes,'GLMP')};
DN = GLMSPopData{popidx+1,strcmp(datatypes,'DN')};
p = GLMSPopData{popidx+1,strcmp(datatypes,'Subunit')};

% Load data parameters
Lcc = cat(1,GLMP.subunit{p}.Lcc,zeros(size(GLMP.subunit{p}.blnspikes)));
Mcc = cat(1,GLMP.subunit{p}.Mcc,zeros(size(GLMP.subunit{p}.blnspikes)));
nsp = cat(1,GLMP.subunit{p}.nspikes,GLMP.subunit{p}.blnspikes);

% Fitting parameters
surftype = 'conicsection';
errortype = 'NegativeBinomial';
vub = [max(GLMP.subunit{p}.nspikes) 100  100 100 4 50  pi 2];
vlb = [min(GLMP.subunit{p}.nspikes)   1 .001 -10 2  0 -pi 0];
vub0 = vub;
vlb0 = vlb;
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs',...
    'display','off','TolFun',10.^-9);

% Set initial guesses and options for 1D fit
angs = unique(GLMP.subunit{p}.theta);
Aguess = max(GLMP.subunit{p}.nspikes) * .8;
orsigguess = 0;
expguess = 3;
blguess = mean(GLMP.subunit{p}.blnspikes);
kappaguess = regress(GLMP.subunit{p}.varnspikes-GLMP.subunit{p}.meannspikes,...
    GLMP.subunit{p}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < 0
    kappaguess = 0;
end
%A = zeros(size(vub));
A = [0 -1 1 0 0 0 0 0];
b = 0;
bidnLL = Inf;

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
for rot = 1:size(angs)
    
    % Indices (+ and -) for a given axis
    L = GLMP.subunit{p}.theta == angs(rot)...
        | GLMP.subunit{p}.theta == angs(rot)-pi...
        | GLMP.subunit{p}.theta == angs(rot)+pi;
    
    % Grab contrasts and responses for each data axis
    Lcc0 = cat(1,GLMP.subunit{p}.Lcc(L),zeros(size(GLMP.subunit{p}.blnspikes)));
    Mcc0 = cat(1,GLMP.subunit{p}.Mcc(L),zeros(size(GLMP.subunit{p}.blnspikes)));
    nsp0 = cat(1,GLMP.subunit{p}.nspikes(L),GLMP.subunit{p}.blnspikes);

    % Fill in parameter guesses
    sigguess = mean(GLMP.subunit{p}.rho(L));
    paramsGuess0 = [Aguess sigguess sigguess orsigguess expguess blguess angs(rot) kappaguess];

    % Fit function to single spoke as a guess for larger function
    vub0(end-1) = angs(rot);
    vlb0(end-1) = angs(rot);
    [fitparams,~] = fmincon('FitNakaRushtonFunJPW',paramsGuess0,...
        A,b,[],[],vlb0,vub0,[],options,[Lcc0 Mcc0],nsp0,...
        surftype,errortype);
    
    % Fit all of the data using axis fit as the initial guess
    [f1,fval,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',fitparams,...
        A,b,[],[],vlb,vub,[],options,[Lcc Mcc],nsp,...
        surftype,errortype);
    
    % If better than previous best fit, replace parmeter values
    if fval < bidnLL
        params = f1;
        nLL = fval;
        hessval = hess;
    end
end

% Display surface and pts
x = linspace(-max(GLMP.subunit{p}.rho),max(GLMP.subunit{p}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],surftype);
surface = reshape(surface,size(xx));
figure(1); clf; hold on; grid on;
surf0 = surfc(xx,yy,surface);
set(surf0,'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);
pts0 = plot3(Lcc,Mcc,nsp,'k*');

disparams = params;
disparams(end-1) = params(end-1)/pi*180;
disparams(2:4) = 1./params(2:4)

    

