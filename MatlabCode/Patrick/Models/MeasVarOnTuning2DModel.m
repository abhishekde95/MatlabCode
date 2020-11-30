% This script models the responses of an LNLN neuron tuned to various 
% directions in the LM plane.  We want to know how the accuracy changes 
% with preferred direction.

function MeasVarOnTuning2DModel()
close all

% Set up figure
figure(1000); clf;
set(gcf,'units','pixels','pos',[50 100 700 700],'NumberTitle','off',...
    'Name','Variance of Tuning Model');
modelfig = get(gcf,'UserData');
modelfig.conpanel = uipanel('pos',[.525 .525 .45 .45],'parent',gcf,'title','Control Panel');
modelfig.surfpanel = uipanel('Pos',[.025 .525 .45 .45],'Parent',gcf);
modelfig.fitspanel = uipanel('pos',[.025 .025 .95 .45],'parent',gcf);
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% User-defined variables
surfpanel.realparams = [20 nan nan nan 2 .1 nan .5];
surfpanel.surftype = 'conicsection_xy';
surfpanel.errortype = 'NegativeBinomial';

% Set up controls for analysis
conpanel.uicontrols.nsamps = uicontrol('style','edit','parent',modelfig.conpanel,...
    'units','normalized','pos',[.4 .8 .2 .07],'string',10,'fontsize',10);
conpanel.labels.nsamps = uicontrol('Parent',modelfig.conpanel,'Units','normalized',...
    'pos',[.2 .88 .6 .1],'HorizontalAlignment','center',...
    'style','text','string','Datasets per "Real Tuning Direction"','FontSize',10);
conpanel.uicontrols.npresstim = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.4 .55 .2 .07],'string',5,'fontsize',10);
conpanel.labels.npresstim = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.2 .63 .6 .1],'style','text','HorizontalAlignment','center',...
    'string','Samples per (L,M) Value','fontsize',10);
conpanel.uicontrols.nRnds = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.4 .3 .2 .07],'string',3,'fontsize',10);
conpanel.labels.nRnds = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.2 .38 .6 .1],'style','text','HorizontalAlignment','center',...
    'string','Rounds per Dataset','fontsize',10);
conpanel.uicontrols.startanalysis = uicontrol('parent',modelfig.conpanel,'style','pushbutton',...
    'units','normalized','pos',[.3 .1 .4 .15],'string','Start Analysis','fontsize',10,...
    'backgroundColor',[.2 1 .2],'callback',@StartAnalysis);
if ismac
    conpanel.library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    conpanel.library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Set up surfpanel
surfpanel.axes = axes('parent',modelfig.surfpanel,'units','normalized','pos',[.1 .1 .8 .8]);

% Set up fitspanel
fitspanel.axes.LL = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .1 .8 .8],'box','on','xlim',[-pi/2 pi/2]);

%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

%load([conpanel.library 'VTM_conicsect_nbinvar'])
try
    % Organize and express in degrees
    angs = data.angs./pi*180;
    angdiffs = data.angdiff./pi*180;
    fitmean = mean(angdiffs,2);
    fitstd = std(angdiffs,[],2);
    fitvar = var(angdiffs,[],2);
    fitstderr = stderr(angdiffs')';
    
    % Plot
    axes(fitspanel.axes.LL); cla; hold on; grid on;
    plot(angs,angdiffs,'ko')
    shadedErrorBar(angs,fitmean,fitstd,'r-')
    xlim([min(angs) max(angs)])
    
catch
    return
end

end

function StartAnalysis(~,~)
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
%surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% GL struct stores experimental variables and tracks the progress of the model.
gl.nPres = str2double(get(conpanel.uicontrols.npresstim,'string'));
gl.nSamps = str2double(get(conpanel.uicontrols.nsamps,'string'));
gl.nRnds = str2double(get(conpanel.uicontrols.nRnds,'string'));
gl.allAngs = linspace(-pi/2,pi/2,73);
%gl.allAngs = linspace(-pi/2,pi/2,30);
data.angs = gl.allAngs;


% Rotate through each angle and number of samples
for rot = 1:numel(gl.allAngs)
    gl.currentAng = rot;
    for sampn = 1:gl.nSamps
        
        gl.currentSamp = sampn;
        
        disp(['Angle ' num2str(gl.currentAng) ' of ' num2str(numel(gl.allAngs))])
        disp(['Sample ' num2str(gl.currentSamp) ' of ' num2str(gl.nSamps)])
        
        disp('Choosing Stim...')
        ChooseLMStimuli;
        
        disp('Modeling Responses...')
        CreateModelSurface;
        
        disp('Fitting Model Data...')
        TwoDFit
        
        % Load in figure variables
        fitspanel = get(modelfig.fitspanel,'UserData');
        angdiff = fitspanel.params(end-1) - gl.allAngs(gl.currentAng);
        if angdiff > pi/2
            angdiff = angdiff - pi;
        elseif angdiff < -pi/2
            angdiff = angdiff + pi;
        end
        data.angdiff(gl.currentAng,gl.currentSamp) = angdiff;
        data.params{gl.currentAng,gl.currentSamp} = fitspanel.params;
        data.LL(gl.currentAng,gl.currentSamp) = fitspanel.LL;
                
%         disp('Saving data...')
%         if ismac
%             library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
%         elseif ispc
%             library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
%         end
%         save([library 'VTM_conicsect_nbinvar'],'data')
        
    end
end

disp('Saving data...')
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
save([library 'VTM_PC'],'data')


% Organize and express in degrees
angs = data.angs./pi*180;
angdiffs = data.angdiff./pi*180;
fitmean = mean(angdiffs,2);
fitvar = var(angdiffs,[],2);
fitstderr = stderr(angdiffs')';

% Plot
axes(fitspanel.axes.LL); cla; hold on; grid on;
%shadedErrorBar(angs,fitstderr,'r-')
plot(angs,fitmean,'r--')
xlim([min(angs) max(angs)])
plot(angs,angdiffs,'ko')

disp('Analysis finished.')

end

function ChooseLMStimuli()
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
colmaxcc = .09;
lummaxcc = .7;

% Construct Polar Grid
if gl.nRnds == 2
    rhospace = rhospace/2;
elseif gl.nRnds == 3
    rhospace = rhospace/2;
    thetaspace = thetaspace/2;
else
    disp('Hardcoded rounds 2 and 3, but not others yet...')
    keyboard
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;
surfpanel.Lcc = Lcc;
surfpanel.Mcc = Mcc;

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function CreateModelSurface()
global gl 

% Load Figure Variables
modelfig = get(gcf,'UserData');
%conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
%fitspanel = get(modelfig.fitspanel,'UserData');

% Plug in new "real" PD
params = surfpanel.realparams;
params(end-1) = gl.allAngs(gl.currentAng);
[pdx,pdy] = pol2cart(params(end-1),1);
proj = [pdx pdy] * [surfpanel.Lcc surfpanel.Mcc]';
params(2:3) = 1/(max(proj)/2);
[pdx,pdy] = pol2cart(params(end-1)-pi/2,1);
proj = [pdx pdy] * [surfpanel.Lcc surfpanel.Mcc]';
params(4) = 1/(max(proj)/2);
nsp = ComputeNakaRushtonJPW(params(1:7),[surfpanel.Lcc surfpanel.Mcc],surfpanel.surftype);

% Use as draw from negative binomial distribution
kappa = surfpanel.realparams(end);
mu = nsp;
sigsq = mu + kappa * mu.^2;
p = (sigsq - mu) ./ sigsq;
r = mu.^2 ./ (sigsq - mu);
surfpanel.nsp = nbinrnd(r,1-p);

% figure(1); clf; grid on; hold on;
% plot3(surfpanel.Lcc,surfpanel.Mcc,nsp,'r*')
% plot3(surfpanel.Lcc,surfpanel.Mcc,nbinrnd(r,1-p),'ko')

% Save Figure Variables
%set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
%set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function TwoDFit()

% Load Figure Variables
modelfig = get(gcf,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Load data parameters
Lcc = surfpanel.Lcc;
Mcc = surfpanel.Mcc;
nsp = surfpanel.nsp;
angs = linspace(-pi/2,pi/2,9);
GOF = nan(numel(angs),1);

% Set up some variables
[~,rho] = cart2pol(Lcc,Mcc);
ub = [max(nsp)        300 300 300 10 max(nsp)  pi 5];
lb = [min(nsp) 1/max(rho)   0   0  0       0  -pi 0];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
Aguess = max(nsp)*.8;
sigguess = max(rho)/2;
expguess = 2;
blguess = min(nsp);
kappaguess = 1;
params = nan(numel(angs),numel(lb));
guessIdx = fullfact([numel(sigguess) numel(angs)]); 

% Fit using a variety of initial guesses
for rot = 1:size(guessIdx,1)
    paramsGuess = [Aguess 1/sigguess(guessIdx(rot,1)) 1/sigguess(guessIdx(rot,1)) 1/sigguess(guessIdx(rot,1))...
        expguess blguess angs(guessIdx(rot,2)) kappaguess];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,a,b,[],[],lb,ub,[],options,[Lcc Mcc],nsp,surfpanel.surftype,surfpanel.errortype);
    params(rot,:) = f1;
    GOF(rot) = fval;
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Display surface and pts
x = linspace(-max(Lcc),max(Lcc),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],surfpanel.surftype);
surface = reshape(surface,size(xx));
axes(surfpanel.axes); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(nsp);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
h = surfc(xx,yy,surface);
set(h(1),'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);
colormap cool
%draw now

% Save results
fitspanel.params = params1;
fitspanel.LL = -GOF(bestIdx);

% Save Figure Variables
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

