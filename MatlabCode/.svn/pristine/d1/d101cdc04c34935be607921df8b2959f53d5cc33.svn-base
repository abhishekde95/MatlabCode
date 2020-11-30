% This script is a model of a V1 cell tuned to various directions in the LM
% plane.  We want to know how various sampling schemes may bias the
% results.
function MeasVarOnTuningModel()
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

% Set up controls for analysis
conpanel.uicontrols.nsamps = uicontrol('style','edit','parent',modelfig.conpanel,...
    'units','normalized','pos',[.1 .8 .35 .08],'string',3,'fontsize',10);
conpanel.labels.nsamps = uicontrol('Parent',modelfig.conpanel,'Units','normalized',...
    'pos',[.1 .88 .35 .1],'HorizontalAlignment','center',...
    'style','text','string','Datasets per "Real Tuning Direction"','FontSize',10);
conpanel.uicontrols.npresstim = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.55 .8 .35 .08],'string',5,'fontsize',10);
conpanel.labels.npresstim = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.55 .88 .35 .1],'style','text','HorizontalAlignment','center',...
    'string','Samples per (L,M) Value','fontsize',10);
conpanel.uicontrols.nRnds = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.1 .6 .35 .08],'string',3,'fontsize',10);
conpanel.labels.nRnds = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.1 .68 .35 .07],'style','text','HorizontalAlignment','center',...
    'string','Rounds per Dataset','fontsize',10);
conpanel.uicontrols.lumcc = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.1 .4 .35 .08],'string',.9,'fontsize',10);
conpanel.labels.lumcc = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.1 .48 .35 .07],'style','text','HorizontalAlignment','center',...
    'string','Max Lum Contrast','fontsize',10);
conpanel.uicontrols.colcc = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.6 .4 .35 .08],'string',.09,'fontsize',10);
conpanel.labels.colcc = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.6 .48 .35 .07],'style','text','HorizontalAlignment','center',...
    'string','Max Color Contrast','fontsize',10);

conpanel.uicontrols.startanalysis = uicontrol('parent',modelfig.conpanel,'style','pushbutton',...
    'units','normalized','pos',[.3 .1 .4 .15],'string','Start Analysis','fontsize',10,...
    'backgroundColor',[.2 1 .2],'callback',@StartAnalysis);

% Set up surfpanel
surfpanel.axes = axes('parent',modelfig.surfpanel,'units','normalized','pos',[.1 .1 .8 .8],...
    'box','on');

% Set up fitspanel
fitspanel.axes.LL = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .6 .8 .325],'box','on','xlim',[-pi/2 pi/2]);
xlabel('True Angle'); ylabel('Difference from Fit Angle'); title('Error in Fit')
fitspanel.axes.Hessian = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .1 .8 .325],'box','on','xlim',[-pi/2 pi/2]);
xlabel('True Angle'); ylabel('Confidence in Fit'); title('Confidence in Fit')

%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);
    

end

function StartAnalysis(~,~)
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
%surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

set(conpanel.uicontrols.startanalysis,'string','Analysis in Progress...','fontsize',10,...
    'backgroundColor',[1 .5 .2],'callback',[]);

% GL struct stores experimental variables and tracks the progress of the
% model.
gl.nPres = str2double(get(conpanel.uicontrols.npresstim,'string'));
gl.nSamps = str2double(get(conpanel.uicontrols.nsamps,'string'));
gl.nRnds = str2double(get(conpanel.uicontrols.nRnds,'string'));
gl.lummaxcc = str2double(get(conpanel.uicontrols.lumcc,'string'));
gl.colmaxcc = str2double(get(conpanel.uicontrols.colcc,'string'));
gl.allAngs = linspace(-pi/2,pi/2,17);
gl.allAngs = 0; % This is a hack to make a single direction

% Preallocate some space
fitData = nan(numel(gl.allAngs),gl.nSamps);
hessData = nan(numel(gl.allAngs),gl.nSamps);

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
        OneDFit
        
        % Load in figure variables
        fitspanel = get(modelfig.fitspanel,'UserData');
        angdiff = fitspanel.fitAng - gl.allAngs(gl.currentAng);
        if angdiff > pi/2
            angdiff = angdiff - pi;
        elseif angdiff < -pi/2
            angdiff = angdiff + pi;
        end
        fitData(gl.currentAng,gl.currentSamp) = angdiff;
        hessData(gl.currentAng,gl.currentSamp) = 1/fitspanel.hessian;
    end
end

% Plot difference in prediction and varance
angs = gl.allAngs./pi*180;
fitmean = mean(fitData./pi*180,2);
fitvar = var(fitData./pi*180,[],2);
axes(fitspanel.axes.LL); hold on;
xlim([min(angs) max(angs)])
errorbar(angs,fitmean,fitvar);

% Plot confidence in the fit
meanhess = mean(hessData./pi*180,2);
varhess = var(hessData./pi*180,[],2);
axes(fitspanel.axes.Hessian); cla; hold on;
errorbar(angs,meanhess,varhess);
xlim([min(angs) max(angs)])

set(conpanel.uicontrols.startanalysis,'string','Start Analysis','fontsize',10,...
    'backgroundColor',[.2 1 .2],'callback',@StartAnalysis);

disp('Finished with Analysis.')

%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);
    

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

% Construct Polar Grid
if mod(gl.nRnds,2) == 0
    thetaspace = thetaspace / 2^((gl.nRnds-2)/2);
    rhospace = rhospace / 2^(gl.nRnds/2);
elseif mod(gl.nRnds,2) == 1
    thetaspace = thetaspace / 2^((gl.nRnds-1)/2);
    rhospace = rhospace / 2^((gl.nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = gl.lummaxcc*gl.colmaxcc./sqrt((gl.colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(gl.lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;

surfpanel.Lcc = Lcc;
surfpanel.Mcc = Mcc;

% % making some adjustments for figure creation
% figure(1); clf; hold on;
% plot(Lcc,Mcc,'ko')
% xlabel('Lcc'); ylabel('Mcc'); title('Model Surface')
% axis equal square

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
surfpanel = get(modelfig.surfpanel,'UserData');

% Retrieve stimuli and parameters
ang = gl.allAngs(gl.currentAng);
L = surfpanel.Lcc;
M = surfpanel.Mcc;

% For setting sigma
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = [L M] * rotMat;
sig = max(projpts(:,1))/2;

% Set up some variables
params = [50 sig 100 3 .1 ang];

% Responses to stimuli
nsp = ComputeNakaRushtonJPW(params,[L M],'surface7');
nsp = poissrnd(nsp);
surfpanel.nsp = nsp;

% Plot figure
axes(surfpanel.axes); cla; hold on;
plot3(L,M,nsp,'ko')
xlabel('Lcc'); ylabel('Mcc'); title('Model Surface')
axis equal square

% Save Figure Variables
set(modelfig.surfpanel,'UserData',surfpanel);
set(gcf,'UserData',modelfig);

end


function OneDFit()

% Load Figure Variables
modelfig = get(gcf,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Load data parameters
Lcc = surfpanel.Lcc;
Mcc = surfpanel.Mcc;
nsp = surfpanel.nsp;

% Set up some variables
angs = linspace(-pi/2,pi/2,9);
GOF = nan(numel(angs),1);
vub = [max(nsp)   50   50   5     max(nsp)    pi];
vlb = [.01       .01  .01   .5      0        -pi];
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','lbfgs',...
    'DiffMinChange',.05,'display','off');
Aguess = max(nsp);
sigguess = [.03 .3 50];
expguess = 2;
blguess = min(nsp)+1;
params = nan(numel(angs),numel(vlb));
guessIdx = fullfact([numel(sigguess) numel(angs)]); 

% Fit using a variety of initial guesses
for rot = 1:size(guessIdx,1)
    paramsGuess = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,1)) expguess blguess angs(guessIdx(rot,2))];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,[],[],[],[],vlb,vub,[],options,[Lcc Mcc],nsp,'surface7','Poisson');
    params(rot,:) = f1;
    GOF(rot) = fval;
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);
greg = @(a) FitNakaRushtonFunJPW(a,[Lcc Mcc],nsp,'surface7','Poisson');
[hess,~] = hessian(greg,params1);
if hess(end,end) < 0
    disp('Error with fit. Hess has a negative value...')
    keyboard
end

% Display surface and pts
x = linspace(-max(Lcc),max(Lcc),50);
y = linspace(-max(Mcc),max(Mcc),50);
[xx,yy] = meshgrid(x,y);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
surface = reshape(surface,size(xx));
axes(surfpanel.axes); cla; hold on;
h = surfc(xx,yy,surface);
set(h,'edgecolor','none')
alpha(.5);
plot3(Lcc,Mcc,nsp,'k*');
set(gca,'view',[0 90]);

% Save results
fitspanel.fitAng = params1(end);
fitspanel.hessian = hess(end,end);

% Save Figure Variables
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end



