% This script is a model of a V1 cell tuned to various directions in the LM
% plane.  We want to know how various sampling schemes may bias the
% results.
function VRFig3c()

% Set up figure
figure(1003); clf;
set(gcf,'units','normalized','pos',[.15 .15 .7 .7],'NumberTitle','off',...
    'Name','Nonlinearly Transformed RS Distribution');
modelfig.conpanel = uipanel('parent',gcf,'pos',[.01 .01 .19 .98],'title','Control Panel');
modelfig.fitspanel = uipanel('parent',gcf,'pos',[.21 .01 .78 .59]);
modelfig.surfpanel = uipanel('parent',gcf,'pos',[.21 .61 .78 .38]);

% Set up controls for analysis
conpanel.uicontrols.nsamps = uicontrol('style','edit','parent',modelfig.conpanel,...
    'units','normalized','pos',[.1 .8 .35 .08],'string',100,'fontsize',10);
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
conpanel.uicontrols.nAngs = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.55 .6 .35 .08],'string',33,'fontsize',10);
conpanel.labels.nAngs = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.55 .68 .35 .07],'style','text','HorizontalAlignment','center',...
    'string','Total Tested Angles','fontsize',10);
conpanel.uicontrols.maxcolcc = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.1 .4 .35 .08],'string',.09,'fontsize',10);
conpanel.labels.maxcolcc = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.1 .48 .35 .07],'style','text','HorizontalAlignment','center',...
    'string','Max Chromatic Contrast','fontsize',10);
conpanel.uicontrols.scale = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.55 .4 .35 .08],'string',4,'fontsize',10);
conpanel.labels.scale = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.55 .48 .35 .07],'style','text','HorizontalAlignment','center',...
    'string','Nonlinear Scale','fontsize',10);
conpanel.uicontrols.startanalysis = uicontrol('parent',modelfig.conpanel,'style','pushbutton',...
    'units','normalized','pos',[.2 .1 .6 .15],'string','Start Analysis','fontsize',10,...
    'backgroundColor',[.2 1 .2],'callback',@StartAnalysis);

% Set up surfpanel
surfpanel.axes(1) = axes('parent',modelfig.surfpanel,'units','normalized',...
    'pos',[.025 .1 .2 .8],'box','on','xtick',[],'Ytick',[]); hold on;
surfpanel.axes(2) = axes('parent',modelfig.surfpanel,'units','normalized',...
    'pos',[.275 .1 .2 .8],'box','on','xtick',[],'Ytick',[]); hold on;
surfpanel.axes(3) = axes('parent',modelfig.surfpanel,'units','normalized',...
    'pos',[.525 .1 .2 .8],'box','on','xtick',[],'Ytick',[]); hold on;
surfpanel.axes(4) = axes('parent',modelfig.surfpanel,'units','normalized',...
    'pos',[.775 .1 .2 .8],'box','on','xtick',[],'Ytick',[]); hold on;

% Set up fitspanel
fitspanel.axes.error = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.05 .15 .9 .7],'box','on','xlim',[-180 180]);
xlabel('True Angle'); ylabel('Error in Fit Angle (Deg)'); title('Error in Fit')

%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(gcf,'UserData',modelfig);


end

function StartAnalysis(~,~)
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
set(conpanel.uicontrols.startanalysis,'string','Modeling in Progress...','fontsize',10,...
    'backgroundColor',[1 .5 .2],'callback',[]);

% GL struct stores experimental variables and tracks the progress of the model.
gl.nPres = str2double(get(conpanel.uicontrols.npresstim,'string'));
gl.nSamps = str2double(get(conpanel.uicontrols.nsamps,'string'));
gl.nRnds = str2double(get(conpanel.uicontrols.nRnds,'string'));
gl.maxcolcc = str2double(get(conpanel.uicontrols.maxcolcc,'string'));
gl.nAngs = str2double(get(conpanel.uicontrols.nAngs,'string'));
gl.allAngs = linspace(-pi/4,3*pi/4,gl.nAngs);
%gl.allAngs = 0;

% Preallocate some space
gl.RWAfit = nan(numel(gl.allAngs),gl.nSamps);
gl.RWAcorrectfit = nan(numel(gl.allAngs),gl.nSamps);
gl.LLfit = nan(numel(gl.allAngs),gl.nSamps);

% Select all the stimuli to be shown. Same for all trials.
disp('Choosing Stim...')
ChooseLMStimuli;
Calculatec50s;

% Calculate RWA and LL for all samples
for rot = 1:numel(gl.allAngs)
    gl.currentAng = rot;    
    for sampn = 1:gl.nSamps
                      
        gl.currentSamp = sampn;
        disp(['Angle ' num2str(gl.currentAng) ' of ' num2str(numel(gl.allAngs))])
        disp(['Sample ' num2str(gl.currentSamp) ' of ' num2str(gl.nSamps)])
  
        % Generate responses for each sample dataset
        disp('Modeling Responses...')
        CreateModelSurface;
        
        % Calculate pref axis based on RWA
        disp('Estimating preferred direction with distended RWA...')
        weightresp = repmat(gl.nsp,[1 2]) .* [gl.Lcc gl.Mcc];
        RWA = sum(weightresp)./numel(gl.Lcc);
        [fitAng,~] = cart2pol(RWA(1),RWA(2));
        gl.RWAfit(gl.currentAng,gl.currentSamp) = fitAng - gl.allAngs(gl.currentAng);
        
        % Calculate RWA on whitened data
        disp('Estimating preferred direction with RS RWA...')
        whitemat = sqrtm(inv(cov([gl.Lcc gl.Mcc])));
        whitepts = [gl.Lcc gl.Mcc] * whitemat';
        WLcc = whitepts(:,1);
        WMcc = whitepts(:,2);
        weightresp = repmat(gl.nsp,[1 2]) .* [WLcc WMcc];
        RWA = sum(weightresp)./numel(WLcc);
        tRWA = RWA * whitemat;
        [fitAng,~] = cart2pol(tRWA(1),tRWA(2));
        gl.RWAcorrectfit(gl.currentAng,gl.currentSamp) = fitAng - gl.allAngs(gl.currentAng);
        
        % Calculate pref axis based on LL
        disp('Estimating preferred direction by LL...')
        fitAng = OneDFit();
        gl.LLfit(gl.currentAng,gl.currentSamp) = fitAng - gl.allAngs(gl.currentAng);
        
    end
end

% Calculate error mean and std
gl.RWAfitmean = mean(gl.RWAfit./pi*180,2);
gl.RWAfitstd = std(gl.RWAfit./pi*180,[],2);
gl.RWAcorrectfitmean = mean(gl.RWAcorrectfit./pi*180,2);
gl.RWAcorrectfitstd = std(gl.RWAcorrectfit./pi*180,[],2);
gl.LLfitmean = mean(gl.LLfit./pi*180,2);
gl.LLfitstd = std(gl.LLfit./pi*180,[],2);

% Fill in display axes
PlotFigs()
        
% Turn start button back on
set(conpanel.uicontrols.startanalysis,'string','Start Analysis','fontsize',10,...
    'backgroundColor',[.2 1 .2],'callback',@StartAnalysis);

% Save Data
data.angs = gl.allAngs;
data.RWAraw.means = gl.RWAfitmean;
data.RWAraw.std = gl.RWAfitstd;
data.RWAcorrect.means = gl.RWAcorrectfitmean;
data.RWAcorrect.std = gl.RWAcorrectfitstd;
data.LL.means = gl.LLfitmean;
data.LL.std = gl.LLfitstd; 
if ismac
    save('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig3cData(NLtransRS)','data');
else ispc
    save('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig3cData(NLtransRS)','data');
end

% Say goodbye
disp('Fig 3c Modeling Completed.')


end


function fitAng = OneDFit()
global gl

% Set up some variables
angs = linspace(-pi,pi,9);
GOF = nan(numel(angs),1);
vub = [max(gl.nsp)*1.1   50   50   5     max(gl.nsp)    pi];
vlb = [.01           .01  .01   .5      0        -pi];
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','lbfgs',...
    'DiffMinChange',.0001,'display','off');
Aguess = max(gl.nsp)*.9;
sigguess = [.1 50];
expguess = 2;
blguess = min(gl.nsp)+1;
params0 = nan(numel(angs),numel(vlb));
guessIdx = fullfact([numel(sigguess) numel(angs)]); 
A = [0 1 -1 0 0 0];
b = 0; 

% Fit using a variety of initial guesses
for rot = 1:size(guessIdx,1)
    paramsGuess = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,1)) expguess blguess angs(guessIdx(rot,2))];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,A,b,[],[],vlb,vub,[],options,[gl.Lcc gl.Mcc],gl.nsp,'surface7','Poisson');
    params0(rot,:) = f1;
    GOF(rot) = fval;
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params0(bestIdx,:);

% Save results
fitAng = params1(end);

end


function PlotFigs()
global gl

% Load variables
modelfig = get(gcf,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Plot select angles without the noise
whichangs = round(linspace(1,numel(gl.allAngs),6));
whichangs = whichangs(2:end-1);
for n = 1:numel(whichangs)
    
    % Retrieve stimuli and parameters
    ang = gl.allAngs(whichangs(n));
    
    % For setting sigma
    rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
    projpts = [gl.Lcc gl.Mcc] * rotMat;
    sig = max(projpts(:,2))/2;
    
    % Set up some variables
    gl.params(2) = sig;
    gl.params(end) = ang;
    scalefac = 20;
    uniquepts = unique([gl.Lcc gl.Mcc],'rows');

    % Responses to stimuli
    nsp = ComputeNakaRushtonJPW(gl.params,uniquepts,'surface7');
    
    %axes(surfpanel.axes(n)); cla; grid on; hold on; 
    for i = 1:size(uniquepts,1)
        h = plot(surfpanel.axes(n),uniquepts(i,1),uniquepts(i,2),'ko');
        set(h,'MarkerFaceColor','white','MarkerSize',(((nsp(i)-min(nsp))/max(nsp))+.2).*scalefac,'MarkerEdgeColor','black');
    end
    X = linspace(min(gl.Lcc),max(gl.Lcc),100);
    Y = linspace(min(gl.Mcc),max(gl.Mcc),100);
    [XX,YY] = meshgrid(X,Y);
    ZZ = ComputeNakaRushtonJPW(gl.params,[XX(:) YY(:)],'surface7');
    ZZ = reshape(ZZ,size(XX));
    contour(surfpanel.axes(n),XX,YY,ZZ);
    text(min(gl.Lcc)*.9,max(gl.Mcc)*.5,['ang = ' num2str(ang/pi*180)])
    
    % Plot standard deviation of angle
    angub = ang + gl.LLfitmean(whichangs(n))/180*pi + gl.LLfitstd(whichangs(n))/180*pi;
    anglb = ang + gl.LLfitmean(whichangs(n))/180*pi - gl.LLfitstd(whichangs(n))/180*pi;
    [x1,y1] = pol2cart(angub,max(gl.Lcc)*2);
    [x2,y2] = pol2cart(anglb,max(gl.Lcc)*2);
    patch([0 x1 x2],[0 y1 y2],[.5 .5 .8],'g','parent',surfpanel.axes(n));
    alpha(.4)
    xlim(surfpanel.axes(n),[min(gl.Lcc)*1.1 max(gl.Lcc)*1.1]);
    ylim(surfpanel.axes(n),[min(gl.Mcc)*1.1 max(gl.Mcc)*1.1]);

end

% Plot difference in prediction and varance
angs = gl.allAngs./pi*180;
axes(fitspanel.axes.error); cla; hold on;
h = shadedErrorBar(angs,gl.RWAfitmean',gl.RWAfitstd,'r-');
set(h.patch,'FaceAlpha',.4)
h = shadedErrorBar(angs,gl.LLfitmean',gl.LLfitstd,'g-');
set(h.patch,'FaceAlpha',.4)
h = shadedErrorBar(angs,gl.RWAcorrectfitmean',gl.RWAcorrectfitstd,'b-');
set(h.patch,'FaceAlpha',.4)
xlim([min(angs) max(angs)])
maxy = max(abs(cat(1,(gl.RWAfitstd+gl.RWAfitmean),(gl.LLfitstd+gl.LLfitmean),gl.RWAcorrectfitstd+gl.RWAcorrectfitmean)));
ylim([-maxy*1.1 maxy*1.1])

end

function Calculatec50s()
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');

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
rhothetalist = [thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))];
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
Lcc = tempLcc .* gl.maxcc;
Mcc = tempMcc .* gl.maxcc;

% Linear transform of RS Dataset
gl.linscale = str2double(get(conpanel.uicontrols.scale,'string'));
scalemat = [gl.linscale 0; 0 1] * [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
newpts = [Lcc Mcc] * scalemat;

for n = 1:numel(gl.allAngs)
    rotMat = [cos(gl.allAngs(n)) sin(gl.allAngs(n)); -sin(gl.allAngs(n)) cos(gl.allAngs(n))];
    projpts = newpts * rotMat;
    gl.sig(n) = max(projpts(:,2))/2;
end

end


function ChooseLMStimuli()
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
colmaxcc = str2double(get(conpanel.uicontrols.maxcolcc,'string'));
scale = str2double(get(conpanel.uicontrols.scale,'string'));

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

% Nonlinear transform of RS Dataset
lummaxcc = colmaxcc * scale;
scale = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(lummaxcc.*sin(thetas-pi/4)).^2);
newpts = [tempLcc tempMcc] .* repmat(scale,[1 2]);
gl.Lcc = newpts(:,1);
gl.Mcc = newpts(:,2);

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end


function CreateModelSurface()
global gl 

% Load Figure Variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');

% Retrieve stimuli and parameters
ang = gl.allAngs(gl.currentAng);

% For setting sigma
% colmaxcc = str2double(get(conpanel.uicontrols.maxcolcc,'string'));
% linscale = str2double(get(conpanel.uicontrols.scale,'string'));
% lummaxcc = colmaxcc * linscale;
% ellipse = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(ang-pi/4)).^2 ...
%     +(lummaxcc.*sin(ang-pi/4)).^2);
% sig = ellipse/2;
sig = gl.sig(gl.currentAng);

% Set up some variables
gl.params = [50 sig 100 3 .1 ang];

% Responses to stimuli
nsp = ComputeNakaRushtonJPW(gl.params,[gl.Lcc gl.Mcc],'surface7');
gl.nsp = poissrnd(nsp);

% Save Figure Variables
%set(gcf,'UserData',modelfig);

end

