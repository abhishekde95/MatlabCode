% This script creates fake neural responses and compares the fits of 1D
% and 2D models.
function  Likelihood_Ratio_Test()
global modeldata
%%
close all
clear all

% Create figure
figure(200); clf;
set(gcf,'position',[150 150 1200 800],....
    'numbertitle','off','name','Model of 1D vs 2D Fits')
ModelFig = get(gcf,'UserData');

%Set up model data
modeldata.params = [];

% Set up panels
ModelFig.ActualPanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.01 .31 .32 .68],'title','Modeled Neural Responses');
ModelFig.OneDPanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.34 .31 .32 .68],'title','1D Fit');
ModelFig.TwoDPanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.67 .31 .32 .68],'title','2D Fit');
ModelFig.ControlPanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.01 .01 .98 .29]);

% Set up axes
actualpanel.axes = axes('parent',ModelFig.ActualPanel,'units','normalized',...
    'position',[.15 .4 .75 .55],'cameraposition',[-0.7 -1.5 350],...
    'xlim',[-.1 .1],'ylim',[-.1 .1]);
xlabel('Lcc'); ylabel('Mcc'); zlabel('Response (sp/s)'); 
title('Modeled Responses'); grid on;

oneDpanel.axes = axes('parent',ModelFig.OneDPanel,'units','normalized',...
    'position',[.15 .4 .75 .55],'cameraposition',[-0.7 -1.5 350],...
    'xlim',[-.1 .1],'ylim',[-.1 .1]);
xlabel('Lcc'); ylabel('Mcc'); zlabel('Fit to Responses (sp/s)'); 
title('1D Model Fit'); grid on;

twoDpanel.axes = axes('parent',ModelFig.TwoDPanel,'units','normalized',...
    'position',[.15 .4 .75 .55],'cameraposition',[-0.7 -1.5 350],...
    'xlim',[-.1 .1],'ylim',[-.1 .1]);
xlabel('Lcc'); ylabel('Mcc'); zlabel('Fit to Responses (sp/s)'); 
title('2D Model Fit'); grid on;

% Set up parameter panels
actualpanel.parampanel = uipanel('parent',ModelFig.ActualPanel,'units','normalized',...
    'position',[.01 .01 .98 .3],'title','Model Parameters');
oneDpanel.parampanel = uipanel('parent',ModelFig.OneDPanel,'units','normalized',...
    'position',[.01 .01 .98 .3],'title','1D Fit Parameters');
twoDpanel.parampanel = uipanel('parent',ModelFig.TwoDPanel,'units','normalized',...
    'position',[.01 .01 .98 .3],'title','2D Fit Parameters');


%%% Actual Panel %%% 
% Set up inputs
actualpanel.parnames.a = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .8 .19 .15],'string','Max FR','fontsize',13);
actualpanel.parnames.sig = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .6 .19 .15],'string','Sigma','fontsize',13);
actualpanel.parnames.exp = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .4 .19 .15],'string','Exponent','fontsize',13);
actualpanel.parnames.bl = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .2 .19 .15],'string','Baseline','fontsize',13);
actualpanel.parnames.ang = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 0 .19 .15],'string','Angle (deg)','fontsize',13);

actualpanel.parvals.A = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.51 .81 .18 .18],'string',100);

actualpanel.parvals.s1 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.21 .61 .18 .18],'string',.045);
actualpanel.parvals.s2 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.41 .61 .18 .18],'string',.045);
actualpanel.parvals.s3 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.61 .61 .18 .18],'string',.045);
actualpanel.parvals.s4 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.81 .61 .18 .18],'string',.045);

actualpanel.parvals.e1 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.21 .41 .18 .18],'string',2);
actualpanel.parvals.e2 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.41 .41 .18 .18],'string',2);
actualpanel.parvals.e3 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.61 .41 .18 .18],'string',2);
actualpanel.parvals.e4 = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.81 .41 .18 .18],'string',2);

actualpanel.parvals.bl = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.51 .21 .18 .18],'string',2);

actualpanel.parvals.ang = uicontrol('parent',actualpanel.parampanel,'units','normalized',...
    'style','edit','position',[.51 .01 .18 .18],'string',45);


%%% 1D and 2D Panels %%%

% Set up display of fit parameters
oneDpanel.parnames.a = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .8 .19 .15],'string','Max FR','fontsize',13);
twoDpanel.parnames.a = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .8 .19 .15],'string','Max FR','fontsize',13);
oneDpanel.parnames.sig = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .6 .19 .15],'string','Sigma','fontsize',13);
twoDpanel.parnames.sig = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .6 .19 .15],'string','Sigma','fontsize',13);
oneDpanel.parnames.exp = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .4 .19 .15],'string','Exponent','fontsize',13);
twoDpanel.parnames.exp = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .4 .19 .15],'string','Exponent','fontsize',13);
oneDpanel.parnames.bl = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .2 .19 .15],'string','Baseline','fontsize',13);
twoDpanel.parnames.bl = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 .2 .19 .15],'string','Baseline','fontsize',13);
oneDpanel.parnames.ang = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 0 .19 .15],'string','Angle (deg)','fontsize',13);
twoDpanel.parnames.ang = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.01 0 .19 .15],'string','Angle (deg)','fontsize',13);

% For displaying 1D parameter values
oneDpanel.pardisp.A = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .81 .18 .14],'fontsize',10);
oneDpanel.pardisp.s = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .61 .3 .14],'fontsize',10);
oneDpanel.pardisp.e = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .41 .3 .14],'fontsize',10);
oneDpanel.pardisp.bl = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .21 .18 .14],'fontsize',10);
oneDpanel.pardisp.ang = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .01 .18 .14],'fontsize',10);

% For displaying 2D parameter values
twoDpanel.pardisp.A = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .8 .18 .14],'fontsize',10);
twoDpanel.pardisp.s = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .6 .3 .2],'fontsize',10);
twoDpanel.pardisp.e = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .4 .3 .2],'fontsize',10);
twoDpanel.pardisp.bl = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 .2 .18 .14],'fontsize',10);
twoDpanel.pardisp.ang = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.21 0 .18 .14],'fontsize',10);

% For displaying fit stats
oneDpanel.parnames.LL = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.61 .71 .18 .15],'string','LL = ','fontsize',13);
oneDpanel.parnames.df = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.61 .51 .18 .15],'string','DF = ','fontsize',13);
twoDpanel.parnames.LL = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.61 .71 .18 .15],'string','LL = ','fontsize',13);
twoDpanel.parnames.df = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.61 .51 .18 .15],'string','DF = ','fontsize',13);
twoDpanel.parnames.pval = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.61 .31 .18 .15],'string','p-val = ','fontsize',13);

oneDpanel.pardisp.LL = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.81 .71 .18 .14],'fontsize',13);
twoDpanel.pardisp.LL = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.81 .71 .18 .14],'fontsize',13);
oneDpanel.pardisp.df = uicontrol('parent',oneDpanel.parampanel,'units','normalized',...
    'style','text','position',[.81 .51 .18 .14],'fontsize',13);
twoDpanel.pardisp.df = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.81 .51 .18 .14],'fontsize',13);
twoDpanel.pardisp.pval = uicontrol('parent',twoDpanel.parampanel,'units','normalized',...
    'style','text','position',[.81 .31 .18 .14],'fontsize',13);


%%% Set up control panel %%%

conpanel.gopanel = uipanel('parent',ModelFig.ControlPanel,'units','normalized',...
    'position',[.01 .01 .32 .98],'title','Model Controls');
conpanel.run = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'style','pushbutton','string','Run','position',[.75 .1 .2 .25],...
    'backgroundcolor',[.2 1 .2],'Callback',@RunAnalysis);
conpanel.nLoops = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'style','edit','position',[.45 .825 .2 .15],'string',1);
conpanel.nRnds = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'style','edit','position',[.45 .625 .2 .15],'string',3);
conpanel.nPres = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'style','edit','position',[.45 .425 .2 .15],'string',5);
conpanel.maxLum = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'style','edit','position',[.45 .225 .2 .15],'string',.09);
conpanel.maxCol = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'style','edit','position',[.45 .025 .2 .15],'string',.09);

% Display control names
conpanel.parnames.nLoops = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'position',[.01 .8 .39 .15],'style','text','fontsize',11,'string','# of Loops');
conpanel.parnames.nRnds = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'position',[.01 .6 .39 .15],'style','text','fontsize',11,'string','# of Rounds');
conpanel.parnames.nPres = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'position',[.01 .4 .39 .15],'style','text','fontsize',11,'string','# of Presentations');
conpanel.parnames.maxCol = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'position',[.01 .2 .39 .15],'style','text','fontsize',11,'string','Max Luminance CC');
conpanel.parnames.maxLum = uicontrol('parent',conpanel.gopanel,'units','normalized',...
    'position',[.01 0 .39 .15],'style','text','fontsize',11,'string','Max Chromatic CC');

% Set up axes for p-vals
conpanel.histax = axes('parent',ModelFig.ControlPanel,'position',[.375 .2 .275 .7],...
    'fontsize',12,'xlim',[0 1]);
xlabel('p-values'); box on; grid on; hold on;

% Save User Variables
set(ModelFig.ActualPanel,'UserData',actualpanel);
set(ModelFig.OneDPanel,'UserData',oneDpanel);
set(ModelFig.TwoDPanel,'UserData',twoDpanel);
set(ModelFig.ControlPanel,'UserData',conpanel);
set(gcf,'UserData',ModelFig)


end

function RunAnalysis(~,~)
global modeldata


% Get User Data
ModelFig = get(gcf,'UserData');
actualpanel = get(ModelFig.ActualPanel,'UserData');
oneDpanel = get(ModelFig.OneDPanel,'UserData');
twoDpanel = get(ModelFig.TwoDPanel,'UserData');
conpanel = get(ModelFig.ControlPanel,'UserData');

% Get User-Defined variables
nLoops = str2double(get(conpanel.nLoops,'string'));
nRnds = str2double(get(conpanel.nRnds,'string'));
nPres = str2double(get(conpanel.nPres,'string'));
lumcc = .09;
colcc = .09;

% Get parameter values
A = str2double(get(actualpanel.parvals.A,'string'));
s1 = str2double(get(actualpanel.parvals.s1,'string'));
s2 = str2double(get(actualpanel.parvals.s2,'string'));
s3 = str2double(get(actualpanel.parvals.s3,'string'));
s4 = str2double(get(actualpanel.parvals.s4,'string'));
e1 = str2double(get(actualpanel.parvals.e1,'string'));
e2 = str2double(get(actualpanel.parvals.e2,'string'));
e3 = str2double(get(actualpanel.parvals.e3,'string'));
e4 = str2double(get(actualpanel.parvals.e4,'string'));
bl = str2double(get(actualpanel.parvals.bl,'string'));
ang = (str2double(get(actualpanel.parvals.ang,'string')))/180*pi; %Deg 2 rad

params = [A s1 s2 s3 s4 e1 e2 e3 e4 bl ang];

if numel(params) < 11
    disp('Not enough parameters. Set unused parameters to zero.')
    return
end

modeldata.params = params;

disp('Generating model responses...')

% Build a normalized polar grid
thetaspace = pi/4;
rhospace = .5;
thetas = shiftdim(0:(thetaspace/nRnds):(2*pi-(thetaspace/nRnds)));
rhos = shiftdim(rhospace/nRnds:rhospace/nRnds:1);

% Enumerate conditions
FFIdx = fullfact([numel(thetas) numel(rhos)]);
LMlist = repmat([thetas(FFIdx(:,1)) rhos(FFIdx(:,2))],nPres,1);

% Transform from polar to cartesian coordinates
[Lcc,Mcc] = pol2cart(LMlist(:,1),LMlist(:,2));
    
% Scale Cone Contrast Units
scale = lumcc*colcc./sqrt((colcc.*cos(LMlist(:,1)-pi/4)).^2 ...
    +(lumcc.*sin(LMlist(:,1)-pi/4)).^2);
Lcc = Lcc .* scale;
Mcc = Mcc .* scale;

% Pass L/M values into model
resps = ComputeNakaRushtonJPW(params,[Lcc(:) Mcc(:)],'surface5');

% Set up some parameters that only need one declaration
loopdata = cell(nLoops,1);
options = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
angs = linspace(0,pi,9);
angs(end) = [];
nrots = numel(angs);

%%% Start of Loops %%%
% Each loop is a new draw from the same underlying distribution from above
% (resps).  Each new draw is fit with a 1D and 2D function, and a
% log-liklihood ratio is used to determine which surface is more
% approriate. The p-value is plotted for each loop.

for loop = 1:nLoops
    
    disp(['Fitting sample ' num2str(loop) ' of ' num2str(nLoops) '...'])
    
    % Generate new model responses
    Presps = poissrnd(resps);
    loopdata{loop}.Lcc = Lcc;
    loopdata{loop}.Mcc = Mcc;
    loopdata{loop}.resps = Presps;
    
    % Fit with a surface for clarity
    x = linspace(min(Lcc),max(Lcc),50);
    [xx,yy] = meshgrid(x,x);
    zz = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],'surface5');
    zz = reshape(zz,size(xx));
    
    % Plot responses
    axes(actualpanel.axes); cla; hold on;
    actualpanel.surf = surfc(xx,yy,zz);
    alpha(.3);
    actualpanel.pts = plot3(Lcc,Mcc,Presps,'k*');
    
    %%% Fit 1D Model %%%
    disp('Fitting a 1D model...')
    
    % Set up some variables
    GOF = nan(nrots,1);
    vlb = [0    .001 .001 .001  .001    0  -pi];
    vub = [1000  50   50   10    10    100  pi];
    params0 = nan(nrots,numel(vlb));
    
    % Generate guesses for params
    Aguess = max(Presps);
    sigguess = mean([0 max(Lcc)]);
    expguess = 2;
    blguess = min(Presps);
    
    % Rotate through angles
    for rot = 1:nrots
        
        % Generating an initial guess
        paramsGuess = [Aguess sigguess sigguess expguess expguess blguess angs(rot)];

        % Fit all rotated points
        [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,[],[],[],[],vlb,vub,[],options,[Lcc Mcc],Presps,'surface6','poisson');
        params0(rot,:) = f1;
        GOF(rot) = fval;
        
    end
    
    % Find the best fitting surface
    [~,bestIdx] = min(GOF);
    params1 = params0(bestIdx,:);
    
    % Display surface and pts
    surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface6');
    surface = reshape(surface,size(xx));
    axes(oneDpanel.axes); cla; hold on;
    oneDpanel.surf = surfc(xx,yy,surface);
    alpha(.5); 
    oneDfit.pts = plot3(Lcc,Mcc,Presps,'k*');
    
    % Save parameter values
    oneDpanel.parvals = params1;
    oneDpanel.LL = -GOF(bestIdx);
    loopdata{loop}.oneD.LL = oneDpanel.LL;
    loopdata{loop}.oneD.params = params1;
    
    % Print parameters
    set(oneDpanel.pardisp.A,'string',num2str(rnddec(params1(1),3)));
    set(oneDpanel.pardisp.s,'string',num2str([rnddec(params1(2),3) rnddec(params1(3),3)]));
    set(oneDpanel.pardisp.e,'string',num2str([rnddec(params1(4),3) rnddec(params1(5),3)]));
    set(oneDpanel.pardisp.bl,'string',num2str(rnddec(params1(6),3)));
    set(oneDpanel.pardisp.ang,'string',num2str(rnddec((params1(end)/pi*180),3)));
    set(oneDpanel.pardisp.LL,'string',num2str(rnddec(-GOF(bestIdx),3)));
    set(oneDpanel.pardisp.df,'string',numel(params1));
   

    %%% Fit 2D model %%%
    disp('Fitting a 2D model...')
    
    % Set up some variables
    GOF = nan(nrots,1);
    vlb = [0    .001 .001 .001 .001 .001 .001 .001 .001  0  -pi];
    vub = [1000  50   50   50   50   10   10   10   10  100  pi];
    params0 = nan(nrots,numel(vlb));
    
    % Rotate through angles
    for rot = 1:nrots
        
        % Generating an initial guess
        paramsGuess = [Aguess sigguess sigguess sigguess sigguess...
            expguess expguess expguess expguess blguess angs(rot)];
        
        % Fit all rotated points
        [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,[],[],[],[],vlb,vub,[],options,[Lcc Mcc],Presps,'surface5','poisson');
        params0(rot,:) = f1;
        GOF(rot) = fval;
        
    end
    
    % Find the best fitting surface
    [~,bestIdx] = min(GOF);
    params1 = params0(bestIdx,:);
    
    % Display surface and pts
    surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface5');
    surface = reshape(surface,size(xx));
    axes(twoDpanel.axes); cla; hold on;
    twoDpanel.surf = surfc(xx,yy,surface);
    alpha(.5); 
    twoDpanel.pts = plot3(Lcc,Mcc,Presps,'k*');
    
    % Save parameter values
    twoDpanel.parvals = params1;
    twoDpanel.LL = -GOF(bestIdx);
    loopdata{loop}.twoD.params = params1;
    loopdata{loop}.twoD.LL = twoDpanel.LL;
    
    % Print parameters
    set(twoDpanel.pardisp.A,'string',num2str(rnddec(params1(1),3)));
    set(twoDpanel.pardisp.s,'string',num2str([rnddec(params1(2),3) rnddec(params1(3),3)...
        rnddec(params1(4),3) rnddec(params1(5),3)]));
    set(twoDpanel.pardisp.e,'string',num2str([rnddec(params1(6),3) rnddec(params1(7),3)
        rnddec(params1(8),3) rnddec(params1(9),3)]));
    set(twoDpanel.pardisp.bl,'string',num2str(rnddec(params1(10),3)));
    set(twoDpanel.pardisp.ang,'string',num2str(rnddec((params1(end)/pi*180),3)));
    set(twoDpanel.pardisp.LL,'string',num2str(rnddec(-GOF(bestIdx),3)));
    set(twoDpanel.pardisp.df,'string',numel(params1));
    
    
    %%% Log-likelihood Ratio Test %%%
    
    % Prepare parameters
    nfp1 = numel(oneDpanel.parvals);
    nfp2 = numel(twoDpanel.parvals);
    diffLL = 2 * (twoDpanel.LL - oneDpanel.LL);
    twoDpanel.pval = 1 - chi2cdf(diffLL,nfp2-nfp1);
    
    % Print out and save p-val
    set(twoDpanel.pardisp.pval,'string',num2str(rnddec(twoDpanel.pval,3)));
    loopdata{loop}.pval = twoDpanel.pval;
    
    % Plot histogram of of p-vals
    axes(conpanel.histax); cla; hold on;
    temp = cat(1,loopdata{:});
    pvals = cat(1,temp.pval);
    edges = 0:.025:1;
    n = histc(pvals,edges);
    bar(edges,n)
    
%     % Save User Variables
%     set(ModelFig.ActualPanel,'UserData',actualpanel);
%     set(ModelFig.OneDPanel,'UserData',oneDpanel);
%     set(ModelFig.TwoDPanel,'UserData',twoDpanel);
%     set(gcf,'UserData',ModelFig)
    
    modeldata.loopdata = loopdata;
end

save modeldata4 modeldata
disp('Fin.')

end
