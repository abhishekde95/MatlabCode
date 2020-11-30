function GLMSPopGUI_Neurometric(~,~)

figure(1992); clf;
set(1992,'pos',[200 200 1200 600],'numbertitle','off','name','Population Tuning');

% Set up panels
NMFig.conpanel = uipanel(1992,'units','normalized','pos',[.01 .51 .24 .48]);
NMFig.dispanel = uipanel(1992,'units','normalized','pos',[.26 .01 .24 .98]);
NMFig.diffpanel = uipanel(1992,'units','normalized','pos',[.51 .01 .48 .98]);
NMFig.neuralpanel = uipanel(1992,'units','normalized','pos',[.01 .01 .24 .49]);

% Set up conpanel
conpanel.uicontrols.reanalyzeall = uicontrol('parent',NMFig.conpanel,...
    'style','pushbutton','string','Reanalyze All Cells','units','normalized',...
    'pos',[.05 .01 .425 .18],'Foregroundcolor',[1 0 0],'callback',@ReanalAll);
conpanel.uicontrols.analindiv = uicontrol('parent',NMFig.conpanel,...
    'style','pushbutton','string','Surface Analysis','units','normalized',...
    'pos',[.525 .01 .425 .18],'callback',@callanalysis);
conpanel.uicontrols.overview = uicontrol('parent',NMFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .2 .425 .18],'callback',@LoadOverview);
conpanel.uicontrols.reanalyze = uicontrol('parent',NMFig.conpanel,...
    'style','pushbutton','string','Reanalyze Cell','units','normalized',...
    'pos',[.525 .2 .425 .18],'callback',@Reanal);
conpanel.table = uitable('parent',NMFig.conpanel,...
    'units','normalized','pos',[.025 .4  .95 .59]);

% Set up dispanel
dispanel.axes.neurometric = axes('parent',NMFig.dispanel,'units','normalized',...
    'pos',[.15 .6 .7 .35]);
dispanel.axes.psychometric = axes('parent',NMFig.dispanel,'units','normalized',...
    'pos',[.15 .1 .7 .35]);
neuralpanel.axes = axes('parent',NMFig.neuralpanel,'units','normalized',...
     'pos',[.15 .1 .7 .8]);
 
% Set up diff panel
diffpanel.axes.pmnmratio = axes('parent',NMFig.diffpanel,'units','normalized',...
    'pos',[.1 .55 .35 .35]);
diffpanel.axes.eightythresh = axes('parent',NMFig.diffpanel,'units','normalized',...
    'pos',[.1 .1 .35 .35]);
diffpanel.axes.proj = axes('parent',NMFig.diffpanel,'units','normalized',...
    'pos',[.55 .55 .35 .35]);
diffpanel.axes.sse = axes('parent',NMFig.diffpanel,'units','normalized',...
    'pos',[.55 .1 .35 .35]);

% Save user data
set(NMFig.conpanel,'userdata',conpanel)
set(NMFig.dispanel,'userdata',dispanel)
set(NMFig.diffpanel,'userdata',diffpanel)
set(NMFig.neuralpanel,'userdata',neuralpanel)
set(1992,'userdata',NMFig)

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

UnpackPopulationData()
DispTable()

end

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
NMFig = get(1992,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.popidx;
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end

function Reanal(~,~)
global GLMSPopData GLMSD GLMP
% Load figure and pop variables
NMFig = get(1992,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

% Which dataset
if isempty(conpanel.psychidx);
    disp('Must Select Dataset First')
    return
end

% Reload nex file, unpack, analyze, and save
datatypes = GLMSPopData(1,:);
neuralfile = conpanel.table.Data(conpanel.psychidx,1);
disp(['Reanalyzing ' char(neuralfile) '...']);
if ismac
    library1 = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
    library2 = ['/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/' char(neuralfile) '/'];
elseif ispc
    library1 = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    library2 = ['C:\Users\jpweller\Dropbox\Patrick\GLMSD Data\' char(neuralfile) '\'];
end

% Unpack psychophysics datasets
% [~,filelist,~] = xlsread([library2 char(neuralfile) '.xlsx']);
% tempdata = [];
% rawdata = [];
% for f = 1:size(filelist,1)
%     datafile = filelist(f,:);
%     if ~strcmp(datafile{1}(1),'%')        
%         tempdata = nex2stro([char(library2) char(datafile)]);
%         rawdata = strocat(rawdata,tempdata);
%     end
%      
%     % Fill in some missing parameters resulting from strocat
%     if f == size(filelist,1)
%        rawdata.sum.exptParams.rf_x = tempdata.sum.exptParams.rf_x;
%        rawdata.sum.exptParams.rf_y = tempdata.sum.exptParams.rf_y;
%        rawdata.sum.exptParams.filename = char(tempdata.sum.exptParams.filename)';
%        rawdata.sum.exptParams.subunit = tempdata.sum.exptParams.subunit;
%     end
%     clear tempdata
% end

% Organize raw data
%GLMSD = OrganizeRawGLMSDData(rawdata);
GLMSD = GLMSPopData{conpanel.popidx,strcmp(datatypes,'GLMSD')};
GLMP = GLMSD.GLMSdata;

% Fit neuro and psycho
%FitPsychometricFun()
FitNeurometricFun()
ComparePsychNeuro()

% Save pop data
GLMSPopData{conpanel.popidx,strcmp(datatypes,'GLMSD')} = GLMSD;
save([library1 'GLMSPopData'],'GLMSPopData')

b.Indices = conpanel.psychidx;
cellselect([],b);

disp(['Datafile ' char(neuralfile) ' Reanalyzed.'])

end

function ReanalAll(~,~,~)

disp('Reanalyzing All Datasets...')

% Load figure and pop variables
NMFig = get(1992,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

% Rotate through selected index
psychidx = find(conpanel.psycho.L);
for n = 1:size(conpanel.table.Data,1)
    conpanel.psychidx = n;
    conpanel.popidx = psychidx(n);
    set(NMFig.conpanel,'userdata',conpanel)
    set(1992,'userdata',NMFig)
    Reanal()
end

end


function UnpackPopulationData()
global GLMSPopData

% Load figure variables
NMFig = get(1992,'userdata');
conpanel = get(NMFig.conpanel,'userdata');
dispanel = get(NMFig.dispanel,'userdata');
datatypes = GLMSPopData(1,:);
neurofiles = GLMSPopData(:,strcmp(datatypes,'Datafile'));

% Which neuro files have psycho data?
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMSD Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMSD Data\';
end
[~,filelist,~] = xlsread([library  'GLMSD_FileNames.xlsx']);
conpanel.psycho.L = ismember(neurofiles,filelist);

% nd params
nd = GLMSPopData(conpanel.psycho.L,strcmp(datatypes,'nD'));
conpanel.neuro.oneD.L = strcmp(nd,'1D');
conpanel.neuro.twoD.L = strcmp(nd,'2D');

% Monkey
filenames = cat(1,GLMSPopData{conpanel.psycho.L,strcmp(datatypes,'Datafile')});
monk = filenames(:,1);
conpanel.monk.M = monk == 'M';
conpanel.monk.N = monk == 'N';

% Save user data
set(NMFig.conpanel,'userdata',conpanel)
set(NMFig.dispanel,'userdata',dispanel)
set(1992,'userdata',NMFig)

end

function DispTable()
global GLMSPopData

% Load figure and pop variables
NMFig = get(1992,'userdata');
conpanel = get(NMFig.conpanel,'userdata');

% Unpack data
datatypes = GLMSPopData(1,:);
data(:,1) = GLMSPopData(conpanel.psycho.L,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(conpanel.psycho.L,strcmp(datatypes,'Tuning'));
data(:,3) = GLMSPopData(conpanel.psycho.L,strcmp(datatypes,'nD'));

% Repackage data and name columns
colname = {'Datafile' 'Tuning' 'nD'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% Display table
set(conpanel.table,'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'BackgroundColor',[1 1 1],...
    'cellselectioncallback',@cellselect);

% Save user data
set(NMFig.conpanel,'userdata',conpanel)
set(1992,'userdata',NMFig)

end


function cellselect(~,b)
global GLMSPopData GLMSD

% Load figure and pop variables
NMFig = get(1992,'userdata');
conpanel = get(NMFig.conpanel,'userdata');
dispanel = get(NMFig.dispanel,'userdata');
diffpanel = get(NMFig.diffpanel,'userdata');
neuralpanel = get(NMFig.neuralpanel,'userdata');

% Set aside the index (+1 for referencing GLMSPopData)
if ~isempty(b.Indices)
    psychidx = find(conpanel.psycho.L);
    conpanel.psychidx = b.Indices(1);
    conpanel.popidx = psychidx(b.Indices(1));
else
    return
end
    
% Unpack variables
datatypes = GLMSPopData(1,:);
surfparams = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Surface Parameters')};
GLMP = GLMSPopData{conpanel.popidx,strcmp(datatypes,'GLMP')};
GLMSD = GLMSPopData{conpanel.popidx,strcmp(datatypes,'GLMSD')};
sub = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Subunit')};
nd = GLMSPopData{conpanel.popidx,strcmp(datatypes,'nD')};
nmparams = GLMSPopData{conpanel.popidx,strcmp(datatypes,'NM Params')};
if strcmp(nd,'1D')
    str = 'oneD';
elseif strcmp(nd,'2D')
    str = 'twoD';
end
AUC = nmparams.(str).AUC;

% Plot Neural Data
params1 = surfparams.(str).parvals;
surftype = surfparams.(str).surftype;
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],surftype);
surface = reshape(surface,size(xx));
axes(neuralpanel.axes); cla; hold on; grid on;
p = surfc(xx,yy,surface);
set(p,'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);
plot3(GLMP.subunit{sub}.Lcc,GLMP.subunit{sub}.Mcc,GLMP.subunit{sub}.nspikes,'k*');
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
title('Fit to Neural Data')

% Plot neurometric surface
params2 = nmparams.(str).parvals;
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params2,[xx(:) yy(:)],surftype);
surface = reshape(surface,size(xx));
axes(dispanel.axes.neurometric); cla; hold on; grid on;
p = surfc(xx,yy,surface);
set(p,'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);
plot3(GLMP.subunit{sub}.uniqueLcc,GLMP.subunit{sub}.uniqueMcc,AUC,'k*');
xlabel('Lcc');
ylabel('Mcc');
zlabel('Prediction of Correct Correct Answer')
title('Neurometric Fit')

% Plot psychometric surface
params3 = GLMSD.fit.parvals;
surftype = GLMSD.fit.surftype;
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params3,[xx(:) yy(:)],surftype);
surface = reshape(surface,size(xx));
axes(dispanel.axes.psychometric); cla; hold on; grid on;
p = surfc(xx,yy,surface);
set(p,'edgecolor','none')
alpha(.5);
set(gca,'view',[0 90]);
plot3(GLMSD.uniqueLcc,GLMSD.uniqueMcc,GLMSD.pCorrect,'k*');
xlabel('Lcc');
ylabel('Mcc');
zlabel('Percentage Correct')
title('Psychometric Fit')

% Plot Psych to Neuro Ratio
thetas = unique(GLMSD.theta);
[maxval,maxidx] = max(nmparams.nmpmIntegralRatio);
axes(diffpanel.axes.pmnmratio); cla; 
t = 0 : .01 : 2 * pi; % this is a trick to control ax lims
P = polar(t, maxval * ones(size(t)));
set(P, 'Visible', 'off'); hold on;
polar([thetas; thetas(1)],[nmparams.nmpmIntegralRatio; nmparams.nmpmIntegralRatio(1)],'k-')
polar([params2(end) params2(end)],[0 maxval],'g--')
polar([thetas(maxidx) thetas(maxidx)],[0 maxval],'r--')
xlabel('Lcc')
ylabel('Mcc')
title('Neurometric : Psychometric Integral Ratio')

% Plot 80% Threshold
axes(diffpanel.axes.eightythresh); cla;
[maxval,maxidx] = max(nmparams.nmpm80ThreshRatio);
t = 0 : .01 : 2 * pi; % this is a trick to control ax lims
P = polar(t, maxval * ones(size(t)));
set(P, 'Visible', 'off'); hold on;
polar([thetas; thetas(1)],[nmparams.nmpm80ThreshRatio; nmparams.nmpm80ThreshRatio(1)],'k-')
polar([params2(end) params2(end)],[0 maxval],'g--')
polar([thetas(maxidx) thetas(maxidx)],[0 maxval],'r--')
xlabel('Lcc')
ylabel('Mcc')
title('Psychometric : Neurometric 80% Thresh Ratio')

% Plot projections of unit vectors
[maxval,maxidx] = max(nmparams.nmpmproj);
axes(diffpanel.axes.proj); cla;
t = 0 : .01 : 2 * pi; % this is a trick to control ax lims
P = polar(t, maxval * ones(size(t)));
set(P, 'Visible', 'off'); hold on;
polar([thetas; thetas(1)],[nmparams.nmpmproj; nmparams.nmpmproj(1)],'k-')
polar([params2(end) params2(end)],[0 maxval],'g--')
polar([thetas(maxidx) thetas(maxidx)],[0 maxval],'r--')
xlabel('Lcc')
ylabel('Mcc')
title('Neurometric : Psychometric Correlation')

% SSE
[~,minidx] = min(nmparams.nmpmsse);
maxval = max(nmparams.nmpmsse);
axes(diffpanel.axes.sse); cla
t = 0 : .01 : 2 * pi; % this is a trick to control ax lims
P = polar(t, maxval * ones(size(t)));
set(P, 'Visible', 'off'); hold on;
polar([thetas; thetas(1)],[nmparams.nmpmsse; nmparams.nmpmsse(1)],'k-')
hold on;
polar([params2(end) params2(end)],[0 maxval],'g--')
polar([thetas(minidx) thetas(minidx)],[0 maxval],'r--')
xlabel('Lcc')
ylabel('Mcc')
title('Neurometric : Psychometric SSE')

% Save user data
set(NMFig.conpanel,'userdata',conpanel)
set(NMFig.dispanel,'userdata',dispanel)
set(NMFig.neuralpanel,'userdata',neuralpanel)
set(1992,'userdata',NMFig)

end

function FitNeurometricFun()
global GLMSPopData GLMP

disp('Fitting a function to Neurometric Data...')

% Load Figure Variables
NMFig = get(1992,'UserData');
conpanel = get(NMFig.conpanel,'UserData');
surfparams = GLMSPopData{conpanel.popidx,strcmp(GLMSPopData(1,:),'Surface Parameters')};
sub = GLMSPopData{conpanel.popidx,strcmp(GLMSPopData(1,:),'Subunit')};
AUC = GLMP.subunit{sub}.AUC;

% Fit 1D Neurometric Functions
[parvals,LL,hessian] = GLMSGUI_AUC(sub,surfparams.oneD.surftype,surfparams.oneD.parvals);
oneD.surftype = surfparams.oneD.surftype;
oneD.parvals = parvals;
oneD.AUC = AUC;
oneD.LL = LL;
oneD.Hessian = hessian;

% Fit 2D Neurometric Functions
[parvals,LL,hessian] = GLMSGUI_AUC(sub,surfparams.twoD.surftype,surfparams.twoD.parvals);
twoD.surftype = surfparams.twoD.surftype;
twoD.parvals = parvals;
twoD.AUC = AUC;
twoD.LL = LL;
twoD.Hessian = hessian;

% Save population data
GLMSPopData{conpanel.popidx,strcmp(GLMSPopData(1,:),'NM Params')} = struct('oneD',oneD,'twoD',twoD);

% Save user varaibles
set(NMFig.conpanel,'UserData',conpanel)
set(1992,'UserData',NMFig)

end



function FitPsychometricFun()
global GLMSPopData GLMSD

disp('Fitting a 2D function to Psychometric Data...')

% Load Figure Variables
NMFig = get(1992,'UserData');
conpanel = get(NMFig.conpanel,'UserData');
dispanel = get(NMFig.dispanel,'UserData');
diffpanel = get(NMFig.diffpanel,'UserData');

% Set up x,y,z vals
xvals = GLMSD.Lcc;
yvals = GLMSD.Mcc;
zvals = GLMSD.AnsCorrect;
conpanel.pcorrect = zvals;
    
% Set up some variables
angs = linspace(0,pi,5);
angs(end) = [];
vlb = [.5 0.001 0.001 0.001 0.001 .001  .5  -pi];
vub = [1   10    10    10    10    10   .5   pi];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
Aguess = 1;
sigguess = [.04 0.4];
expguess = 2;
blguess = .5;
guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(sigguess) numel(sigguess) numel(angs)]);
GOF = nan(size(guessIdx,1),1);
params = nan(size(guessIdx,1),numel(vlb));
surftype = 'surface8';
errortype = 'bernoulli';
    
% Rotate through angles
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    params0 = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,2))...
        sigguess(guessIdx(rot,3)) sigguess(guessIdx(rot,4))...
        expguess blguess angs(guessIdx(rot,5))];
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,surftype,errortype);
    params(rot,:) = f1;
    GOF(rot) = fval;
    
end
    
% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Display and save stats
GLMSD.fit.parvals = params1;
GLMSD.fit.LL = -GOF(bestIdx);
GLMSD.fit.surftype = surftype;

% Save population data
GLMSPopData{conpanel.popidx,strcmp(GLMSPopData(1,:),'GLMSD')} = GLMSD;

% Save user varaibles
set(NMFig.conpanel,'UserData',conpanel)
set(NMFig.dispanel,'UserData',dispanel)
set(NMFig.diffpanel,'UserData',diffpanel)
set(1992,'UserData',NMFig)

end


function ComparePsychNeuro()
global GLMSPopData GLMSD

% Load Figure Variables
NMFig = get(1992,'UserData');
conpanel = get(NMFig.conpanel,'UserData');
dispanel = get(NMFig.dispanel,'UserData');
diffpanel = get(NMFig.diffpanel,'UserData');

% Unpack population variables
datatypes = GLMSPopData(1,:);
nd = GLMSPopData{conpanel.popidx,strcmp(datatypes,'nD')};
nmparams = GLMSPopData{conpanel.popidx,strcmp(datatypes,'NM Params')};
if strcmp(nd,'1D')
   str = 'oneD';
elseif strcmp(nd,'2D')
    str = 'twoD';
end
surftype = nmparams.(str).surftype;
nmparvals = nmparams.(str).parvals;
pmparvals = GLMSD.fit.parvals;
thetas = unique(GLMSD.theta);
diffpanel.thetas = thetas;
diffpanel.pmintegral = nan(size(thetas));
diffpanel.nmintegral = nan(size(thetas));
diffpanel.pm80thresh = nan(size(thetas));
diffpanel.nm80thresh = nan(size(thetas));
diffpanel.nmpmproj = nan(size(thetas));
diffpanel.nmpmsse = nan(size(thetas));

% Rotate through each color direction and fit function
for n = 1:numel(thetas)
   
    % Psychometric and neurometric fit responses 
    pmIdx = find(GLMSD.theta == thetas(n));
    [~,maxidx] = max(abs(GLMSD.Lcc(pmIdx)));
    maxLcc = GLMSD.Lcc(pmIdx(maxidx));
    [~,maxidx] = max(abs(GLMSD.Mcc(pmIdx)));
    maxMcc = GLMSD.Mcc(pmIdx(maxidx));
    x = linspace(0,maxLcc,100)';
    y = linspace(0,maxMcc,100)';
    pmresp = ComputeNakaRushtonJPW(pmparvals,[x y],'surface8');
    nmresp = ComputeNakaRushtonJPW(nmparvals,[x y],surftype);
    
    % Integrate under the curve
    diffpanel.pmintegral(n) = trapz(pmresp-.5);
    diffpanel.nmintegral(n) = trapz(nmresp-.5);
    
    % Find the contrast that results in 80% correct
    [~,r] = cart2pol(x,y);
    diffpanel.pm80thresh(n) = interp1(pmresp,r,.8,'spline');
    if (nmresp(end) - nmresp(1)) > .01
        diffpanel.nm80thresh(n) = interp1(nmresp,r,.8,'spline');
    else
        diffpanel.nm80thresh(n) = 100;
    end
    
    % Projection
    pmlen = sqrt(sum((pmresp-.5).^2));
    normpmresp = (pmresp-.5)./pmlen;
    normnmresp = (nmresp-.5)./pmlen;
    diffpanel.nmpmproj(n) = normnmresp' * normpmresp;
    
    % SSE
    diffpanel.nmpmsse(n) = sum((pmresp - nmresp).^2);
    
end
diffpanel.nm2pmIntegralRatio = diffpanel.nmintegral ./ diffpanel.pmintegral;
diffpanel.pm2nm80ThreshRatio = diffpanel.pm80thresh ./ diffpanel.nm80thresh;

% Save population data
nmparams.pmintegral = diffpanel.pmintegral;
nmparams.nmintegral = diffpanel.nmintegral;
nmparams.nmpmIntegralRatio = diffpanel.nm2pmIntegralRatio;
nmparams.pm80thresh = diffpanel.pm80thresh;
nmparams.nm80thresh = diffpanel.nm80thresh;
nmparams.nmpm80ThreshRatio = diffpanel.pm2nm80ThreshRatio;
nmparams.nmpmproj = diffpanel.nmpmproj;
nmparams.nmpmsse = diffpanel.nmpmsse;
GLMSPopData{conpanel.popidx,strcmp(datatypes,'NM Params')} = nmparams;

% Save user varaibles
set(NMFig.conpanel,'UserData',conpanel)
set(NMFig.dispanel,'UserData',dispanel)
set(NMFig.diffpanel,'UserData',diffpanel)
set(1992,'UserData',NMFig)

end

