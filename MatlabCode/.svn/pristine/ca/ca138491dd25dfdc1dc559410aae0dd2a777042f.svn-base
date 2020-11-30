function GLMSPopGUI_CardinalHypothesis()

BuildFig()
PopTun()
ScatterPvals()

end


%%% Interactive %%%

function cellselect(~,b)
global GLMSPopData

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

% Record which datafile was selected
conpanel.selectedidx = b.Indices(1);

% Load saved data
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
data = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Cardinal Hypothesis')};
conpanel.surftype = data.surftype;
conpanel.errortype = data.errortype;
oneDfits.sig.lpmpval = data.pvals.oneD.lpm;
oneDfits.sig.lmmpval = data.pvals.oneD.lmm;
oneDfits.LL = data.LL.oneD;
oneDfits.params = data.params.oneD;
twoDfits.sig.lpmpval = data.pvals.twoD.lpm;
twoDfits.sig.lmmpval = data.pvals.twoD.lmm;
twoDfits.LL = data.LL.twoD;
twoDfits.params = data.params.twoD;

% Load data
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);

%%% 1D Surfaces %%%
% Display L+M surface and stim
surface = ComputeNakaRushtonJPW(oneDfits.params.lpm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(oneDfits.axes.lpm); cla; hold on; grid on; axis square;
oneDfits.surf.lpm = surfc(xx,yy,surface);
set(oneDfits.surf.lpm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Display L-M surface and stim
surface = ComputeNakaRushtonJPW(oneDfits.params.lmm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(oneDfits.axes.lmm); cla; hold on; grid on; axis square;
oneDfits.surf.lmm = surfc(xx,yy,surface);
set(oneDfits.surf.lmm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Plot unconstrained surf and stim
surface = ComputeNakaRushtonJPW(oneDfits.params.free,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(oneDfits.axes.freepd); cla; hold on; grid on; axis square;
oneDfits.surf.free = surfc(xx,yy,surface);
set(oneDfits.surf.free(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

%%% 2D Surfaces %%%
% Display L+M surface and stim
surface = ComputeNakaRushtonJPW(twoDfits.params.lpm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(twoDfits.axes.lpm); cla; hold on; grid on; axis square;
twoDfits.surf.lpm = surfc(xx,yy,surface);
set(twoDfits.surf.lpm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Display L-M surface and stim
surface = ComputeNakaRushtonJPW(twoDfits.params.lmm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(twoDfits.axes.lmm); cla; hold on; grid on; axis square;
twoDfits.surf.lmm = surfc(xx,yy,surface);
set(twoDfits.surf.lmm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Plot unconstrained surf and stim
surface = ComputeNakaRushtonJPW(twoDfits.params.free,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(twoDfits.axes.freepd); cla; hold on; grid on; axis square;
twoDfits.surf.free = surfc(xx,yy,surface);
set(twoDfits.surf.free(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

%Replot pop tun and scatterplot
PopTun()
ScatterPvals()

% Highlight histogram box with selected data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

if conpanel.oneDL(conpanel.selectedidx) % 1D hist
    axes(oneDfits.axes.pvalscatter)
    edges = oneDfits.hist.XBinEdges;
    h = histogram2(oneDfits.sig.lpmpval,oneDfits.sig.lmmpval,edges,edges);
    set(h,'FaceColor','none','EdgeColor','none')
    [x,y] = find(h.Values);
    vals = repmat([oneDfits.sig.lpmpval oneDfits.sig.lmmpval],oneDfits.hist.Values(x,y),1);
    h = histogram2(vals(:,1),vals(:,2),edges,edges);
    set(h,'FaceColor',conpanel.oneDcol)
else % 2D hist
    axes(twoDfits.axes.pvalscatter)
    edges = twoDfits.hist.XBinEdges;
    h = histogram2(twoDfits.sig.lpmpval,twoDfits.sig.lmmpval,edges,edges);
    set(h,'FaceColor','none','EdgeColor','none')
    [x,y] = find(h.Values);
    vals = repmat([twoDfits.sig.lpmpval twoDfits.sig.lmmpval],twoDfits.hist.Values(x,y),1);
    h = histogram2(vals(:,1),vals(:,2),edges,edges);
    set(h,'FaceColor',conpanel.twoDcol)
end

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

end

function PopTun()

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

% Organize into 1D and 2D
conpanel.ndthresh = str2double(conpanel.uicontrols.ndthresh.String);
conpanel.oneDL = conpanel.normLLdiff < conpanel.ndthresh;

% Distinguish cells that fail hypothesis from those that don't
conpanel.pvalthresh = str2double(conpanel.uicontrols.pvalthresh.String);
conpanel.pvals.oneD.lpmL = conpanel.pvals.oneD.lpm >= conpanel.pvalthresh;
conpanel.pvals.oneD.lmmL = conpanel.pvals.oneD.lmm >= conpanel.pvalthresh;
conpanel.pvals.twoD.lpmL = conpanel.pvals.twoD.lpm >= conpanel.pvalthresh;
conpanel.pvals.twoD.lmmL = conpanel.pvals.twoD.lmm >= conpanel.pvalthresh;

%%% Plot preferred axes of 1D cells %%%
axes(oneDfits.axes.poptun); cla;

% Plot all 1D tuning
L = conpanel.oneDL;
angs = conpanel.params.oneD.free(L,end-1);
oneDfits.poptunhist.all = rose(angs,20); hold on;
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13)); % Deleting angle labels and outermost rho label
end
set(oneDfits.poptunhist.all,'color',conpanel.oneDcol,'LineWidth',2)

% Plot subs that pass hypothesis test
L = conpanel.oneDL & (conpanel.pvals.oneD.lpmL | conpanel.pvals.oneD.lmmL);
angs = conpanel.params.oneD.free(L,end-1);
oneDfits.poptunhist.pass = rose(angs,20); hold on;
set(oneDfits.poptunhist.pass,'color',[0 0 0],'LineWidth',.25)
axis tight
title('1D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')


%%% Plot preferred axes of 2D cells %%%
axes(twoDfits.axes.poptun); cla;

% Plot all 1D tuning
L = ~conpanel.oneDL;
angs = conpanel.params.twoD.free(L,end-1);
twoDfits.poptunhist.all = rose(angs,20); hold on;
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13)); % Deleting angle labels and outermost rho label
end
set(twoDfits.poptunhist.all,'color',conpanel.twoDcol,'LineWidth',2)

% Plot subs that pass hypothesis test
L = ~conpanel.oneDL & (conpanel.pvals.twoD.lpmL | conpanel.pvals.twoD.lmmL);
angs = conpanel.params.twoD.free(L,end-1);
twoDfits.poptunhist.pass = rose(angs,20); hold on;
set(twoDfits.poptunhist.pass,'color',[0 0 0],'LineWidth',.25)
axis tight
title('1D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

end

function ScatterPvals()

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

% Plot 1D pvals
axes(oneDfits.axes.pvalscatter); cla; hold on; box on;
edges = 0:.05:1;
oneDfits.hist = histogram2(conpanel.pvals.oneD.lpm(conpanel.oneDL),...
    conpanel.pvals.oneD.lmm(conpanel.oneDL),edges,edges);
xlabel('L+M P-value')
ylabel('L-M P-value')
zlabel('Subunit count')

% Plot 2D pvals
axes(twoDfits.axes.pvalscatter); cla; hold on; box on;
edges = 0:.05:1;
twoDfits.hist = histogram2(conpanel.pvals.twoD.lpm(~conpanel.oneDL),...
    conpanel.pvals.twoD.lmm(~conpanel.oneDL),edges,edges);
xlabel('L+M P-value')
ylabel('L-M P-value')
zlabel('Subunit count')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

end

function resetplots(~,~)

PopTun()
ScatterPvals()

end


%%% Workhorse functions %%%

function reanal(~,~)
global GLMSPopData GLMP sub

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'Userdata');
twoDfits = get(cardhypFig.twoDfits,'Userdata');

% Pull out data from pop structure and set up some options/variables
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
surfparams = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Surface Parameters')};
conpanel.surftype = 'conicsection_xy';
conpanel.errortype = 'NegativeBinomial';
conpanel.options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',10.^-9);

% 1D params
oneDfits.params.lpm = [];
oneDfits.params.lmm = [];
oneDfits.params.free = surfparams.oneD.parvals;
oneDfits.LL.lpm = -Inf;
oneDfits.LL.lmm = -Inf;
oneDfits.LL.free = surfparams.oneD.LL;
oneDfits.ub = [max(GLMP.subunit{sub}.nspikes)                            300 300 0 10 max(GLMP.subunit{sub}.nspikes)/2  pi    5];
oneDfits.lb = [max(GLMP.subunit{sub}.nspikes)/2 1/max(GLMP.subunit{sub}.rho)   0 0  1                             .001 -pi .001];

% 2D params
twoDfits.params.lpm = [];
twoDfits.params.lmm = [];
twoDfits.params.free = surfparams.twoD.parvals;
twoDfits.LL.lpm = -Inf;
twoDfits.LL.lmm = -Inf;
twoDfits.LL.free = surfparams.twoD.LL;
twoDfits.ub = [max(GLMP.subunit{sub}.nspikes)                            300 300 300 10 max(GLMP.subunit{sub}.nspikes)/2  pi    5];
twoDfits.lb = [max(GLMP.subunit{sub}.nspikes)/2 1/max(GLMP.subunit{sub}.rho)   0 -10  1                             .001 -pi .001];

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

% Fit surfaces
disp(['Analyzing datafile #' num2str(conpanel.selectedidx) ': ' GLMP.datafile ' sub # ' num2str(sub)])
disp(['Fitting using ' conpanel.surftype])
lpmfit1D()
lmmfit1D()
lpmfit2D()
lmmfit2D()
freefits()
comparefits()
PopTun()
ScatterPvals()

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'Userdata');
twoDfits = get(cardhypFig.twoDfits,'Userdata');

% Save results in Pop structure
disp('Saving results...')
data.surftype = conpanel.surftype;
data.errortype = conpanel.errortype;
data.params.oneD = oneDfits.params;
data.params.twoD = twoDfits.params;
data.pvals.oneD.lpm = oneDfits.sig.lpmpval;
data.pvals.twoD.lpm = twoDfits.sig.lpmpval;
data.pvals.oneD.lmm = oneDfits.sig.lmmpval;
data.pvals.twoD.lmm = twoDfits.sig.lmmpval;
data.LL.oneD = oneDfits.LL;
data.LL.twoD = twoDfits.LL;
GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Cardinal Hypothesis')} = data;
save([conpanel.library 'GLMSPopData'],'GLMSPopData')

disp(['Finished Reanalyzing datafile #' num2str(conpanel.selectedidx) ': ' GLMP.datafile ' sub # ' num2str(sub)])

b.Indices = conpanel.selectedidx;
cellselect([],b)

end

function reanalAll(~,~)
global GLMSPopData

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'Userdata');
twoDfits = get(cardhypFig.twoDfits,'Userdata');

for n = 1:(size(GLMSPopData,1)-1)
    
    conpanel.selectedidx = n;
    set(cardhypFig.conpanel,'userdata',conpanel)
    set(gcf,'userdata',cardhypFig)
    reanal()
    
end

disp('Finished reanalyzing datafiles.')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

end

function lpmfit1D()
global GLMP sub

disp('Fitting a 1D L+M Surface...')

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');

% Load data
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Set initial guesses and options for 1D fit
angs = [pi/4 -3*pi/4];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
vub = oneDfits.ub;
vlb = oneDfits.lb;

% Parameter guesses
L = GLMP.subunit{sub}.theta == angs(1)...
    | GLMP.subunit{sub}.theta == angs(2);
Aguess = max(GLMP.subunit{sub}.nspikes) * .8;
sig1guess = max(GLMP.subunit{sub}.rho(L))/2;
sig2guess = sig1guess*2;
orsigguess = 0;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
for n = 1:numel(angs)
    
    % Customize guess and bounds
    paramsGuess = [Aguess 1/sig1guess 1/sig2guess 1/orsigguess expguess blguess angs(n) kappaguess];
    vlb(end-1) = angs(n);
    vub(end-1) = angs(n);
    
    % Fit all of the data using axis fit as the initial guess
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        a,b,[],[],vlb,vub,[],conpanel.options,[Lcc Mcc],nsp,...
        conpanel.surftype,conpanel.errortype);
        
    % Record results
    if -fval > oneDfits.LL.lpm
        oneDfits.params.lpm = f1;
        oneDfits.LL.lpm = -fval;
    end
    
end

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(oneDfits.params.lpm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(oneDfits.axes.lpm); cla; hold on; grid on; axis square;
oneDfits.surf.lpm = surfc(xx,yy,surface);
set(oneDfits.surf.lpm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(gcf,'userdata',cardhypFig)

end

function lpmfit2D()
global GLMP sub

disp('Fitting a 2D L+M Surface...')

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

% Load data
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Set initial guesses and options for 1D fit
angs = [pi/4 -3*pi/4];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
vub = twoDfits.ub;
vlb = twoDfits.lb;

% Parameter guesses
L = GLMP.subunit{sub}.theta == angs(1)...
    | GLMP.subunit{sub}.theta == angs(2);
Aguess = max(GLMP.subunit{sub}.nspikes) * .8;
sig1guess = max(GLMP.subunit{sub}.rho(L))/2;
sig2guess = sig1guess*2;
orsigguess = 0;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
for n = 1:numel(angs)
    
    % Customize guess and bounds
    paramsGuess = [Aguess 1/sig1guess 1/sig2guess 1/orsigguess expguess blguess angs(n) kappaguess];
    vlb(end-1) = angs(n);
    vub(end-1) = angs(n);
    
    % Fit all of the data using axis fit as the initial guess
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        a,b,[],[],vlb,vub,[],conpanel.options,[Lcc Mcc],nsp,...
        conpanel.surftype,conpanel.errortype);
        
    % Record results
    if -fval > twoDfits.LL.lpm
        twoDfits.params.lpm = f1;
        twoDfits.LL.lpm = -fval;
    end
    
end

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(twoDfits.params.lpm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(twoDfits.axes.lpm); cla; hold on; grid on; axis square;
twoDfits.surf.lpm = surfc(xx,yy,surface);
set(twoDfits.surf.lpm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)


end

function lmmfit1D()
global GLMP sub

disp('Fitting a 1D L-M Surface...')

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');

% Load data
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Set initial guesses and options for 1D fit
angs = [-pi/4 3*pi/4];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
vub = oneDfits.ub;
vlb = oneDfits.lb;

% Parameter guesses
L = GLMP.subunit{sub}.theta == angs(1)...
    | GLMP.subunit{sub}.theta == angs(2);
Aguess = max(GLMP.subunit{sub}.nspikes) * .8;
sig1guess = max(GLMP.subunit{sub}.rho(L))/2;
sig2guess = sig1guess*2;
orsigguess = 0;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
for n = 1:numel(angs)
    
    % Customize guess and bounds
    paramsGuess = [Aguess sig1guess sig2guess orsigguess expguess blguess angs(n) kappaguess];
    vlb(end-1) = angs(n);
    vub(end-1) = angs(n);
    
    % Fit all of the data using axis fit as the initial guess
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        a,b,[],[],vlb,vub,[],conpanel.options,[Lcc Mcc],nsp,...
        conpanel.surftype,conpanel.errortype);
        
    % Record results
    if -fval > oneDfits.LL.lmm
        oneDfits.params.lmm = f1;
        oneDfits.LL.lmm = -fval;
    end
    
end

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(oneDfits.params.lmm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(oneDfits.axes.lmm); cla; hold on; grid on; axis square;
oneDfits.surf.lmm = surfc(xx,yy,surface);
set(oneDfits.surf.lmm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(gcf,'userdata',cardhypFig)

end

function lmmfit2D()
global GLMP sub

disp('Fitting a 2D L-M Surface...')

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

% Load data
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Set initial guesses and options for 1D fit
angs = [-pi/4 3*pi/4];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
vub = twoDfits.ub;
vlb = twoDfits.lb;

% Parameter guesses
L = GLMP.subunit{sub}.theta == angs(1)...
    | GLMP.subunit{sub}.theta == angs(2);
Aguess = max(GLMP.subunit{sub}.nspikes) * .8;
sig1guess = max(GLMP.subunit{sub}.rho(L))/2;
sig2guess = sig1guess*2;
orsigguess = 0;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
if kappaguess < vlb(end)
    kappaguess = vlb(end);
end

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
for n = 1:numel(angs)
    
    % Customize guess and bounds
    paramsGuess = [Aguess 1/sig1guess 1/sig2guess 1/orsigguess expguess blguess angs(n) kappaguess];
    vlb(end-1) = angs(n);
    vub(end-1) = angs(n);
    
    % Fit all of the data using axis fit as the initial guess
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        a,b,[],[],vlb,vub,[],conpanel.options,[Lcc Mcc],nsp,...
        conpanel.surftype,conpanel.errortype);
        
    % Record results
    if -fval > twoDfits.LL.lmm
        twoDfits.params.lmm = f1;
        twoDfits.LL.lmm = -fval;
    end
    
end

% Display surface and pts
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(twoDfits.params.lmm,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));
axes(twoDfits.axes.lmm); cla; hold on; grid on; axis square;
twoDfits.surf.lmm = surfc(xx,yy,surface);
set(twoDfits.surf.lmm(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)


end

function freefits()
global GLMP sub

disp('Calculating free fit values.')

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');

% Pull out variables and saved fit
Lcc = GLMP.subunit{sub}.Lcc;
Mcc = GLMP.subunit{sub}.Mcc;
nsp = GLMP.subunit{sub}.nspikes;
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);

%%% 1D Fit %%%
surface = ComputeNakaRushtonJPW(oneDfits.params.free,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));

% Plot bubble plot and surf
axes(oneDfits.axes.freepd); cla; hold on; grid on; axis square;
oneDfits.surf.free = surfc(xx,yy,surface);
set(oneDfits.surf.free(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

%%% 2D Fit %%%
surface = ComputeNakaRushtonJPW(oneDfits.params.free,[xx(:) yy(:)],conpanel.surftype);
surface = reshape(surface,size(xx));

% Plot bubble plot and surf
axes(twoDfits.axes.freepd); cla; hold on; grid on; axis square;
twoDfits.surf.free = surfc(xx,yy,surface);
set(twoDfits.surf.free(1),'edgecolor','none')
alpha(.5);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')


% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

end

function comparefits()

disp('Comparing fits and calculating P-values...')

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
oneDfits = get(cardhypFig.oneDfits,'UserData');
twoDfits = get(cardhypFig.twoDfits,'UserData');
conpanel.pvalthresh = str2double(conpanel.uicontrols.pvalthresh.String);

% Compare 1D free PD with L+M
diffLL = 2 * (oneDfits.LL.free - oneDfits.LL.lpm); % larger model nLL - smaller model LL
oneDfits.sig.lpmpval = 1 - chi2cdf(diffLL,1);
oneDfits.sig.lpm = oneDfits.sig.lpmpval < conpanel.pvalthresh;

% Compare 1D free PD with L-M
diffLL = 2 * (oneDfits.LL.free - oneDfits.LL.lmm); % larger model nLL - smaller model LL
oneDfits.sig.lmmpval = 1 - chi2cdf(diffLL,1);
oneDfits.sig.lmm = oneDfits.sig.lmmpval < conpanel.pvalthresh;

% Compare 2D free PD with L+M
diffLL = 2 * (twoDfits.LL.free - twoDfits.LL.lpm); % larger model nLL - smaller model LL
twoDfits.sig.lpmpval = 1 - chi2cdf(diffLL,1);
twoDfits.sig.lpm = twoDfits.sig.lpmpval < conpanel.pvalthresh;

% Compare 2D free PD with L-M
diffLL = 2 * (twoDfits.LL.free - twoDfits.LL.lmm); % larger model nLL - smaller model LL
twoDfits.sig.lmmpval = 1 - chi2cdf(diffLL,1);
twoDfits.sig.lmm = twoDfits.sig.lmmpval < conpanel.pvalthresh;

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'UserData',twoDfits)
set(gcf,'userdata',cardhypFig)

end


%%% Setup %%%

function BuildFig()

% Build figure
figure(55588); clf;
set(gcf,'units','normalized','pos',[.05 .1 .9 .8],'numbertitle','off',...
    'name','Population Cardinal PDs Hypothesis Test')

% Set up panels
cardhypFig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .24 .98]);
cardhypFig.oneDfits = uipanel(gcf,'units','normalized','pos',[.26 .01 .35 .98],'title','1D Fits');
cardhypFig.twoDfits = uipanel(gcf,'units','normalized','pos',[.62 .01 .35 .98],'title','2D Fits');

% Control Panel
conpanel.uicontrols.renal = uicontrol('style','pushbutton',...
    'parent',cardhypFig.conpanel,'units','normalized','pos',[.05 .07 .4 .05],...
    'string','Reanalyze Subunit','fontsize',12,'callback',@reanal,...
    'backgroundcolor',[.2 1 0]);
conpanel.uicontrols.renalAll = uicontrol('style','pushbutton',...
    'parent',cardhypFig.conpanel,'units','normalized','pos',[.55 .07 .4 .05],...
    'string','Reanalyze All','fontsize',12,'callback',@reanalAll,...
    'backgroundcolor',[1 .2 0]);
conpanel.uicontrols.overview = uicontrol('parent',cardhypFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .01 .4 .05],'fontsize',12,'callback',@LoadOverview);
conpanel.uicontrols.analGUI = uicontrol('parent',cardhypFig.conpanel,...
    'style','pushbutton','string','Load Control Panel','units','normalized',...
    'pos',[.55 .01 .4 .05],'fontsize',12,'callback',@LoadControlPanel);
conpanel.uicontrols.ndthresh = uicontrol(cardhypFig.conpanel,'style','edit',...
    'units','normalized','pos',[.6 .15 .3 .04],'string',.1,...
    'callback',@resetplots);
conpanel.uicontrols.ndthreshlabel = uicontrol(cardhypFig.conpanel,'style','text',...
    'units','normalized','pos',[.1 .15 .45 .04],'string','Norm LL 1D/2D threshold:');
conpanel.uicontrols.pvalthresh = uicontrol(cardhypFig.conpanel,'style','edit',...
    'units','normalized','pos',[.6 .2 .3 .04],'string',.05,...
    'callback',@resetplots);
conpanel.uicontrols.pvalthreshlabel = uicontrol(cardhypFig.conpanel,'style','text',...
    'units','normalized','pos',[.1 .2 .45 .04],'string','P-value threshold:');
conpanel.selectedidx = [];
conpanel.oneDcol = [1 .5 0];
conpanel.twoDcol = [0 0 1];

%%% Fitspanel %%%
% 1D
oneDfits.axes.freepd = axes('parent',cardhypFig.oneDfits,'units','normalized',...
    'pos',[.05 .6875 .4 .25]); box on; grid on; hold on;
title('Unconstrained PD Fit')
oneDfits.axes.lpm = axes('parent',cardhypFig.oneDfits,'units','normalized',...
    'pos',[.05 .375 .4 .25]); box on; grid on; hold on;
title('L+M Fit')
oneDfits.axes.lmm = axes('parent',cardhypFig.oneDfits,'units','normalized',...
    'pos',[.05 .0625 .4 .25]); box on; grid on; hold on;
title('L-M Fit')
oneDfits.axes.pvalscatter = axes('parent',cardhypFig.oneDfits,'units','normalized',...
    'pos',[.55 .1 .4 .35]); box on; grid on; hold on;
title('1D P-val Scatterplot')
oneDfits.axes.poptun = axes('parent',cardhypFig.oneDfits,'units','normalized',...
    'pos',[.55 .6 .4 .35]); box on;
title('Pop 1D Tuning')

% 2D
twoDfits.axes.freepd = axes('parent',cardhypFig.twoDfits,'units','normalized',...
    'pos',[.05 .6875 .4 .25]); box on; grid on; hold on;
title('Unconstrained PD Fit')
twoDfits.axes.lpm = axes('parent',cardhypFig.twoDfits,'units','normalized',...
    'pos',[.05 .375 .4 .25]); box on; grid on; hold on;
title('L+M Fit')
twoDfits.axes.lmm = axes('parent',cardhypFig.twoDfits,'units','normalized',...
    'pos',[.05 .0625 .4 .25]); box on; grid on; hold on;
title('L-M Fit')
twoDfits.axes.pvalscatter = axes('parent',cardhypFig.twoDfits,'units','normalized',...
    'pos',[.55 .1 .4 .35]); box on; grid on; hold on;
title('2D P-val Scatterplot')
twoDfits.axes.poptun = axes('parent',cardhypFig.twoDfits,'units','normalized',...
    'pos',[.55 .6 .4 .35]); box on;
title('Pop 2D Tuning')

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(cardhypFig.oneDfits,'userdata',oneDfits)
set(cardhypFig.twoDfits,'userdata',twoDfits)
set(gcf,'userdata',cardhypFig)

UnpackPopData

end

function UnpackPopData()
global GLMSPopData

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'UserData');
%oneDfits = get(cardhypFig.oneDfits,'UserData');

% Grab saved population data
if ismac
    conpanel.library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    conpanel.library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([conpanel.library 'GLMSPopData.mat'])

% Unpack data
datatypes = GLMSPopData(1,:);
data = cell(size(GLMSPopData,1)-1,4);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'Tuning'));
data(:,4) = GLMSPopData(2:end,strcmp(datatypes,'normLLdiff'));

% Put some of this data aside, too
conpanel.tuning = [GLMSPopData{2:end,strcmp(datatypes,'Tuning')}]';
conpanel.normLLdiff = [GLMSPopData{2:end,strcmp(datatypes,'normLLdiff')}]';

% Repackage data and name columns
colname = {'Datafile' 'Sub' 'Tuning' 'normLLdiff'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% Display table
conpanel.table = uitable('parent',cardhypFig.conpanel,...
    'units','normalized','pos',[.01 .7  .98 .29],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'BackgroundColor',[1 1 1],...
    'CellSelectionCallback',@cellselect);

% Unpack saved fits
data = [GLMSPopData{2:end,strcmp(datatypes,'Cardinal Hypothesis')}];
conpanel.surftype = data(1).surftype;
conpanel.errortype = data(1).errortype;
t = [data.pvals];
t1 = [t.oneD];
t2 = [t.twoD];
conpanel.pvals.oneD.lpm = [t1.lpm]';
conpanel.pvals.oneD.lmm = [t1.lmm]';
conpanel.pvals.twoD.lpm = [t2.lpm]';
conpanel.pvals.twoD.lmm = [t2.lmm]';
t = [data.LL];
t1 = [t.oneD];
t2 = [t.twoD];
conpanel.LL.oneD.lpm = [t1.lpm]';
conpanel.LL.oneD.lmm = [t1.lmm]';
conpanel.LL.oneD.free = [t1.free]';
conpanel.LL.twoD.lpm = [t2.lpm]';
conpanel.LL.twoD.lmm = [t2.lmm]';
conpanel.LL.twoD.free = [t2.free]';
t = [data.params];
t1 = [t.oneD];
t2 = [t.twoD];
conpanel.params.oneD.lpm = cat(1,t1.lpm);
conpanel.params.oneD.lmm = cat(1,t1.lmm);
conpanel.params.oneD.free = cat(1,t1.free);
conpanel.params.twoD.lpm = cat(1,t2.lpm);
conpanel.params.twoD.lmm = cat(1,t2.lmm);
conpanel.params.twoD.free = cat(1,t2.free);

% Save user data
set(cardhypFig.conpanel,'userdata',conpanel)
set(gcf,'userdata',cardhypFig)

end


%%% Call other functions %%%

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
if isempty(idx)
    disp('Must select a datafile before overview can be shown.')
    return
end
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMSGUI_Overview();

end

function LoadControlPanel(~,~)
global GLMSPopData GLMP DN

% Load user data
cardhypFig = get(gcf,'userdata');
conpanel = get(cardhypFig.conpanel,'userdata');

% Load selected data
datanames = GLMSPopData(1,:);
idx = conpanel.selectedidx+1;
if isempty(idx)
    disp('Must select a datafile before Control Panel can be loaded.')
    return
end
DN = GLMSPopData{idx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{idx,strcmp(datanames,'GLMP')};
GLMS_AnalysisGUI(DN,GLMP);

end

