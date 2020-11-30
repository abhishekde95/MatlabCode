function GLMSPopGUI_ConfidenceIntervals()
% Brute force testing of surfaces with fixed angles to determine confidence intervals

% Main script
setupFig()
resetAll()

end

function resetAll(~,~)

UnpackPopData()
exclusioncriteria()
PlotTuningHists()
%roseHistConfs
scatterplotConf()

end


%%% Workhorse functions %%%

function roseHistConfs()

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
fitspanel = get(confintFig.fitspanel,'UserData');
surfpanel = get(confintFig.surfpanel,'UserData');
poppanel = get(confintFig.poppanel,'userdata');

x = linspace(-pi,pi,poppanel.nbins*2+1);
angbins = x(1:2:end);
angcents = x(2:2:end)';
angfits = poppanel.params(poppanel.oneD.L,end-1);
[~,angidx] = histc(angfits,angbins);

confs = poppanel.bfconfs(poppanel.oneD.L);
confmeans = nan(size(angcents))';
confstds = nan(size(angcents))';
for n = 1:poppanel.nbins
    L = angidx == n;
    confmeans(n) = mean(confs(L));
    confstds(n) = std(confs(L));
end

% Generate smoothing gaussian
binsize = 30/180*pi; % in rad
bins = -pi:binsize:pi;
gaussSize = 20/180*pi; % in rad
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing

% Bin and smooth PDs
angPSTH = histc(angfits,bins);
angPSTH = (angPSTH./numel(angfits)) ./ binsize;
smoothangPSTH = conv(angPSTH,gaussfilter,'same');

% bin and smooth confidence intervals and stds
smoothconfmeans = conv(confmeans,gaussfilter,'same');
smoothconfstds = conv(confstds,gaussfilter,'same');

figure(444); clf;
polar(angcents,smoothangPSTH,'k-'); hold on;
polar(angfits,confs,'ro')


keyboard
figure(445); clf; hold on;
plot(angcents',confmeans,'ko--')
%shadedErrorBar(angcents,confmeans,confstds)
%shadedErrorBar(bins',smoothangPSTH)
plot(angfits,confs,'ro')
xlim([-pi pi])
ylim([0 max(confmeans+confstds)])

figure(45); clf; 
polar(cat(1,angcents,angcents(1)),cat(1,confmeans,confmeans(1)),'k-'); hold on;
polar(cat(1,angcents,angcents(1)),cat(1,confmeans-confstds,confmeans(1)-confstds(1)),'r')
polar(cat(1,angcents,angcents(1)),cat(1,confmeans+confstds,confmeans(1)+confstds(1)),'r')



end

function analyzeCell(~,~)
global GLMSPopData GLMP sub

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
poppanel = get(confintFig.poppanel,'userdata');
fitspanel = get(confintFig.fitspanel,'UserData');
surfpanel = get(confintFig.surfpanel,'Userdata');

% If no datafile is selected, punt
if isempty(conpanel.selectedidx)
    disp('Must select datafile before analysis can be performed...')
    return
else
    datatypes = GLMSPopData(1,:);
    GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
    disp(['Calculating confidence intervals on datafile #' num2str(conpanel.selectedidx)...
        ': ' GLMP.datafile ' sub #' num2str(sub) '...'])
end

% Load pop data
fitspanel.params1D = poppanel.oneD.params(conpanel.selectedidx,:);
fitspanel.params2D = poppanel.twoD.params(conpanel.selectedidx,:);

disp('Fitting 1D and 2D surfaces...')

% Set up some fitting specifics
fitspanel.ub.oneD.uc = [max(GLMP.subunit{sub}.nspikes)                            300 0 0 10 max(GLMP.subunit{sub}.nspikes)/2  pi    5];
fitspanel.lb.oneD.uc = [max(GLMP.subunit{sub}.nspikes)/2 1/max(GLMP.subunit{sub}.rho) 0 0  1                             .001 -pi .001];
fitspanel.ub.twoD.uc = [max(GLMP.subunit{sub}.nspikes)                            300 0 300 10 max(GLMP.subunit{sub}.nspikes)/2  pi    5];
fitspanel.lb.twoD.uc = [max(GLMP.subunit{sub}.nspikes)/2 1/max(GLMP.subunit{sub}.rho) 0 -10  1                             .001 -pi .001];
fitspanel.ub.oneD.bc = [max(GLMP.subunit{sub}.nspikes)                            300 300 0 10 max(GLMP.subunit{sub}.nspikes)/2  pi    5];
fitspanel.lb.oneD.bc = [max(GLMP.subunit{sub}.nspikes)/2 1/max(GLMP.subunit{sub}.rho)   0 0  1                             .001 -pi .001];
fitspanel.ub.twoD.bc = [max(GLMP.subunit{sub}.nspikes)                            300 300 300 10 max(GLMP.subunit{sub}.nspikes)/2  pi    5];
fitspanel.lb.twoD.bc = [max(GLMP.subunit{sub}.nspikes)/2 1/max(GLMP.subunit{sub}.rho)   0 -10  1                             .001 -pi .001];
fitspanel.params.oneD.uc = nan(fitspanel.nangs,numel(fitspanel.ub.oneD.uc));
fitspanel.params.twoD.uc = nan(fitspanel.nangs,numel(fitspanel.ub.twoD.uc));
fitspanel.params.oneD.bc = nan(fitspanel.nangs,numel(fitspanel.ub.oneD.bc));
fitspanel.params.twoD.bc = nan(fitspanel.nangs,numel(fitspanel.ub.twoD.bc));
fitspanel.LL.oneD.uc = ones(fitspanel.nangs,1)*-Inf;
fitspanel.LL.twoD.uc = ones(fitspanel.nangs,1)*-Inf;
fitspanel.LL.oneD.bc = ones(fitspanel.nangs,1)*-Inf;
fitspanel.LL.twoD.bc = ones(fitspanel.nangs,1)*-Inf;


% Load data
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Set initial guesses and options for 1D fit
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);

% Incorporate best fit params
[~,startrot] = min(abs(fitspanel.angs-fitspanel.params1D(end-1)));
angnums = (1:numel(fitspanel.angs))';
angnums = circshift(angnums,-(startrot-1));

% Other params guesses
Aguess = max(GLMP.subunit{sub}.nspikes) * .8;
sig1guess = 1/max(GLMP.subunit{sub}.rho)/2;
sig2guess = sig1guess/2;
orsigguess = Inf;
expguess = 3;
blguess = mean(GLMP.subunit{sub}.blnspikes);
kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
    GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
paramsGuess = [Aguess sig1guess sig2guess orsigguess expguess blguess 0 kappaguess];

% Rotate through axes for which we have data and fit each one. This serves
% as the initial guess for the larger surface.
for n = 1:numel(angnums)
    rot = angnums(n);
    
    strs = {'uc' 'bc'};
    disp(['Fitting surface # ' num2str(n) ' of ' num2str(numel(angnums))])
    for i = 1:2
        
        %%% 1D Fits %%%
        % Bound at this angle
        vub = fitspanel.ub.oneD.(strs{i});
        vlb = fitspanel.lb.oneD.(strs{i});
        vub(end-1) = fitspanel.angs(rot);
        vlb(end-1) = fitspanel.angs(rot);
        
        % Generic guess
        [x,y] = pol2cart(fitspanel.angs(rot),1);
        proj = [x y] * [Lcc Mcc]';
        paramsGuess(2) = 1/max(proj)/2;
        if i == 1
            paramsGuess(3) = Inf;
        else
            paramsGuess(3) = paramsGuess(2)/2;
        end
        paramsGuess(4) = Inf;
        paramsGuess(end-1) = fitspanel.angs(rot);
        [params,nll,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
            fitspanel.surftype,fitspanel.errortype);
        if -nll > fitspanel.LL.oneD.(strs{i})(rot)
            fitspanel.params.oneD.(strs{i})(rot,:) = params;
            fitspanel.LL.oneD.(strs{i})(rot) = -nll;
            fitspanel.Hess.oneD.(strs{i}){rot} = hess;
        end
        
        % Use previous fit as guess
        if n == 1
            paramsGuess = fitspanel.params1D;
        else
            paramsGuess = fitspanel.params.oneD.(strs{i})(angnums(n-1),:);
        end
        paramsGuess(end-1) = fitspanel.angs(rot);
        [params,nll,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
            fitspanel.surftype,fitspanel.errortype);
        if -nll > fitspanel.LL.oneD.(strs{i})(rot)
            fitspanel.params.oneD.(strs{i})(rot,:) = params;
            fitspanel.LL.oneD.(strs{i})(rot) = -nll;
            fitspanel.Hess.oneD.(strs{i}){rot} = hess;
        end
        
        
        %%% 2D Fits %%%
        % Bound at this angle
        vub = fitspanel.ub.twoD.(strs{i});
        vlb = fitspanel.lb.twoD.(strs{i});
        vub(end-1) = fitspanel.angs(rot);
        vlb(end-1) = fitspanel.angs(rot);
        
        % Generic guess
        [x,y] = pol2cart(fitspanel.angs(rot),1);
        proj = [x y] * [Lcc Mcc]';
        sig1guess = 1/max(proj)/2;
        paramsGuess(2) = sig1guess;
        if i == 1
            paramsGuess(3) = 0;
        else
            paramsGuess(3) = sig1guess/2;
        end
        paramsGuess(4) = 1/.5;
        paramsGuess(end-1) = fitspanel.angs(rot);
        [params,nll,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
            fitspanel.surftype,fitspanel.errortype);
        if -nll > fitspanel.LL.twoD.(strs{i})(rot)
            fitspanel.params.twoD.(strs{i})(rot,:) = params;
            fitspanel.LL.twoD.(strs{i})(rot) = -nll;
            fitspanel.Hess.twoD.(strs{i}){rot} = hess;
        end
        
        % Test best 1D fit on 2D constraints
        paramsGuess = fitspanel.params.oneD.(strs{i})(rot,:);
        [params,nll,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
            fitspanel.surftype,fitspanel.errortype);
        if -nll > fitspanel.LL.twoD.(strs{i})(rot)
            fitspanel.params.twoD.(strs{i})(rot,:) = params;
            fitspanel.LL.twoD.(strs{i})(rot) = -nll;
            fitspanel.Hess.twoD.(strs{i}){rot} = hess;
        end
        
        % Use previous 2D fit
        if n ~= 1
            paramsGuess = fitspanel.params.twoD.(strs{i})(angnums(n-1),:);
            if i == 1
                paramsGuess(3) = 0;
            end
            [params,nll,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
                a,b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
                fitspanel.surftype,fitspanel.errortype);
            if -nll > fitspanel.LL.twoD.(strs{i})(rot)
                fitspanel.params.twoD.(strs{i})(rot,:) = params;
                fitspanel.LL.twoD.(strs{i})(rot) = -nll;
                fitspanel.Hess.twoD.(strs{i}){rot} = hess;
            end
        end
        
        %%% 1D Recheck %%%
        
        % Bound at this angle
        vub = fitspanel.ub.oneD.(strs{i});
        vlb = fitspanel.lb.oneD.(strs{i});
        vub(end-1) = fitspanel.angs(rot);
        vlb(end-1) = fitspanel.angs(rot);
        
        % Test 2D fit on 1D constraints
        paramsGuess = fitspanel.params.twoD.(strs{i})(rot,:);
        if i == 1
            paramsGuess(3) = 0;
        end
        paramsGuess(4) = 0;
        [params,nll,~,~,~,~,hess] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],vlb,vub,[],fitspanel.options,[Lcc Mcc],nsp,...
            fitspanel.surftype,fitspanel.errortype);
        if -nll > fitspanel.LL.oneD.(strs{i})(rot)
            fitspanel.params.oneD.(strs{i})(rot,:) = params;
            fitspanel.LL.oneD.(strs{i})(rot) = -nll;
            fitspanel.Hess.oneD.(strs{i}){rot} = hess;
        end
        
        % Calculate confidence intervals
        oneDbound = max(fitspanel.LL.oneD.(strs{i})) - chi2inv(.95,1)/2;
        twoDbound = max(fitspanel.LL.twoD.(strs{i})) - chi2inv(.95,1)/2;
        fitspanel.bfconfdeg.oneD.(strs{i}) = sum(fitspanel.LL.oneD.(strs{i}) >= oneDbound) * (360/(numel(angnums)-1));
        fitspanel.bfconfdeg.twoD.(strs{i}) = sum(fitspanel.LL.twoD.(strs{i}) >= twoDbound) * (360/(numel(angnums)-1));
        [~,oneDidx] = max(fitspanel.LL.oneD.(strs{i}));
        [~,twoDidx] = max(fitspanel.LL.twoD.(strs{i}));
        fitspanel.hessconfdeg.oneD.(strs{i}) = sqrt(1/fitspanel.Hess.oneD.(strs{i}){oneDidx}(end-1,end-1))/pi*180;
        fitspanel.hessconfdeg.twoD.(strs{i}) = sqrt(1/fitspanel.Hess.twoD.(strs{i}){twoDidx}(end-1,end-1))/pi*180;
        
    end
    
end

% Save user data
set(confintFig.conpanel,'userdata',conpanel)
set(confintFig.poppanel,'userdata',poppanel)
set(confintFig.fitspanel,'userdata',fitspanel)
set(confintFig.surfpanel,'userdata',surfpanel)
set(gcf,'userdata',confintFig)

% Save results in pop data
datatypes = GLMSPopData(1,:);
data.angs = fitspanel.angs;
data.params = fitspanel.params;
data.LL = fitspanel.LL;
data.Hess = fitspanel.Hess;
data.hessconfdeg = fitspanel.hessconfdeg;
data.bfconfdeg = fitspanel.bfconfdeg;
load([conpanel.library 'GLMSPopData'],'GLMSPopData')
GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Confidence Intervals')} = data;
save([conpanel.library 'GLMSPopData'],'GLMSPopData')

resetAll()

disp(['Finishied fitting confidence intervals on datafile #' num2str(conpanel.selectedidx) ':'...
    GLMP.datafile 'sub #' num2str(sub) '.'])


end

function ReanalAll(~,~)
global GLMSPopData

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
%fitspanel = get(confintFig.fitspanel,'UserData');
%surfpanel = get(confintFig.surfpanel,'UserData');

for n = 2%:(size(GLMSPopData,1)-1)
    
    conpanel.selectedidx = n;
    set(confintFig.conpanel,'UserData',conpanel)
    analyzeCell()
    
end

end

function scatterplotConf()

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
poppanel = get(confintFig.poppanel,'userdata');

% Plot scatterplot
% axes(conpanel.axes.scatterplot); cla; hold on; grid on;
% plot(poppanel.oneD.bfconfs(poppanel.oneD.L),poppanel.oneD.hessconfs(poppanel.oneD.L),...
%     'o','color',conpanel.oneDcol)
% plot(poppanel.oneD.bfconfs(~poppanel.oneD.L),poppanel.twoD.hessconfs(~poppanel.oneD.L),...
%     'o','color',conpanel.twoDcol)
% xlabel('Brute Force (deg)')
% ylabel('Hessian (deg)')
% xlim([0 360])
% ylim([0 360])

% axes(conpanel.axes.scatterplot)
% if poppanel.diffnormLL(conpanel.selectedidx) < conpanel.normLLthresh
%     plot(poppanel.oneD.bfconfs(conpanel.selectedidx),poppanel.oneD.hessconfs(conpanel.selectedidx),...
%         '*','color',conpanel.oneDcol,'markersize',12,'linewidth',2)
% else
%     plot(poppanel.twoD.bfconfs(conpanel.selectedidx),poppanel.twoD.hessconfs(conpanel.selectedidx),...
%         '*','color',conpanel.twoDcol,'markersize',12,'linewidth',2)
% end

% Save user data
set(confintFig.conpanel,'userdata',conpanel)
set(confintFig.poppanel,'userdata',poppanel)
set(gcf,'userdata',confintFig)

end


%%% Interactive functions %%%

function cellselect(~,b)
global GLMSPopData GLMP sub

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
fitspanel = get(confintFig.fitspanel,'UserData');
surfpanel = get(confintFig.surfpanel,'UserData');
poppanel = get(confintFig.poppanel,'userdata');

% Find selected dataset and display
if isfield(b,'Indices')
    conpanel.selectedidx = b.Indices(1);
end
conpanel.selectedidx = b.Indices(1);
datatypes = GLMSPopData(1,:);

% Pull out saved data and put into figure structure
data = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Confidence Intervals')};
fitspanel.angs = data.angs;
fitspanel.params = data.params;
fitspanel.LL = data.LL;
fitspanel.Hess = data.Hess;
fitspanel.hessconfdeg = data.hessconfdeg;
fitspanel.bfconfdeg = data.bfconfdeg;

% Plot LLs
axes(fitspanel.axes); cla; hold on;
if poppanel.oneD.unichrom.L(conpanel.selectedidx)
    ll1 = fitspanel.LL.oneD.uc;
    bound1 = max(fitspanel.LL.oneD.uc) - chi2inv(.95,1)/2;
    disp('1D fits are unipolar')
else
    ll1 = fitspanel.LL.oneD.bc;
    bound1 = max(fitspanel.LL.oneD.bc) - chi2inv(.95,1)/2;
    disp('1D fits are bipolar')
end
if poppanel.twoD.unichrom.L(conpanel.selectedidx)
    ll2 = fitspanel.LL.twoD.uc;
    bound2 = max(fitspanel.LL.twoD.uc) - chi2inv(.95,1)/2;
    disp('2D fits are unipolar')
else
    ll2 = fitspanel.LL.twoD.bc;
    bound2 = max(fitspanel.LL.twoD.bc) - chi2inv(.95,1)/2;
    disp('2D fits are bipolar')
end
plot([-pi pi],[bound1 bound1],'--','color',conpanel.oneDcol)
plot([-pi pi],[bound2 bound2],'--','color',conpanel.twoDcol)
p1 = plot(fitspanel.angs,ll1,'o-','color',conpanel.oneDcol);
p2 = plot(fitspanel.angs,ll2,'o-','color',conpanel.twoDcol);
set(p1,'ButtonDownFcn',@surfselect)
set(p2,'ButtonDownFcn',@surfselect)
xlim([-pi pi])

axes(surfpanel.axes.oneD); cla; box on; grid on; hold on;
title('1D Fit')
axes(surfpanel.axes.twoD); cla; box on; grid on; hold on;
title('2D Fit')

% Plot bubble plot
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% 1D bubble plot
axes(surfpanel.axes.oneD); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mn,'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% 2D bubble plot
axes(surfpanel.axes.twoD); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mn,'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')

% Save user data
set(confintFig.conpanel,'userdata',conpanel)
set(confintFig.fitspanel,'userdata',fitspanel)
set(confintFig.surfpanel,'userdata',surfpanel)
set(gcf,'userdata',confintFig)



scatterplotConf()
%rosehistconfs()

end

function surfselect(~,b)
global GLMSPopData

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
fitspanel = get(confintFig.fitspanel,'UserData');
surfpanel = get(confintFig.surfpanel,'Userdata');
poppanel = get(confintFig.poppanel,'userdata');

% Load pop data
datatypes = GLMSPopData(1,:);
GLMP = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{conpanel.selectedidx+1,strcmp(datatypes,'Subunit')};
Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);

% Which surface was selected?
[~,surfidx] = min(abs(fitspanel.angs - b.IntersectionPoint(1)));

% Is the selected cell unchrom or bichrom?
if poppanel.oneD.unichrom.L(conpanel.selectedidx)
    str1d = 'uc';
else
    str1d = 'bc';
end
if poppanel.twoD.unichrom.L(conpanel.selectedidx)
    str2d = 'uc';
else
    str2d = 'bc';
end

% Display 1d surface and pts
params1 = fitspanel.params.oneD.(str1d)(surfidx,:);
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(surfpanel.axes.oneD); cla; hold on; grid on;
ticks = linspace(min(Lcc),max(Lcc),5);
set(gca,'XTick',ticks,'YTick',ticks,'xlim',[min(xx(:)) max(xx(:))],'ylim',[min(yy(:)) max(yy(:))]);
surfpanel.surf.oneD = surf(gca,xx,yy,surface);
set(surfpanel.surf.oneD(1),'edgecolor','none')
alpha(.3);
contour3(xx,yy,surface,'linewidth',2);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mn,'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
colormap('cool')

% Display 2d surface and pts
params2 = fitspanel.params.twoD.(str2d)(surfidx,:);
x = linspace(-max(GLMP.subunit{sub}.rho),max(GLMP.subunit{sub}.rho),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params2,[xx(:) yy(:)],fitspanel.surftype);
surface = reshape(surface,size(xx));
axes(surfpanel.axes.twoD); cla; hold on; grid on;
ticks = linspace(min(Lcc),max(Lcc),5);
set(gca,'XTick',ticks,'YTick',ticks,'xlim',[min(xx(:)) max(xx(:))],'ylim',[min(yy(:)) max(yy(:))]);
surfpanel.surf.twoD = surf(gca,xx,yy,surface);
set(surfpanel.surf.twoD(1),'edgecolor','none')
alpha(.3);
contour3(xx,yy,surface,'linewidth',2);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(GLMP.subunit{sub}.meannspikes);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mn,'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','white')
end
xlabel('Lcc');
ylabel('Mcc');
zlabel('# of spikes')
colormap('cool')

% Replot fitspanel
axes(fitspanel.axes); cla; hold on; grid on;
p1 = plot(fitspanel.angs,fitspanel.LL.oneD.(str1d),'o-','color',conpanel.oneDcol);
p2 = plot(fitspanel.angs,fitspanel.LL.twoD.(str2d),'o-','color',conpanel.twoDcol);
plot(fitspanel.angs(surfidx),fitspanel.LL.oneD.(str1d)(surfidx),'*','color',conpanel.oneDcol)
plot(fitspanel.angs(surfidx),fitspanel.LL.twoD.(str2d)(surfidx),'*','color',conpanel.twoDcol)
bound = max(fitspanel.LL.oneD.(str1d)) - chi2inv(.95,1)/2;
plot([-pi pi],[bound bound],'--','color',conpanel.oneDcol)
bound = max(fitspanel.LL.twoD.(str2d)) - chi2inv(.95,1)/2;
plot([-pi pi],[bound bound],'--','color',conpanel.twoDcol)
set(p1,'ButtonDownFcn',@surfselect)
set(p2,'ButtonDownFcn',@surfselect)

% Save user data
set(confintFig.conpanel,'userdata',conpanel)
set(confintFig.fitspanel,'userdata',fitspanel)
set(confintFig.surfpanel,'userdata',surfpanel)
set(gcf,'userdata',confintFig)

end


%%% Call other analyses %%%

function LoadOverview(~,~)
global GLMSPopData GLMP DN

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'userdata');

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
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'userdata');

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


%%% Setup %%%

function setupFig()

% Build figure
figure(55589); clf;
set(gcf,'units','normalized','pos',[.1 .1 .8 .8],'numbertitle','off',...
    'name','Population Confidence Intervals')

% Set up panels
confintFig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .24 .98]);
confintFig.poppanel = uipanel(gcf,'units','normalized','pos',[.26 .01 .24 .98]);
confintFig.fitspanel = uipanel(gcf,'units','normalized','pos',[.51 .41 .48 .58]);
confintFig.surfpanel = uipanel(gcf,'units','normalized','pos',[.51 .01 .48 .39]);

%%% Control Panel %%%
conpanel.startanalysis = uicontrol('style','pushbutton',...
    'parent',confintFig.conpanel,'units','normalized','pos',[.05 .01 .4 .04],...
    'string','Reanalyze Cell','fontsize',12,'callback',@analyzeCell,...
    'backgroundcolor',[.2 1 0]);
conpanel.uicontrols.reanalAll = uicontrol('parent',confintFig.conpanel,...
    'style','pushbutton','string','Reanalyze All','units','normalized',...
    'pos',[.55 .01 .4 .04],'fontsize',12,'callback',@ReanalAll,...
    'backgroundcolor',[1 .2 0]);
conpanel.uicontrols.overview = uicontrol('parent',confintFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .06 .4 .04],'fontsize',12,'callback',@LoadOverview);
conpanel.uicontrols.analGUI = uicontrol(confintFig.conpanel,...
    'style','pushbutton','string','Load Control Panel','units','normalized',...
    'pos',[.55 .06 .4 .04],'fontsize',12,'callback',@LoadControlPanel);
conpanel.axes.scatterplot = axes('parent',confintFig.conpanel,...
    'pos',[.15 .25 .75 .3]); box on; grid on;
conpanel.normLLthresh = .075;
conpanel.labels.normLLthresh = uicontrol(confintFig.conpanel,'units','normalized',...
    'style','text','pos',[.05 .11 .4 .03],'string','Norm LL Thresh','fontsize',12);
conpanel.uicontrols.normLLthresh = uicontrol(confintFig.conpanel,'units','normalized',...
    'style','edit','pos',[.55 .11 .35 .03],'string',conpanel.normLLthresh,...
    'Callback',@cellselect);
conpanel.selectedidx = [];
conpanel.oneDcol = [1 .5 0];
conpanel.twoDcol = [0 0 1];
conpanel.oneDL = [];

%%% population panel %%%
poppanel.axes.oneD = axes('parent',confintFig.poppanel,'unit','normalized',...
    'pos',[.15 .575 .75 .4]);
poppanel.axes.twoD = axes('parent',confintFig.poppanel,'units','normalized',...
    'pos',[.15 .15 .75 .4]);
strs = {'All 1D Surfs' 'Unichromatic' 'Bichromatic'};
poppanel.uicontrols.which1Dsurf = uicontrol(confintFig.poppanel,...
    'style','popup','units','normalized','pos',[.1 .55 .8 .02],...
    'string',strs,'value',1,'callback',@resetAll);
strs = {'All 2D Surfs' 'Unichromatic hyperbola' 'Bichromatic hyperbola' 'Unichromatic elipse' 'Bichromatic elipse'};
poppanel.uicontrols.which2Dsurf = uicontrol(confintFig.poppanel,...
    'style','popup','units','normalized','pos',[.1 .125 .8 .02],...
    'string',strs,'value',1,'callback',@resetAll);
strs = {'Both Monkeys' 'Nut Data' 'Maui Data'};
poppanel.uicontrols.whichMonk = uicontrol('parent',confintFig.poppanel,...
    'style','popup','units','normalized','pos',[.1 .05 .8 .05],...
    'string',strs,'callback',@resetAll);
poppanel.nbins = 20;
poppanel.oneDL = [];
poppanel.twoDL = [];

%%% Fitspanel %%%
fitspanel.axes = axes('parent',confintFig.fitspanel,'units','normalized',...
    'pos',[.075 .15 .9 .8],'xlim',[-pi pi],'xtick',linspace(-pi,pi,5),...
    'xticklabel',linspace(-180,180,5)); box on; grid on; hold on;
xlabel('Fixed Preferred Direction')
ylabel('Likelihood')
fitspanel.options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',10.^-9,'UseParallel',1);
fitspanel.surftype = 'conicsection_xy';
fitspanel.errortype = 'NegativeBinomial';
fitspanel.nangs = 361; % number of bins from -pi to pi inclusive
fitspanel.angs = linspace(-pi,pi,fitspanel.nangs);

% Surfpanel
surfpanel.axes.oneD = axes('parent',confintFig.surfpanel,'units','normalized',...
    'pos',[.1 .1 .35 .8]); box on; grid on; hold on;
title('1D Fit')
surfpanel.axes.twoD = axes('parent',confintFig.surfpanel,'units','normalized',...
    'pos',[.6 .1 .35 .8]); box on; grid on; hold on;
title('2D Fit')

% Save user data
set(confintFig.conpanel,'userdata',conpanel)
set(confintFig.poppanel,'userdata',poppanel)
set(confintFig.fitspanel,'userdata',fitspanel)
set(confintFig.surfpanel,'userdata',surfpanel)
set(gcf,'userdata',confintFig)

end

function PlotTuningHists()

% Load figure variables
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'userdata');
poppanel = get(confintFig.poppanel,'userdata');
%fitspanel = get(confintFig.fitspanel,'userdata');
%surfpanel = get(confintFig.surfpanel,'userdata');

% Which monkey
if poppanel.uicontrols.whichMonk.Value == 1
    monkL = ones(size(poppanel.monk.NutL));
elseif poppanel.uicontrols.whichMonk.Value == 2
    monkL = poppanel.monk.NutL;
elseif poppanel.uicontrols.whichMonk.Value == 3
    monkL = poppanel.monk.MauiL;
end

% For plotting 1 subunit per cell
%datafiles = conpanel.table.Data(:,1);
datafiles = poppanel.filenames;
%[uniquefile,a,b] = unique(datafiles); uniqL = zeros(size(datafiles)); uniqL(a) = 1; %for 1 sub/cell
uniqL = ones(size(datafiles)); % all subunits

% Which 1D surface
if poppanel.uicontrols.which1Dsurf.Value == 1
    surf1dL = poppanel.oneD.L;
elseif poppanel.uicontrols.which1Dsurf.Value == 2
    surf1dL = poppanel.oneD.unichrom.L;
elseif poppanel.uicontrols.which1Dsurf.Value == 3
    surf1dL = poppanel.oneD.bichrom.L;
end

% Which 2D Surface
if poppanel.uicontrols.which2Dsurf.Value == 1
    surf2dL = ones(size(poppanel.oneD.L));
elseif poppanel.uicontrols.which2Dsurf.Value == 2
    surf2dL = poppanel.twoD.unichrom.hypL;
elseif poppanel.uicontrols.which2Dsurf.Value == 3
    surf2dL = poppanel.twoD.bichrom.hypL;
elseif poppanel.uicontrols.which2Dsurf.Value == 4
    surf2dL = poppanel.twoD.unichrom.eliL;
elseif poppanel.uicontrols.which2Dsurf.Value == 5
    surf2dL = poppanel.twoD.bichrom.eliL;
end

% Which 1d and 2d surfs to show
poppanel.oneDL = poppanel.oneD.L & monkL & surf1dL & uniqL & ~poppanel.excludeL;
poppanel.twoDL = ~poppanel.oneD.L & monkL & surf2dL & uniqL & ~poppanel.excludeL;

% Bin centers, edges, and counts
bins = linspace(-pi,pi,poppanel.nbins*2+1);
poppanel.hist.binedges = bins(1:2:end);
poppanel.hist.bincenters = bins(2:2:end);
poppanel.hist.counts.oneD = histcounts(poppanel.tuning(poppanel.oneDL),poppanel.hist.binedges);
poppanel.hist.counts.twoD = histcounts(poppanel.tuning(poppanel.twoDL),poppanel.hist.binedges);
%radpm = str2double(conpanel.uicontrols.wedgesize.String)/2/180*pi;

% Plot preferred axes of 1D cells
angs = poppanel.tuning(poppanel.oneDL);
axes(poppanel.axes.oneD); cla;
% for n = 1:4
%     h = polar([0 conpanel.cardirs(n)+radpm conpanel.cardirs(n)-radpm 0],...
%         [0 max(poppanel.hist.counts.oneD)*1.1 max(poppanel.hist.counts.oneD)*1.1 0]);
%     set(h,'color',conpanel.carcols(n,:)); hold on;
% end
poppanel.hist.oneD = rose(angs,poppanel.nbins); hold on;
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13)); % Deleting angle labels and outermost rho label
end
axis tight
set(poppanel.hist.oneD,'color',conpanel.oneDcol,'ButtonDownFcn',@SelectTuningDirs)
title('1D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')
t = get(gca,'xlim');
text(t(2)*.6,t(2)*.6,['n = ' num2str(numel(angs))])
%angs = conpanel.tuning(conpanel.oneD.L & monkL & surf1dL & uniqL & ~conpanel.exclude & conpanel.oneD.selectedangs);
%if sum(conpanel.oneD.selectedangs > 0)
%    poppanel.hist.oneDselected = rose2(angs,poppanel.nangs,conpanel.oneDcol,...
%        'edgecolor',conpanel.oneDcol,'ButtonDownFcn',@SelectTuningDirs);
%end

% Plot axes of 2D cells
angs = poppanel.tuning(poppanel.twoDL);
axes(poppanel.axes.twoD); cla;
% for n = 1:4
%     h = polar([0 conpanel.cardirs(n)+radpm conpanel.cardirs(n)-radpm 0],...
%         [0 max(poppanel.hist.counts.twoD)*1.1 max(poppanel.hist.counts.twoD)*1.1 0]);
%     set(h,'color',conpanel.carcols(n,:)); hold on;
% end
poppanel.hist.twoD = rose(angs,poppanel.nbins); hold on; axis tight
t = findall(gca,'type','text');
if size(t,1) > 13
    delete(t(1:13));% Deleting angle labels and outermost rho label
end
set(poppanel.hist.twoD,'color',conpanel.twoDcol,'ButtonDownFcn',@SelectTuningDirs)
title('2D Tuning')
xlabel('L-Cone Contrast')
ylabel('M-Cone Contrast')
t = get(gca,'xlim');
text(t(2)*.6,t(2)*.6,['n = ' num2str(numel(angs))])
% angs = conpanel.tuning(~conpanel.oneD.L & monkL & surf2dL & uniqL & ~conpanel.exclude & conpanel.twoD.selectedangs);
% if sum(conpanel.twoD.selectedangs > 0)
%     poppanel.hist.twoDselected = rose2(angs,poppanel.nangs,conpanel.twoDcol);
%     set(poppanel.hist.twoDselected,'ButtonDownFcn',@SelectTuningDirs);
% end


% Save Data
set(confintFig.poppanel,'UserData',poppanel)
set(confintFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',confintFig)

end

function UnpackPopData()
global GLMSPopData

% Load user data
confintFig = get(gcf,'userdata');
conpanel = get(confintFig.conpanel,'UserData');
poppanel = get(confintFig.poppanel,'userdata');
fitspanel = get(confintFig.fitspanel,'userdata');

% Grab saved population data
if ismac
    conpanel.library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    conpanel.library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    conpanel.library = 'C:\Users\Patty\Dropbox\Patrick\GLMS Data\';
end
load([conpanel.library 'GLMSPopData.mat'])

% Pull out some tuning data from population
datatypes = GLMSPopData(1,:);

% Datafiles
poppanel.filenames = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));

% Which monkey
monks = cat(1,poppanel.filenames{:});
poppanel.monk.NutL = monks(:,1) == 'N';
poppanel.monk.MauiL = monks(:,1) == 'M';
conpanel.monkL = true(size(poppanel.monk.NutL));

% Surf parameters
surfparams = [GLMSPopData{2:end,strcmp(datatypes,'Surface Parameters')}]';
poppanel.surftype = surfparams(1).oneD.surftype;
fitcomps = [surfparams.fitcomps]';
oneD = [fitcomps.oneD]';
parvals = cat(1,oneD.parvals);
poppanel.oneD.unichrom.parvals = parvals(1:2:end,:);
poppanel.oneD.bichrom.parvals = parvals(2:2:end,:);
nLLs = [oneD.nLLs]';
poppanel.oneD.unichrom.LLs = -nLLs(:,1);
poppanel.oneD.bichrom.LLs = -nLLs(:,2);
normLLs = cat(1,oneD.normLL);
poppanel.oneD.unichrom.normLLs = normLLs(1:2:end);
poppanel.oneD.bichrom.normLLs = normLLs(2:2:end);

twoD = [fitcomps.twoD]';
parvals = cat(1,twoD.parvals);
poppanel.twoD.unichrom.parvals = parvals(1:2:end,:);
poppanel.twoD.bichrom.parvals = parvals(2:2:end,:);
nLLs = [twoD.nLLs]';
poppanel.twoD.unichrom.LLs = -nLLs(:,1);
poppanel.twoD.bichrom.LLs = -nLLs(:,2);
normLLs = cat(1,twoD.normLL);
poppanel.twoD.unichrom.normLLs = normLLs(1:2:end);
poppanel.twoD.bichrom.normLLs = normLLs(2:2:end);

% Distinguish 1D from 2D using bichromatic fits
poppanel.oneD.L = (poppanel.twoD.bichrom.normLLs - poppanel.oneD.bichrom.normLLs) < conpanel.normLLthresh;

% Distinguish unichrom from bichrom
poppanel.oneD.unichrom.L = (poppanel.oneD.bichrom.normLLs - poppanel.oneD.unichrom.normLLs) < conpanel.normLLthresh;
poppanel.oneD.bichrom.L = (poppanel.oneD.bichrom.normLLs - poppanel.oneD.unichrom.normLLs) >= conpanel.normLLthresh;
poppanel.twoD.unichrom.L = (poppanel.twoD.bichrom.normLLs - poppanel.twoD.unichrom.normLLs) < conpanel.normLLthresh;
poppanel.twoD.bichrom.L = (poppanel.twoD.bichrom.normLLs - poppanel.twoD.unichrom.normLLs) >= conpanel.normLLthresh;

% Unichrom and bichrom parvals
poppanel.oneD.params = nan(size(poppanel.oneD.unichrom.parvals));
poppanel.oneD.params(poppanel.oneD.unichrom.L,:) = poppanel.oneD.unichrom.parvals(poppanel.oneD.unichrom.L,:);
poppanel.oneD.params(poppanel.oneD.bichrom.L,:) = poppanel.oneD.bichrom.parvals(poppanel.oneD.bichrom.L,:);
poppanel.twoD.params = nan(size(poppanel.oneD.unichrom.parvals));
poppanel.twoD.params(poppanel.twoD.unichrom.L,:) = poppanel.twoD.unichrom.parvals(poppanel.twoD.unichrom.L,:);
poppanel.twoD.params(poppanel.twoD.bichrom.L,:) = poppanel.twoD.bichrom.parvals(poppanel.twoD.bichrom.L,:);

% All parvals
poppanel.params = nan(size(poppanel.oneD.unichrom.parvals));
poppanel.params(poppanel.oneD.L,:) = poppanel.oneD.params(poppanel.oneD.L,:);
poppanel.params(~poppanel.oneD.L,:) = poppanel.twoD.params(~poppanel.oneD.L,:);

% Sort by surface type
orax = poppanel.twoD.params(:,4);
poppanel.twoD.unichrom.hypL = orax < 0 & poppanel.twoD.unichrom.L & ~poppanel.oneD.L;
poppanel.twoD.unichrom.eliL = orax > 0 & poppanel.twoD.unichrom.L & ~poppanel.oneD.L;
poppanel.twoD.bichrom.hypL = orax < 0 & ~poppanel.twoD.unichrom.L & ~poppanel.oneD.L;
poppanel.twoD.bichrom.eliL = orax > 0 & ~poppanel.twoD.unichrom.L & ~poppanel.oneD.L;

% Tuning
poppanel.tuning = nan(size(poppanel.oneD.L));
poppanel.tuning(poppanel.oneD.L) = poppanel.oneD.params(poppanel.oneD.L,end-1);
poppanel.tuning(~poppanel.oneD.L) = poppanel.twoD.params(~poppanel.oneD.L,end-1);
poppanel.diffnormLL = poppanel.twoD.bichrom.normLLs - poppanel.oneD.bichrom.normLLs;

% Load previously recorded confidence interval data
temp = [GLMSPopData{2:end,strcmp(datatypes,'Confidence Intervals')}]';
tempbf = [temp.bfconfdeg]';
temphess = [temp.hessconfdeg]';
tempLL = [temp.LL]';

% LL
oneD = [tempLL.oneD];
poppanel.oneD.unichrom.LL = [oneD.uc]';
poppanel.oneD.bichrom.LL = [oneD.bc]';
twoD = [tempLL.twoD];
poppanel.twoD.unichrom.LL = [twoD.uc]';
poppanel.twoD.bichrom.LL = [twoD.bc]';

% unichromativ vs bichromatic LL
poppanel.oneD.LL = nan(numel(poppanel.oneD.L),fitspanel.nangs);
poppanel.twoD.LL = nan(numel(poppanel.oneD.L),fitspanel.nangs);
poppanel.oneD.LL(poppanel.oneD.unichrom.L,:) = poppanel.oneD.unichrom.LL(poppanel.oneD.unichrom.L,:);
poppanel.oneD.LL(poppanel.oneD.bichrom.L,:) = poppanel.oneD.bichrom.LL(poppanel.oneD.bichrom.L,:);
poppanel.twoD.LL(poppanel.twoD.unichrom.L,:) = poppanel.twoD.unichrom.LL(poppanel.twoD.unichrom.L,:);
poppanel.twoD.LL(poppanel.twoD.bichrom.L,:) = poppanel.twoD.bichrom.LL(poppanel.twoD.bichrom.L,:);

% 1d vs 2d LL
poppanel.LL = nan(numel(poppanel.oneD.L),fitspanel.nangs);
poppanel.LL(poppanel.oneD.L,:) = poppanel.oneD.LL(poppanel.oneD.L,:);
poppanel.LL(~poppanel.oneD.L,:) = poppanel.oneD.LL(~poppanel.oneD.L,:);

% BF conf
oneD = [tempbf.oneD];
poppanel.oneD.unichrom.bfconfs = [oneD.uc]';
poppanel.oneD.bichrom.bfconfs = [oneD.bc]';
twoD = [tempbf.twoD];
poppanel.twoD.unichrom.bfconfs = [twoD.uc]';
poppanel.twoD.bichrom.bfconfs = [twoD.bc]';

% unichromativ vs bichromatic conf
poppanel.oneD.bfconfs = nan(size(poppanel.oneD.L));
poppanel.twoD.bfconfs = nan(size(poppanel.oneD.L));
poppanel.oneD.bfconfs(poppanel.oneD.unichrom.L) = poppanel.oneD.unichrom.bfconfs(poppanel.oneD.unichrom.L);
poppanel.oneD.bfconfs(poppanel.oneD.bichrom.L) = poppanel.oneD.bichrom.bfconfs(poppanel.oneD.bichrom.L);
poppanel.twoD.bfconfs(poppanel.twoD.unichrom.L) = poppanel.twoD.unichrom.bfconfs(poppanel.twoD.unichrom.L);
poppanel.twoD.bfconfs(poppanel.twoD.bichrom.L) = poppanel.twoD.bichrom.bfconfs(poppanel.twoD.bichrom.L);

% 1d vs 2d conf
poppanel.bfconfs = nan(size(poppanel.oneD.L));
poppanel.bfconfs(poppanel.oneD.L) = poppanel.oneD.bfconfs(poppanel.oneD.L);
poppanel.bfconfs(~poppanel.oneD.L) = poppanel.oneD.bfconfs(~poppanel.oneD.L);

% Hess conf
oneD = [temphess.oneD];
poppanel.oneD.unichrom.hessconfs = [oneD.uc]';
poppanel.oneD.bichrom.hessconfs = [oneD.bc]';
twoD = [temphess.twoD];
poppanel.twoD.unichrom.hessconfs = [twoD.uc]';
poppanel.twoD.bichrom.hessconfs = [twoD.bc]';

% unichromativ vs bichromatic conf
poppanel.oneD.hessconfs = nan(size(poppanel.oneD.L));
poppanel.twoD.hessconfs = nan(size(poppanel.oneD.L));
poppanel.oneD.hessconfs(poppanel.oneD.unichrom.L) = poppanel.oneD.unichrom.hessconfs(poppanel.oneD.unichrom.L);
poppanel.oneD.hessconfs(poppanel.oneD.bichrom.L) = poppanel.oneD.bichrom.hessconfs(poppanel.oneD.bichrom.L);
poppanel.twoD.hessconfs(poppanel.twoD.unichrom.L) = poppanel.twoD.unichrom.hessconfs(poppanel.twoD.unichrom.L);
poppanel.twoD.hessconfs(poppanel.twoD.bichrom.L) = poppanel.twoD.bichrom.hessconfs(poppanel.twoD.bichrom.L);

% 1d vs 2d conf
poppanel.hessconfs = nan(size(poppanel.oneD.L));
poppanel.hessconfs(poppanel.oneD.L) = poppanel.oneD.hessconfs(poppanel.oneD.L);
poppanel.hessconfs(~poppanel.oneD.L) = poppanel.oneD.hessconfs(~poppanel.oneD.L);

% Pull out data for table
data = cell(size(GLMSPopData,1)-1,3);
data(:,1) = GLMSPopData(2:end,strcmp(datatypes,'Datafile'));
data(:,2) = GLMSPopData(2:end,strcmp(datatypes,'Subunit'));
data(:,3) = GLMSPopData(2:end,strcmp(datatypes,'normLLdiff'));

% Repackage data and name columns
colname = {'Datafile' 'Sub' 'Norm LL Diff'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};

% Display table
conpanel.table = uitable('parent',confintFig.conpanel,...
    'units','normalized','pos',[.01 .6  .98 .39],...
    'data',data,'columnname',colname,...
    'columnformat',colformat,...
    'BackgroundColor',[1 1 1],...
    'CellSelectionCallback',@cellselect);

% Save user data
set(confintFig.conpanel,'userdata',conpanel)
set(confintFig.poppanel,'userdata',poppanel)
set(gcf,'userdata',confintFig)

end

function exclusioncriteria()
global GLMSPopData

% Load figure variables
confintFig = get(gcf,'userdata');
poppanel = get(confintFig.poppanel,'userdata');

% parmeter bounds
% modmin = str2double(get(conpanel.uicontrols.modmin,'string')); % in sp/s
% kappamax = str2double(get(conpanel.uicontrols.kappamax,'string'));
% confmax = str2double(get(conpanel.uicontrols.confmax,'string'));
modmin = 10; % in sp/s
kappamax = 4;
confmax = 360;

% Load data
datatypes = GLMSPopData(1,:);

% Modulation
modvals = nan(size(GLMSPopData,1)-1,1);
for n = 2:size(GLMSPopData,1)
    GLMP = GLMSPopData{n,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{n,strcmp(datatypes,'Subunit')};
    modvals(n-1) = max(GLMP.subunit{sub}.meannspikes)...
        ./ mean(GLMP.subunit{sub}.stimDur) - mean(GLMP.subunit{sub}.blfr);
end

% Kappa
kappavals = nan(size(GLMSPopData,1)-1,1);
kappavals(poppanel.oneD.L) = poppanel.oneD.params(poppanel.oneD.L,end);
kappavals(~poppanel.oneD.L) = poppanel.twoD.params(~poppanel.oneD.L,end);

% Tuning Std
%confvals = nan(size(GLMSPopData,1)-1,1);
%confint = [GLMSPopData{2:end,strcmp(datatypes,'Confidence Intervals')}]';
%t = [confint.bfconfdeg]';
%onedconf = [t.oneD]';
%confvals(poppanel.oneD.L) = onedconf(poppanel.oneD.L); % including just 1D confs
confvals = poppanel.bfconfs;

% Which don't meet criteria?
poppanel.excludeL = modvals < modmin | kappavals > kappamax | confvals > confmax;
disp([num2str(sum(poppanel.excludeL)) ' datafiles excluded based on current criteria.'])

% Save user data
set(confintFig.poppanel,'userdata',poppanel)
set(gcf,'userdata',confintFig)

end

