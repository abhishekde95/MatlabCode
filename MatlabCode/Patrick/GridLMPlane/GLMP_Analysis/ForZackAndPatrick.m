

%%
% A couple of analyses for Patrick's experiment
% switching over the B-splines
stro = nex2stro('Users/jpatrickweller/Documents/MATLAB/GridLMPlane/Datafiles/S111011003.nex');% Symmetric L-M Trench
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Lcc'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Mcc'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
ntrials = size(stro.trial,1);
LM = Lcc+Mcc;
LvM = 10*(Lcc-Mcc);

bins = linspace(0,2,50);
PSTH = zeros(1,length(bins));
% First PSTHs to find good values for "offset"
for i = 1:ntrials
    spiketimes = stro.ras{i,1}-stimon_t(i);
    PSTH = PSTH +histc(spiketimes',bins);
end

figure; axes; hold on;
bar(bins,PSTH);
gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
bguess = mean(PSTH(bins<.25));
aguess = max(PSTH)-bguess;
muguess = bins(find(PSTH == max(PSTH),1,'first'));
sigmaguess = .1; % terrible. Need something better.
%plot(bins,gauss(bins,[aguess,bguess,muguess,sigmaguess]),'m.')
fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3)
offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);

% Now counting up spikes in a window
ntrials = size(stro.trial,1);
nspikes = zeros(ntrials,1);
for i = 1:ntrials
    nspikes(i) = sum(stro.ras{i,1} > stimon_t(i)+offset(1) & stro.ras{i,1}<stimon_t(i)+offset(2));
end
%%
% Projecting onto a few unit vectors in the LM plane and 
% Fitting the 1-D function (response vs. projection) with 
% regular old 1-D cubic smoothing spline.
nthetas = 1000;
thetas = linspace(0,pi,nthetas);
data = nan*ones(nthetas,1);
figure; axes; hold on;
smoothingparam = .9;
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    projs = [LM LvM]*unitvect;
    lastwarn('');
    pp  = csaps(projs,nspikes,smoothingparam);
    plot(projs,nspikes,'k.');
    fnplt(pp);
    drawnow;
    cla
    data(i) = sum((fnval(pp,projs)-nspikes).^2);
    if ~isempty(lastwarn)
        data(i) = nan;
    end
end
figure; plot(thetas,data);
preftheta = thetas(data == min(data))
% Taking a look at the residuals from a 1-D fit
unitvect = [cos(preftheta) sin(preftheta)]';
projs = [LM LvM]*unitvect;
pp  = csaps(projs,nspikes,smoothingparam);
predfr1D = fnval(pp,projs);
figure; 
plot(predfr1D,nspikes-predfr1D,'k.');
xlabel('fitted value'); ylabel('residual');
[x,y]= meshgrid(linspace(-1,1,50),linspace(-1,1,50));
projs = [x(:) y(:)]*unitvect;
predfr2D = fnval(pp,projs);
figure; axes; hold on; axis square;
surf(x,y,reshape(predfr2D,size(x)));
plot3(LM,LvM,nspikes,'k.');
LMconeweights = 1/sqrt(2)*[1 1;1 -1]*unitvect;
title(['L weight = ',num2str(LMconeweights(1),1),' M weight = ',num2str(LMconeweights(2),1)]);
%%
% Trying a thinplate spline on the sqrt-transformed spikecounts.

figure; axes; hold on;
out2d = tpaps([LM, LvM]', nspikes');
plot3(LM,LvM,nspikes,'k.');
h = fnplt(out2d);
surf(h{1},h{2},h{3})

%out2d = tpaps([LM, LvM]', sqrt(nspikes'));
%plot3(LM,LvM,sqrt(nspikes),'k.');
%h = fnplt(out2d);
%surf(h{1},h{2},h{3}.^2)
axis square;

%out2d_noxform = tpaps([LM, LvM]', nspikes');
%h = fnplt(out2d_noxform);
%h =surf(h{1},h{2},h{3},repmat(5,size(h{1})))
%set(h,'FaceAlpha',.7)

% The square root transformation looks like it's pretty subtle.

% Now some 1-D splines
%outx = tpaps([LM, repmat(mean(LvM),size(LvM,1),1)]', sqrt(nspikes)',1);
%outx = csaps(LM, sqrt(nspikes));
%outy = csaps(LvM, sqrt(nspikes));
%figure;
%subplot(2,1,1);
%fnplt(outx);
%subplot(2,1,2);
%fnplt(outy);
%%
% Bootstrapping test to determine whether the 2-D fit is actually helping at all.
nboot = 2000;
PARAMETRIC = 2; % 0=Nonparametric, 1=Poisson, 2=Gaussian
unitvect = [cos(preftheta) sin(preftheta)]';
projs = [LM LvM]*unitvect;
%projs = [LM LvM]*[unitvect(2) unitvect(1)]';
pp1D  = csaps(projs,nspikes,smoothingparam);
predfr1D = fnval(pp1D,projs);
resid = predfr1D-nspikes;
SS1D = sum(resid.^2);
SS2D = sum((fnval(out2d,[LM LvM]')-nspikes').^2);
[u,a,b] = unique([LM LvM],'rows');
SS = zeros(nboot,2);
for i = 1:nboot
    i
    tmpdata = zeros(size(u,1),1);
    for j = 1:max(b)    
        L = b == j;
        r = resid(L);
        if (PARAMETRIC == 0)
            tmpdata(L) = predfr1D(L) + r(unidrnd(length(r),length(r),1)); % Nonparametric
        elseif (PARAMETRIC == 1)
            tmpdata(L) = poissrnd(predfr1D(L)); % Parametric (Poisson)
        else
            tmpdata(L) = predfr1D(L) + normrnd(mean(r),std(r),size(r)); % Parametric (Gaussian)
        end
    end
    pp  = csaps(projs,tmpdata,smoothingparam);
    predfr = fnval(pp,projs);
    SS(i,1) = sum((predfr-tmpdata).^2);
    out = tpaps([LM, LvM]', tmpdata');
    SS(i,2) = sum((fnval(out,[LM LvM]')-tmpdata').^2);
end
SSratio = SS(:,1)./SS(:,2);  % Large values mean that TP spline fits way better
figure; axes; hold on;
hist(SSratio);
plot(SS1D./SS2D,0,'m*');
p = sum((SS1D./SS2D)<SSratio)/nboot;
title(['p = ',num2str(p)]);

%%
% Fitting individual contrast-response functions by least-squares
% This only works for cells that were tested on a radial grid 
thetas = atan2(LvM, LM);
bins = linspace(-pi,pi,1000);
n = histc(thetas,bins);
uniquethetas = bins(n>(max(n)/2));
fittedparams = [];
figure; axes; hold on;
options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
vlb = [0  0.001 0.001 0];
vub = [1000 1000 100 100];
for i = 1:length(uniquethetas)-1
    unitvect = [cos(uniquethetas(i)) sin(uniquethetas(i))]';
    projs = [LM LvM]*unitvect;
    L = thetas > uniquethetas(i) & thetas < uniquethetas(i+1);
    plot(projs(L),nspikes(L),'k.')
    % Using the thin-plate spline to derive good initial guesses
    params0 = [max(nspikes(L)), prctile(abs(projs(L)),1), 1, 0];
    params1 = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,abs(projs(L)),fnval(out2d,[LM(L) LvM(L)]')'); 
    f = fmincon('MyFitNakaRushtonFun',params1,[],[],[],[],vlb,vub,[],options,abs(projs(L)),nspikes(L)); 
    
    fittedparams(i,:) = f;
    x = linspace(0,max(projs(L)),100);
    y = MyComputeNakaRushton(f,x);
    plot(x,y,'b-');
    title(num2str(uniquethetas(i)));
    drawnow;    
    pause(.5);
    cla;
end

figure;
for i = 1:size(fittedparams,2)
    subplot(3,1,i);
    polar(uniquethetas(1:end-1),fittedparams(:,i)','k.-')
end

%%
% Now fitting a single naka-rushton to the full data set, assuming L+M tuning
options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
vlb = [0  0.001 0.001 0];
vub = [1000 1000 100 100];

params0 = fittedparams(find(uniquethetas<0,1,'last'),:); % Using the L+M spoke as a first guess

f = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,abs(LM),nspikes);
% Note that in the above line we're fitting to the absolute values of LM

% Plotting the fit
figure; axes; hold on;
plot3(LM,LvM,nspikes,'k.');
x = linspace(min(LM),max(LM),100);
z = MyComputeNakaRushton(f,abs(x));
h = surf(repmat(x,length(x),1),repmat(x,length(x),1)',repmat(z,length(z),1))
set(h,'EdgeAlpha',0);


%%
% Fitting a Naka-Rushton to the full data set, allowing different
% A, sigma, and n for different luminance polarities
params0 = [f(1) f(1) f(2) f(2) f(3) f(3) f(4)];
vlb = [0  0 0.001 0.001 0.001 0.001 0 0];
vub = [1000 1000 1000 1000 100 100 100 100];

ff = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,LM,nspikes);
% Plotting the fit
figure; axes; hold on;
plot3(LM,LvM,nspikes,'k.');
x = linspace(min(LM),max(LM),100);
z = MyComputeNakaRushton(ff,x);
h = surf(repmat(x,length(x),1),repmat(x,length(x),1)',repmat(z,length(z),1))
set(h,'EdgeAlpha',0);



%%
% ---------------------------------------
% Everything below this point is obsolete
% ---------------------------------------
%
%%
% Projecting onto a few unit vectors in the LM plane and 
% Fitting the 1-D function (response vs. projection) with a sum
% of weighted Gaussians.
nthetas = 100;
thetas = linspace(0,pi,nthetas);
data = nan*ones(nthetas,1);
nknots = 20;  sigmascale = nknots/70;
X = zeros(length(LM),1);
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    projs = [LM LvM]*unitvect;
 %   knots = linspace(min(projs),max(projs),nknots);
    knots = linspace(.5*min(projs),2*max(projs),nknots);
    sigma = range(projs)*sigmascale;
    X = exp((-(repmat(projs,1,size(knots,2))-repmat(knots,size(projs,1),1)).^2)./sigma.^2);
    %for j = 1:nknots-4
    %    sp = spmak(knots(j:j+4),1);
    %    X(:,j) = spval(sp,projs); 
    %end
    [b,dev,stats] = glmfit(X,nspikes,'poisson','link','identity','constant','off');
    data(i) = dev;
end
figure; subplot(2,1,1); hold on;
plot(thetas,data,'.-')
set(gca,'YScale','log')
subplot(2,1,2); hold on;
unitvect = [cos(thetas(data == min(data))) sin(thetas(data == min(data)))]';
projs = [LM LvM]*unitvect;
sigma = range(projs)*sigmascale;
%knots = linspace(min(projs),max(projs),nknots);
knots = linspace(.5*min(projs),2*max(projs),nknots);
%clear y
%for j = 1:nknots-4
%    sp = spmak(knots(j:j+4),1);
%    y(:,j) = spval(sp,projs);
%end
X = exp((-(repmat(projs,1,size(knots,2))-repmat(knots,size(projs,1),1)).^2)./sigma.^2);
[b,dev,stats] = glmfit(X,nspikes,'poisson','link','identity','constant','off');
x = linspace(min(projs),max(projs),200)';
plot(projs,nspikes,'k.');
%clear y
%for j = 1:nknots-4
%    sp = spmak(knots(j:j+4),1);
%    y(:,j) = spval(sp,x);
%end
y = exp((-(repmat(x,1,size(knots,2))-repmat(knots,size(x,1),1)).^2)./sigma.^2);
plot(x,y*b,'y-','LineWidth',3);
title(num2str(unitvect));

% Now seeing whether there is any structure at all in the orthogonal axis
orthvect = [unitvect(2) -unitvect(1)]';
projs = [LM LvM]*orthvect;
figure; axes; hold on;
plot(projs,nspikes,'k.');
[b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
title(['p = ',num2str(stats.p(2))])
plot([min(projs) max(projs)],b(2)*[min(projs) max(projs)]+b(1),'k-')
% Hmm.. still some structure.

% Are there any projections with no linear structure? 
% This analysis is sketchy because many of the fits do not converge
ps = zeros(length(thetas),1);
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    projs = [LM LvM]*unitvect;
    [b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
    ps(i) = stats.p(2);
end
orthvect = [cos(thetas(ps == max(ps))) sin(thetas(ps == max(ps)))]';
projs = [LM LvM]*orthvect;
figure; axes; hold on;
plot(projs,nspikes,'k.');
[b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
title(['p = ',num2str(stats.p(2))])
plot([min(projs) max(projs)],b(2)*[min(projs) max(projs)]+b(1),'k-')

figure; axes; hold on;
plot3(LM, LvM, nspikes,'k.');
plot([0 unitvect(1)],[0 unitvect(2)],'m-','LineWidth',2)
plot([0 orthvect(1)],[0 orthvect(2)],'g-','LineWidth',2)

% Projecting onto the L-M axis
%projs = [LM LvM]*orthvect;
%figure; axes; hold on;
%plot(projs,nspikes,'k.');
%[b,dev,stats] = glmfit(projs,nspikes,'poisson','link','identity');
%stats.p
%plot([min(projs) max(projs)],b(2)*[min(projs) max(projs)]+b(1),'k-')

% How good a job do we do representing the surface with a single direction?

%%
% How about just fitting a glm using L+M and L-M?
X = [LM, LM.^2, LvM, LvM.^2, LM.*LvM];
%X = [Lcc, Lcc.^2, Mcc, Mcc.^2, Lcc.*Mcc];

[b,dev,stats] = glmfit(X,nspikes,'poisson','link','identity','constant','on');
stats.t


%%
% Polynomial regression with a bunch of rotations of the data
% Looking for a projection where including the orthogonal axis
% with linear and quadratic terms does not increase significance.
thetas = linspace(-pi/4,pi/4,100);
data = [];
for i = 1:length(thetas);
    unitvect = [cos(thetas(i)) sin(thetas(i))]';
    orthvect = [unitvect(2) -unitvect(1)]';
    projs1 = [LM LvM]*unitvect;
    projs2 = [LM LvM]*orthvect;
    [b1,dev1,stats1] = glmfit([projs1 projs1.^2],nspikes,'poisson','link','identity');
    [b2,dev2,stats2] = glmfit([projs2 projs2.^2],nspikes,'poisson','link','identity');
    [b3,dev3,stats3] = glmfit([projs1 projs1.^2 projs2 projs2.^2 projs1.*projs2],nspikes,'poisson','link','identity');
    overall_p =  1-chi2cdf((dev1-dev3),3);
    data = [data; b1' b2' b3' overall_p];
end

% Trying to find a projection with little structure
idx = find(data(:,end) == max(data(:,end)));
unitvect = [cos(thetas(idx)) sin(thetas(idx))]'; % proj onto this vector has *least* structure
orthvect = [unitvect(2) -unitvect(1)]';
projs1 = [LM LvM]*orthvect;
projs2 = [LM LvM]*unitvect;
[x,y] = meshgrid(linspace(-1,1,30),linspace(-1,1,30));
x = x(:); y = y(:);
z1 = [ones(size(x,1),1) x x.^2 y y.^2 x.*y]*data(idx,[7:12])';
z2 = [ones(size(x,1),1) x x.^2]*data(idx,[1:3])';
figure; axes; hold on;
surf(reshape(x,30,30),reshape(y,30,30),reshape(z1,30,30));
plot3(projs2,projs1,nspikes,'k.');
title(['p = ',num2str(data(idx,end))]);
axis square;

%% 
% Looking for non-stationarity
figure; axes; hold on;
plot(nspikes,'k.')
[b,dev,stats] = glmfit([1:length(nspikes)],nspikes,'poisson','link','identity');
plot([1 length(nspikes)],[1 length(nspikes)]*b(2)+b(1),'m-','LineWidth',5);
title(['p=',num2str(stats.p(2))]);