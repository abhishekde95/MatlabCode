% LMTF Isosamp figures
% Section 1: The Watson 1986 model
%
% Section 2: A single LMTF data set (at one retinal location) and a fit
% Section 3: A Gabor of the variety used in LMTF
% Section 4: Scatter plot of retinal positions tested
% Section 5: SVD on model parameters for Apollo: xi_lum, xi_rg, and theta
% Section 6: Plotting summed absolute residuals for each model
% Section 7: Superimposed fits on data from a single retinal location
% Section 8: Luminance and chromatic sensitivity across the visual field
% Section 9: A single 1-D luminance contrast senstivity function with some
% nearby data points.
% Section 9.1: Other slices through the model fit with data points.
% Section 10: 1-D contrast sensitivity functions (two for each observer)
% Section 11: Rasters from a single cell tested with IsoSamp near-threshold
% stimuli
% Section 12: Number of subjects tested at each screen location (pie
% charts)
% Section 13: Analysis of residuals from first round fits
% Section 14: Analysis of residuals from mode 5 fits
% Section 15: ANOVA around monkeys and humans low temporal frequency
% luminance senstivity
% Section 16: Data for Table 2 (model parameters)
% Section 17: Reviewer figure showing STA and responses to L+M
% Section 18: GIF for icon
%%
% Section 1) The Watson (1986) 1D temporal contrast sensitivity model
xi = 2;
zeta = 1;
% filter 1
n1 = 1;
tau1 = 10^-1.67;
% filter 2
n2 = n1+1;
tau2 = 10^(-1.4 +.17);
maxT = 0.12 % in secs
% General parameters
nsamps = 1000;
t = linspace(0,maxT,nsamps); % time in s

% First the impulse response functions
h1 = tau1.^(n1-1)*(factorial(round(n1)-1)).^-1 .*(t./tau1).^(n1-1).*exp(-t./tau1);
h2 = tau2.^(n2-1)*(factorial(round(n2)-1)).^-1 .*(t./tau2).^(n2-1).*exp(-t./tau2);
%h1 = tau1.^(n1-1)*gamma(n1).^-1 .*(t./tau1).^(n1-1).*exp(-t./tau1);
%h2 = tau2.^(n2-1)*gamma(n2).^-1 .*(t./tau2).^(n2-1).*exp(-t./tau2);
h = xi*(h1-zeta*h2);

% Plotting temporal impulse response
figure('position',[886 172 268 633]);
; subplot(3,1,1); hold on;
plot(t,h1,'b--');
plot(t,zeta*h2,'r--');
plot(t,h,'k-');
xlabel('Time (s)');
ylabel('Amplitude');
set(gca,'Xlim',[0 maxT])

% Analytical FT
omega = logspace(0,2,100);
f1 = (1i*2*pi*tau1.*omega+1).^-n1;
f2 = (1i*2*pi*tau2.*omega+1).^-n2;
fana = xi*(f1-zeta*f2);
subplot(3,1,2);
plot(omega,abs(fana),'k-');
xlabel('Freq (Hz)');
ylabel('Sensitivity');
set(gca,'Yscale','log','Xscale','log');

subplot(3,1,3);
plot(omega,1./abs(fana),'k-');
xlabel('Freq (Hz)');
ylabel('Threshold');
set(gca,'Yscale','log','Xscale','log');

%%
% Section 2) a 2D LMTF detection threshold surface with fit
% Right now fitting the model using quickLMTFfit. Could hardcode.

PLOTMODELFIT = 1;
MAKEMOVIE = 1;
MAKEICON = 0;
CCPLTLIM = 0.3; % cone contrast plotting limits
moviefilename = 'SpinningDataWFit';
% Apollo's data
data = [    0.0506    0.0506    1.0000         0   60.0000   30.0000
   -0.0234    0.0234    1.0000         0   60.0000   30.0000
    0.0421    0.0421   25.0000         0   60.0000   30.0000
   -0.0955    0.1768   25.0000    1.0000   60.0000   30.0000
    0.0661         0    4.7979         0   60.0000   30.0000
    0.0000    0.0326    4.7979         0   60.0000   30.0000
    0.1286         0   10.9520         0   60.0000   30.0000
    0.0000    0.1101   10.9520         0   60.0000   30.0000
    0.0486         0    2.1018         0   60.0000   30.0000
    0.0000    0.0396    2.1018         0   60.0000   30.0000
    0.0344    0.0344    7.2489         0   60.0000   30.0000
   -0.0453    0.0453    7.2489         0   60.0000   30.0000
    0.0397    0.0397    3.1756         0   60.0000   30.0000
   -0.0305    0.0305    3.1756         0   60.0000   30.0000
    0.1736         0   17.9706    1.0000   60.0000   30.0000
    0.0000    0.1181   17.9706         0   60.0000   30.0000
    0.0623    0.0623   15.2360         0   60.0000   30.0000
   -0.0955    0.1252   15.2360    1.0000   60.0000   30.0000
    0.0000    0.0361    1.2810         0   60.0000   30.0000
    0.0656         0    1.2810         0   60.0000   30.0000
    0.0436    0.0436    2.4791         0   60.0000   30.0000
   -0.0241    0.0241    2.4791         0   60.0000   30.0000];

figure; axes; hold on;
Loog = data(:,4) == 1;
plot3(data(~Loog,1),data(~Loog,2),data(~Loog,3),'ko','MarkerFaceColor','black');
plot3(-data(~Loog,1),-data(~Loog,2),data(~Loog,3),'ko','MarkerFaceColor','black');
plot3(data(Loog,1),data(Loog,2),data(Loog,3),'ro','MarkerFaceColor','red');
plot3(-data(Loog,1),-data(Loog,2),data(Loog,3),'ro','MarkerFaceColor','red');
set(gca,'Zscale','log'); set(gca,'Xlim',[-CCPLTLIM CCPLTLIM],'Ylim',[-CCPLTLIM CCPLTLIM],'Zlim',[1 25]);
set(gca,'View',[-45 20]);
axis square
set(gca,'Zticklabel',[1 10]);

if PLOTMODELFIT
    fpar = quickLMTFmodelfit(data(:,1:4));
    [xx,yy,zz] = meshgrid(linspace(-max(abs(data(:,1))),max(abs(data(:,1))),20),...
    linspace(-max(abs(data(:,2))),max(abs(data(:,2))),20),...
    linspace(min(data(:,3)),max(data(:,3)),20));

    f1 = @(omega)fpar(1)*abs(((1i*2*pi*10^fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*10.^(fpar(5)+fpar(6)).*omega+1).^-(fpar(3)+fpar(4))));
    f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*10^fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*10.^(fpar(5+6)+fpar(6+6)).*omega+1).^-(fpar(3+6)+fpar(4+6))));

    a = abs(f1(zz)).^-1; % chromatic threshold
    b = abs(f2(zz)).^-1; % luminance threshold
    thtmp = atan2(yy,xx)-fpar(13); % clockwise rotation from [L,M] to [a,b]
    rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    
    V = sqrt(xx.^2+yy.^2)-rtmp;
    FV = isosurface(xx,yy,zz,V,0);
    h = patch(FV,'facecolor','green');
    set(h,'Edgecolor','none'); % camlight headlight; lighting phong
    alpha(.7);
end

if (MAKEMOVIE)
    set(gca,'Xlimmode','manual','ylimmode','manual');
    plot3([0 0],[0 0],[1 25],'k-');
	plot3([CCPLTLIM CCPLTLIM -CCPLTLIM  -CCPLTLIM CCPLTLIM],[-CCPLTLIM CCPLTLIM CCPLTLIM -CCPLTLIM -CCPLTLIM],[1 1 1 1 1],'k-');
    set(gcf,'Color','white')
    axis vis3d
    originalview = get(gca,'View');
    deltaviews = linspace(0,360,120);
    deltaviews(end) = [];
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[deltaviews(i) 0]);
        set(gca,'xtick',[-.2 0 .2], 'ytick',[-.2 0 .2]);
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    close(writerObj);
end


% 96x96 GIF for icon
nframes = 30;
rotations = linspace(0,pi,nframes+1);
rotations(end) = [];
im = uint8(zeros(96,96,1,length(rotations)));
if MAKEICON
    figure;
    ax = axes; set(ax,'Units','pixels','position',[0 0 96 96]); hold on;
    set(gcf,'Units','pixels','position',[0 0 96 96]);
    h = patch(FV,'facecolor','green');
    plot3(data(~Loog,1),data(~Loog,2),data(~Loog,3),'ko','MarkerFaceColor','black');
    plot3(-data(~Loog,1),-data(~Loog,2),data(~Loog,3),'ko','MarkerFaceColor','black');
    plot3(data(Loog,1),data(Loog,2),data(Loog,3),'ro','MarkerFaceColor','red');
    plot3(-data(Loog,1),-data(Loog,2),data(Loog,3),'ro','MarkerFaceColor','red');
    set(h,'Edgecolor','none');  camlight headlight; lighting phong;
    set(gca,'Zscale','log','Visible','off')
    alpha(.7);
    axis vis3d
    for i = 1:length(rotations)
        view(rotations(i)*180/pi,12);
        export_fig test.tif -nocrop -append
    end
    im2gif 'test.tif' -nocrop
end


%%
% Section 3
% An animated Gabor such as used in LMTF

MAKEMOVIE = 1;
STILLFRAMES = 0;
filenamestem = 'GaborMovie';
M = [0.0761    0.1524    0.0218;
    0.0275    0.1582    0.0321;
    0.0024    0.0118    0.1220];  % Something standard; I think this is Dell 4
lambda = 1;
nframespercycle = 4; % 48, 4
ncycles = 36/nframespercycle;
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';
gaborccs = [0.1 0.1 0];  % L, M, S deltas
%gaborccs = [0.03 -0.03 0];  % L, M, S deltas
moviefilename = [filenamestem,num2str(round(ncycles))];
gaborlms = bkgndlms'.*(gaborccs+1);

gaborrgb = inv(M)*gaborlms';
phis = linspace(0,2*pi*ncycles,nframespercycle*ncycles+1);
phis(end) = [];
figure('units','inches','position',[1 1 4 3]);
axes('units','inches','position',[0 0 4 3],'Visible','off','color',bkgndrgb)
envelope = [linspace(0,1,round(length(phis)/4)) ones(1,round(length(phis)/2))];
envelope = [envelope, linspace(1,0,length(phis)-length(envelope))];

if (MAKEMOVIE)
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
end
for j = 1:length(phis)
    im = DrawGaborEdge(bkgndrgb, gaborrgb, [0 0 0], pi, lambda, .15, 1, phis(j), 0, 0, 0, 0, .999, 30);
    image((im-bkgndrgb(1))*envelope(j)+bkgndrgb(1));
    set(gcf,'Color',bkgndrgb);
    set(gca,'XTick',[],'YTick',[],'visible','off');
    axis equal;
    drawnow;
    if (MAKEMOVIE)
        frame = getframe(gca);
        writeVideo(writerObj,frame);
    end
    if STILLFRAMES
        pause
    end
end
if (MAKEMOVIE)
    close(writerObj);
end

%%
% Section 4) Scatter plot of retinal positions tested
textfile = fullfile(nexfilepath,'nexfilelists','Greg','LMTF','ApolloTmpListLMTF.txt');
fnames = fnamesFromTxt2(textfile);
data = zeros(length(fnames),2);
for i = 1:length(fnames)
   stro = nex2stro(findfile(char(fnames{i})));
   data(i,:) = [abs(stro.sum.exptParams.stim_x) stro.sum.exptParams.stim_y];
end
uniqueXYs = unique(data,'rows','stable');
figure; axes; hold on; plot(uniqueXYs(:,1),uniqueXYs(:,2),'ko','MarkerFaceColor','black');
plot(0,0,'k*');
axis equal;

%%
% Section 5) SVD analysis of xi_lum, xi_rg, and theta for Apollo
% Need to load in the file saved in GrantBrainStorming Section 4.15

WHICHSV = 1; % Which singular vector to look at
MAKEMOVIE = 1;
INCLUDEPANUP = 0;
PLOTMODELFIT = 1;
COORDINATES = 'pol'; % 'pol' or 'cart'
moviefilename = 'RotSV';
elevations = linspace(90,25,100);
azimuths = linspace(0,360,200);

load BigApolloTmpListLMTFLMTFBrowserData
% The regressions should really be weighted regression
models = B.models;
labels = {'xi_{LUM}','zeta_{LUM}','n1_{LUM}','delta_n_{LUM}','tau1_{LUM}','kappa_{LUM}','xi_{RG}','zeta_{RG}','n1_{RG}','delta_n_{RG}','tau1_{RG}','kappa_{RG}','theta'};
Lnonzerovar = diag(cov(models')) > 10^-20;
covmat = cov(models(Lnonzerovar,:)');
[zscores,parammeans,paramsds] = zscore(models');
zscores = zscores';
[u,s,v] = svd(zscores(Lnonzerovar,:));
s = diag(s);

% Plotting one of the singular vectors
figure; axes; hold on;
bar(u(:,WHICHSV)); axis square; title(num2str(s(WHICHSV)));
set(gca,'Xtick',1:sum(Lnonzerovar),'XtickLabel',labels(Lnonzerovar)); 
print -dpsc junk;

% Making a spinning movie
if (INCLUDEPANUP)
    views = [zeros(1,length(elevations)); elevations];
else
    views = [];
end
views = [views, [azimuths; repmat(elevations(end),1,length(azimuths))]];

figure; axes; hold on;
set(gca,'Color',[.5 .5 .5]);hold on;
if strcmp(COORDINATES,'pol')
    [phi,r] = cart2pol(abs(B.eccs(:,1))/10,B.eccs(:,2)/10);
    plot3(phi,r,s(1)*v(:,WHICHSV),'ko','MarkerFaceColor','black');
    set(gca,'Ylim',[0 12]);
    if (PLOTMODELFIT)
        [b,bint,~,~,stats] = regress(s(1)*v(:,1),[ones(size(r)) r phi.^2]); % order b0, br, bphi
        [phimesh,rmesh] = meshgrid(-pi/2:pi/50:pi/2,0:.1:max(r));
        predprojs = [ones(size(rmesh(:))), rmesh(:), phimesh(:).^2]*b;
        h = surface(phimesh,rmesh,reshape(predprojs,size(rmesh)),'EdgeColor','none');
        set(h,'Facealpha',.7);

    end
else
    plot3(abs(B.eccs(:,1))/10,B.eccs(:,2)/10,s(1)*v(:,WHICHSV),'ko','MarkerFaceColor','black');
    axis equal
end
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
axis vis3d

if (MAKEMOVIE)
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(views)
        set(gca,'View',views(:,i));
        drawnow;
        frame = getframe(gca);
        writeVideo(writerObj,frame);
    end
    close(writerObj);
end

%%
% Section 6) Scatterplot of model errors for comparisons between simplifed
% and complicated models

load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
% First round vs. second round
figure;
%subplot(2,1,1); hold on;
%plot(A.legacy.firstroundfvs,A.legacy.secondroundfvs,'ko','MarkerFaceColor','black');
%plot([0 20],[0 20],'k--'); axis equal square; 
%xlabel('Error (364 parameters)');
%ylabel('Error (94 parameters)');
%propdeltaerr(1) = (A.legacy.secondroundfvs - A.legacy.firstroundfvs)./A.legacy.firstroundfvs;
properr(:,1) = A.legacy.secondroundfvs./A.legacy.firstroundfvs;
%subplot(2,1,2);
[n,x] = hist(log10(properr(:,1)),10);
bar(x,n);
ticks = 1:.5:ceil(max(properr));
set(gca,'Xtick',log10(ticks),'XTickLabel',ticks)
xlabel('error nested/error full')
axis square;
title(['Geo mean: ',num2str(geomean(properr(:,1)))]);

% Second round vs. third round
figure; 
properr(:,2) = A.legacy.thirdroundfvs./A.legacy.secondroundfvs;
%subplot(2,1,2);
[n,x] = hist(log10(properr(:,2)),10);
bar(x,n);
ticks = 1:.5:ceil(max(properr));
set(gca,'Xtick',log10(ticks),'XTickLabel',ticks)
xlabel('error nested/error full')
axis square;
title(['Geo mean: ',num2str(geomean(properr(:,2)))]);
axis square;

geomean(A.legacy.thirdroundfvs./A.legacy.firstroundfvs)
figure
plot([364 94 17],[mean(A.legacy.firstroundfvs) mean(A.legacy.secondroundfvs) mean(A.legacy.thirdroundfvs)],'ko','MarkerFaceColor','black');
set(gca,'Ylim',[0 7]);
ylabel('Sum(abs(residuals))');
xlabel('# parameters');
%%
% Section 7: Superimposed fits on data from a single retinal location
CCPLTLIM = .3;
MAKEMOVIE = 1;
PLOTMODELBOOL = [0 0 0]; % Which models to plot: 1, 2, and/or 3

moviefilename = 'onefits';
%load C:/Users/ghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat;
stimloc = [0 0]; % initializing
i = 1;
while ~(stimloc(1) == 60 & stimloc(2) == 0)
    stimloc = A.raw{i}(1,[5 6]);
    i = i + 1;
end
data = A.raw{i}(:,1:4);
model1 = A.legacy.firstroundmodels(:,i);
model2 = A.legacy.secondroundmodels(:,i);
model3 = A.legacy.thirdroundmodels(:,i);

figure; axes; hold on;
Loog = data(:,4) == 1;
plot3(data(~Loog,1),data(~Loog,2),data(~Loog,3),'ko','MarkerFaceColor','black');
plot3(-data(~Loog,1),-data(~Loog,2),data(~Loog,3),'ko','MarkerFaceColor','black');
plot3(data(Loog,1),data(Loog,2),data(Loog,3),'ro','MarkerFaceColor','red');
plot3(-data(Loog,1),-data(Loog,2),data(Loog,3),'ro','MarkerFaceColor','red');
set(gca,'Zscale','log'); set(gca,'Xlim',[-CCPLTLIM CCPLTLIM],'Ylim',[-CCPLTLIM CCPLTLIM],'Zlim',[1 25]);
set(gca,'View',[-45 20]);
axis square
set(gca,'Zticklabel',[1 10]);

facecolors = {'green','blue','red'};
for i = find(PLOTMODELBOOL)
    if i == 1
        fpar = model1;
    elseif i == 2
        fpar = model2;
    elseif i == 3
        fpar = model3;
    end
    [xx,yy,zz] = meshgrid(linspace(-max(abs(data(:,1))),max(abs(data(:,1))),20),...
    linspace(-max(abs(data(:,2))),max(abs(data(:,2))),20),...
    linspace(min(data(:,3)),max(data(:,3)),20));

    f1 = @(omega)fpar(1)*abs(((1i*2*pi*10^fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*10.^(fpar(5)+fpar(6)).*omega+1).^-(fpar(3)+fpar(4))));
    f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*10^fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*10.^(fpar(5+6)+fpar(6+6)).*omega+1).^-(fpar(3+6)+fpar(4+6))));

    a = abs(f1(zz)).^-1; % chromatic threshold
    b = abs(f2(zz)).^-1; % luminance threshold
    thtmp = atan2(yy,xx)-fpar(13); % clockwise rotation from [L,M] to [a,b]
    rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    
    V = sqrt(xx.^2+yy.^2)-rtmp;
    FV = isosurface(xx,yy,zz,V,0);

    h = patch(FV,'facecolor',char(facecolors{i}));
    set(h,'Edgecolor','none'); % camlight headlight; lighting phong
    alpha(.2);
end


if (MAKEMOVIE)
    plot3([0 0],[0 0],[1 25],'k-');
	plot3([CCPLTLIM CCPLTLIM -CCPLTLIM  -CCPLTLIM CCPLTLIM],[-CCPLTLIM CCPLTLIM CCPLTLIM -CCPLTLIM -CCPLTLIM],[1 1 1 1 1],'k-');
    set(gcf,'Color','white')
    axis vis3d
    originalview = get(gca,'View');
    deltaviews = linspace(0,360,120);
    deltaviews(end) = [];
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[deltaviews(i) 0]);
        set(gca,'xtick',[-.2 0 .2], 'ytick',[-.2 0 .2]);
        drawnow;
        frame = getframe(gca);
        writeVideo(writerObj,frame);
    end
    close(writerObj);
end

%%
% Section 8
% Model fit of luminance and chromatic sensitivity across the visual field
load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
params = U.legacy.mode0models;
eccs = U.eccs;
TF = 6.18;
xi_lum_coeffs = params(1,:);
xi_rg_coeffs = params(7,:);
[~,~,tmp] = tf_fiterr2([1; params(2:6,1); 1; params(8:13,1)],[1 1 TF 0]); % Height of fitted LUM TCSF at "TF" Hz assuming xi_lum = 1
XiLUMScaleFactor = 1/tmp; % sensitivity peak, setting xi_lum to 1. Third output argument from ft_fiterr2 is predicted *threshold*.
[~,~,tmp] = tf_fiterr2([1; params(2:6,1); 1; params(8:13,1)],[1 -1 TF 0]); % Height of fitted LUM TCSF at "TF" Hz assuming xi_lum = 1
XiRGScaleFactor = 1/tmp; % sensitivity peak, setting xi_lum to 1. Third output argument from ft_fiterr2 is predicted *threshold*.

[phi_ecc,r_ecc] = cart2pol(abs(eccs(:,1))/10,eccs(:,2)/10);
b_lum = regress(log10(xi_lum_coeffs'*XiLUMScaleFactor),[ones(numel(phi_ecc),1),r_ecc(:),r_ecc(:).*sin(2*phi_ecc(:)),r_ecc(:).*cos(2*phi_ecc(:))]);
b_rg = regress(log10(xi_rg_coeffs'*XiRGScaleFactor),[ones(numel(phi_ecc),1),r_ecc(:),r_ecc(:).*sin(2*phi_ecc(:)),r_ecc(:).*cos(2*phi_ecc(:))]);

[x,y] = meshgrid(linspace(0,10,25),linspace(-10,10,50));
[phi,r] = cart2pol(x,y);
X = [ones(numel(r),1), r(:), r(:).*sin(2.*phi(:)), r(:).*cos(2.*phi(:))];
log_lum_sens = reshape(X*b_lum, size(r));
log_rg_sens = reshape(X*b_rg, size(r));

% Both LUM and RG scaled identically
figure;
h = [];
for i = 1:2
    subplot(1,2,i);
    if (i == 1)
        h(i) = surf(x,y,10.^log_lum_sens);
    else
        h(i) = surf(x,y,10.^log_rg_sens);
    end
    axis tight;
    axis equal;
    set(h,'EdgeAlpha',0);
    set(gca,'View',[0 90]); 
    caxis([min(10.^[log_lum_sens(:); log_rg_sens(:)]) max(10.^[log_lum_sens(:); log_rg_sens(:)])]);
end
set(gcf,'Renderer','painters')
colormap(hot);

% Insets show slope of contrast sensitivity decrease as a function of
% phi
phi = linspace(-pi/2,pi/2,100);
X1 = [ones(numel(phi),1), ones(numel(phi),1), sin(2.*phi(:)), cos(2.*phi(:))];
X2 = [ones(numel(phi),1), 2*ones(numel(phi),1), 2*sin(2.*phi(:)), 2*cos(2.*phi(:))];

log_lum_sens1 = reshape(X1*b_lum, size(phi)); % This is the xi_lum parameter
log_rg_sens1 = reshape(X1*b_rg, size(phi)); % This is the xi_rg parameter

log_lum_sens2 = reshape(X2*b_lum, size(phi)); % This is the xi_lum parameter
log_rg_sens2 = reshape(X2*b_rg, size(phi)); % This is the xi_rg parameter


% % Testing. Calculating the slope in two different ways. Checks out.
% X3 = [ones(numel(phi),1), 3*ones(numel(phi),1), 3*sin(2.*phi(:)), 3*cos(2.*phi(:))];
% xi_lum3 = reshape(X3*b_lum, size(phi)); % This is the xi_lum parameter
% xi_rg3 = reshape(X3*b_rg, size(phi)); % This is the xi_rg parameter
% log_lum_sens3 = xi_lum3*abs(((1i*2*pi*10^params(4).*tf+1).^-params(2))-params(1)*((1i*2*pi*10^(params(4)+params(5)).*tf+1).^-(params(2)+params(3))));
% log_rg_sens3 = xi_rg3*abs(((1i*2*pi*10^params(9).*tf+1).^-params(7))-params(6)*((1i*2*pi*10^(params(9)+params(10)).*tf+1).^-(params(7)+params(8))));
% (log_lum_sens3'-log_lum_sens2') - (log_lum_sens2'-log_lum_sens1') % zero up to round off error.

figure; axes; 
SensDecayFactor = 10.^(log_lum_sens2-log_lum_sens1);
plot(phi,SensDecayFactor,'k-','Linewidth',2);
axis square; set(gca,'xlim',[-pi/2 pi/2],'Xtick',[-pi/2 -pi/4 0 pi/4 pi/2]);
ylabel('contrast sensitivity scaling/�');
set(gca,'YTick',[.9:.025:.975],'Ylim',[.9 .975],'Tickdir','out');

figure; axes; 
SensDecayFactor = 10.^(log_rg_sens2-log_rg_sens1);
plot(phi,SensDecayFactor,'r-','Linewidth',2);
axis square; set(gca,'xlim',[-pi/2 pi/2],'Xtick',[-pi/2 -pi/4 0 pi/4 pi/2]);
set(gca,'YTick',[.825:.025:.95],'Ylim',[.825 .95],'Tickdir','out');
ylabel('contrast sensitivity scaling/�');

% Sanity checks
% Computing what the sensitivity should be 10 degrees out from fovea
33*exp((-(1-.92)*10))
33*exp((-(1-.845)*10))
23*exp((-(1-.95)*10))
23*exp((-(1-.915)*10))

%%
% Section 9
% An example luminance temporal contrast senstivity function with data
RFX = 5;
RFY = 0;
STIMTYPE = 'LUM';
SMALLAXES = 1;
if strcmp(STIMTYPE,'RG')
    angle = 3*pi/4;
else
    angle = pi/4;
end
angletol= pi/7; % arbitrary
isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
% "A" is the structure of LMTF model parameters for Apollo
behavioralmodel = LMTF_global_to_local_model(A.legacy.mode5params, RFX, RFY, 5);
if strcmp(STIMTYPE,'LUM')
    behavioral_tcsf = @(omega)(behavioralmodel(1)*abs(((1i*2*pi*10^behavioralmodel(5).*omega+1).^-behavioralmodel(3))-behavioralmodel(2)*((1i*2*pi*10^(behavioralmodel(5)+behavioralmodel(6)).*omega+1).^-(behavioralmodel(3)+behavioralmodel(4)))));
elseif strcmp(STIMTYPE,'RG')
    behavioral_tcsf = @(omega)(behavioralmodel(7)*abs(((1i*2*pi*10^behavioralmodel(11).*omega+1).^-behavioralmodel(9))-behavioralmodel(8)*((1i*2*pi*10^(behavioralmodel(11)+behavioralmodel(12)).*omega+1).^-(behavioralmodel(9)+behavioralmodel(10)))));
end
w = logspace(log10(1),log10(60),200);
L = all(A.eccs == repmat(10*[RFX RFY],size(A.eccs,1),1),2);
raw = A.raw{L};
% Getting raw data near the luminance axis
Loog = logical(raw(:,4));
[th,r] = cart2pol(raw(:,1),raw(:,2));
L = th<angle+angletol/2 & th>angle-angletol/2;
L = L & ~Loog;

figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
plot(w,behavioral_tcsf(w),'y-','LineWidth',2);
set(gca,'Xscale','log','Yscale','log','Ylim',[.1 1000]);
if SMALLAXES
    set(gca,'Xlim',[1 30],'Ylim',[1 1000]);
end
h = plot(raw(L,3),1./r(L),'yo');
set(h,'MarkerFaceColor','yellow')

%%
% Section 9.1 
% A few other sections.
% theta, X, Y, 

isosamppath = which('IsoSampOnline');
isosamppath(find(isosamppath==filesep,1,'last'):end) = [];
load ([isosamppath,filesep,'private',filesep','data',filesep,'LMTF.mat']);
globalmodel = A.legacy.mode5params;

% Getting all the raw data points
rawdata = [];
for i=1:length(A.raw)
    rf = A.eccs(i,:);
    rawtmp = A.raw{i};
    rawdata = [rawdata; repmat(rf,size(rawtmp,1),1),rawtmp(:,1:4)]; % RFX RFY Lcc Mcc TF Loog
end
[thetas, rs] = cart2pol(rawdata(:,4),rawdata(:,3));


% % First, theta
% TF = 5;
% TFmarginfactor = .7;
% RFX = 10;
% RFY = 0;
% behavioralmodel = LMTF_global_to_local_model(A.legacy.mode5params, RFX, RFY, 5);
% tmp = [cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))' repmat(8,100,1)];
% [~, ~, preds] = tf_fiterr2(behavioralmodel, tmp)
% 
% figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
% set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
% plot(preds.*tmp(:,1), preds.*tmp(:,2),'w-','LineWidth',2);
% L = abs(rawdata(:,1)) == RFX*10 & rawdata(:,2) == RFY*10 & rawdata(:,5) > TF*TFmarginfactor & rawdata(:,5) < TF/TFmarginfactor & rawdata(:,6) ~= 1;
% plot(rawdata(L,3),rs(L),'wo','MarkerFaceColor','white');
% plot(-rawdata(L,3),-rs(L),'wo','MarkerFaceColor','white');
% axis equal
% 
% % Horizontal and vertical meridians
% TF = 5;
% theta = pi/4;
% thetamargin = pi/7;
% RFpositions = A.eccs(A.eccs(:,2) == 0,:)/10;
% RFpositions = [linspace(2,10,10)' zeros(10,1)];
% preds = zeros(size(RFpositions,1),1);
% for i = 1:size(RFpositions, 1)
%     behavioralmodel = LMTF_global_to_local_model(A.legacy.mode5params, RFpositions(i,1), RFpositions(i,2), 5);
%     [~, ~, pred] = tf_fiterr2(behavioralmodel, [cos(theta) sin(theta) TF 0])
%     preds(i) = pred;
% end
% figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
% set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
% plot(RFpositions(:,1),1./preds,'w-','LineWidth',2);
% L = rawdata(:,2) == 0 & rawdata(:,5) > TF*TFmarginfactor & rawdata(:,5) < TF/TFmarginfactor & rawdata(:,6) ~= 1 & thetas > theta-thetamargin &  thetas < theta+thetamargin ;
% plot(rawdata(L,1)/10,1./rs(L),'wo','MarkerFaceColor','white');

% TCSFs
w_lum = logspace(log10(1),log10(60),200);
w_rg = w_lum(w_lum<20);
angletol= pi/7; % arbitrary
for i = 1:size(A.eccs,1)
    RFX = A.eccs(i,1);
    RFY = A.eccs(i,2);
    figure; axes; hold on; set(gcf,'Color',[0 0 0],'InvertHardCopy','off');
    set(gca,'Color',[0 0 0],'TickDir','out','YColor',[1 1 1],'XColor',[1 1 1]);
    behavioralmodel = LMTF_global_to_local_model(A.legacy.mode5params, RFX/10, RFY/10, 5);
    for j = 1:2 % LUM, RG
        if j == 1
            angle = pi/4; % LUM
        else
            angle = 3*pi/4; % RG
        end
        L = all(A.eccs == repmat([RFX RFY],size(A.eccs,1),1),2);
        raw = A.raw{L};
        % Getting raw data near the luminance axis
        Loog = logical(raw(:,4));
        [th,r] = cart2pol(raw(:,1),raw(:,2));
        L = th<angle+angletol/2 & th>angle-angletol/2;
        L = L & ~Loog;
        if j == 1
            behavioral_tcsf = @(omega)(behavioralmodel(1)*abs(((1i*2*pi*10^behavioralmodel(5).*omega+1).^-behavioralmodel(3))-behavioralmodel(2)*((1i*2*pi*10^(behavioralmodel(5)+behavioralmodel(6)).*omega+1).^-(behavioralmodel(3)+behavioralmodel(4)))));
            plot(w_lum,behavioral_tcsf(w_lum),'y-','LineWidth',2);
            h = plot(raw(L,3),1./r(L),'yo');
            set(h,'MarkerFaceColor','yellow');
        else
            behavioral_tcsf = @(omega)(behavioralmodel(7)*abs(((1i*2*pi*10^behavioralmodel(11).*omega+1).^-behavioralmodel(9))-behavioralmodel(8)*((1i*2*pi*10^(behavioralmodel(11)+behavioralmodel(12)).*omega+1).^-(behavioralmodel(9)+behavioralmodel(10)))));
            plot(w_rg,behavioral_tcsf(w_rg),'r-','LineWidth',2);
            h = plot(raw(L,3),1./r(L),'ro');
            set(h,'MarkerFaceColor',[1 .2 .2],'MarkerEdgeColor',[1 .2 .2]);
        end
    end
    set(gca,'Xscale','log','Yscale','log','Ylim',[.1 1000]);
    title([num2str(RFX),', ',num2str(RFY)],'Color','white');
    set(gcf,'Renderer','painters')
end
%%
% 

%%
% Section 10
% 1-D contrast sensitivity functions for humans and monkeys
TWOPLOTS = 1;
mode = 5;
RF = [5 0];
angletol = pi/7;
load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat
SIDs = {'A','U','E','G'};
colors = [1 .5 0; .75 0 .75; 0 0 .75; 0 0 0];
lumtfs = logspace(log10(1),log10(60),100);
rgtfs = logspace(log10(1),log10(20),100);
figure;
h = [];
for j = 1:length(SIDs)
    s = eval(SIDs{j});
    globalparams = getfield(s.legacy,['mode',num2str(mode),'params']);
    model = LMTF_global_to_local_model(globalparams, RF(1), RF(2), mode);
    disp(model(13))
    lumf1 = @(omega)(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
    rgf1 = @(omega)(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
    
    subplot(1,2,1); hold on;
    h(j) = plot(lumtfs,lumf1(lumtfs),'-','LineWidth',2,'Color',colors(j,:));

    % Plotting raw data
    L = all(s.eccs == repmat(10*[RF(1) RF(2)],size(s.eccs,1),1),2);
    raw = s.raw{L};
    Loog = logical(raw(:,4));
    [th,r] = cart2pol(raw(:,1),raw(:,2));
    % Getting raw data near the luminance axis
    L = th<model(13)+angletol/2 & th>model(13)-angletol/2;
    L = L & ~Loog;
    plot(raw(L,3), 1./r(L),'o','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:));

    subplot(1,2,2); hold on;
    plot(rgtfs,rgf1(rgtfs),'-','LineWidth',2,'Color',colors(j,:));

    % Getting raw data near the chromatic axis
    L = th<3*pi/4+angletol/2 & th>3*pi/4-angletol/2;
    L = L & ~Loog;
    plot(raw(L,3), 1./r(L),'o','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:));
end
for j = 1:2
    subplot(1,2,j)
    title(['RF @ (',num2str(RF),')']);
    set(gca,'xscale','log','Yscale','log');
    set(gca,'Xlim',[1 20],'Ylim',[5 40]);
    xlabel('Temporal frequency (Hz)');
    ylabel('Contrast sensitivity');
    if i == 2
        legend(h,SIDs);
    end
end

%%
% Section 11
% Example rasters
filename = 'A061517005.nex'; % example magnocell
spikeNum = 'sig001a';
stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Greg',filesep,'Apollo']));

Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'rew_t'));

dur = mean(stimoff_t-stimon_t);
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikeNum);
spikes = stro.ras(:,spikeidx);
uniquestim = IsoSampGetDPrime(stro,1);
Lstimtype = (sign(uniquestim(:,1)) == sign(uniquestim(:,2))); % includes blank

% Rasters for all (Lstimtype) conditions
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
figure; axes; hold on; 
set(gca,'TickDir','out');
counter = 0;
for j = find(Lstimtype)'
    L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
    for i = find(L)'
        tmpspikes = spikes{i}-stimon_t(i);
        tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
        nspikestot = length(tmpspikes);
        plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
        counter = counter + 1;
    end
    plot([offset(1) dur+offset(2)],counter*[1 1],'b:');
    h = text(-.12,counter-sum(L)/2,num2str(uniquestim(j,3),2),'FontSize',12,'HorizontalAlignment','right')
end
set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[],'Box','on');
set(gcf,'Renderer','painters')

%%
% Section 12
% Pie charts showing the number of subjects tested at each retinal location
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
location_query = 'SELECT rfX, rfY, subjID FROM LMTF WHERE recDate > ''2016-05-27'' AND recDate < ''2017-02-03'' AND quality = 1 ORDER BY `LMTF`.`subjID` ASC';
alldata = fetch(conn, location_query);
close(conn);
allsubjects = cell2mat(alldata(:,3));
rfxy = cell2mat(alldata(:,[1 2]));
uniquexys = unique(rfxy,'rows');

% Setting up dummy variables for subjects
subjects = {'A','U','E','G'};
colors = [.9 .4 .2;0 0 1; 1 0 1; 0 0 0];
n_subjects = length(subjects);
n_files = zeros(size(uniquexys,1),length(subjects));

for i = 1:size(uniquexys,1)
    L = rfxy(:,1) == uniquexys(i,1) & rfxy(:,2) == uniquexys(i,2);
    for j = 1:length(subjects)
        n_files(i,j) = sum(L & allsubjects == subjects{j});
    end
end

figure; axes; hold on;
r = 5;
tmp = [cos(linspace(0,2*pi,100))' sin(linspace(0,2*pi,100))'];
for i = 1:size(uniquexys,1)
    Lsubject = n_files(i,:)> 0;
    nsubjects = sum(Lsubject);
    %plot(uniquexys(i,1)+tmp(:,1)*r, uniquexys(i,2)+tmp(:,2)*r,'k-');
    for j = 1:nsubjects
        subjectidxs = find(Lsubject);
        tmp = linspace((j-1)/nsubjects,j/nsubjects,50);
        tmpx = [0 r*cos(2*pi*tmp)];
        tmpy = [0 r*sin(2*pi*tmp)];
        h = patch(uniquexys(i,1)+tmpx,uniquexys(i,2)+tmpy,colors(subjectidxs(j),:));
        set(h,'EdgeColor',colors(subjectidxs(j),:))
    end
end
axis equal

sum(n_files > 0) 
%%
% Section 13
% Analysis of residuals from first round fits.

Lidx = 1;
Midx = 2;
TFidx = 3;
Xidx = 4;
Yidx = 5;
threshidx = 6;
predidx = 7;

load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
subjects = {'A','U','E','G'};
colors = [.9 .4 .2;0 0 1; 1 0 1; 0 0 0];
figure; axes; hold on;
for i = 1:length(subjects)
    s = eval(subjects{i});
    neccs = length(s.raw);
    data = [];
    for j = 1:neccs
        Loog = logical(s.raw{j}(:,4));
        [~,r,f] = tf_fiterr2(s.legacy.firstroundmodels(:,j), s.raw{j});
        data = [data; s.raw{j}(~Loog,[1:3,5:6]) r(~Loog) f(~Loog)];
    end  
    resid = log10(data(:,threshidx))-log10(data(:,predidx)); 
    h = plot(log10(data(:,predidx)), resid,'o','MarkerSize',2);
    set(h,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
    
    % Looking at mean and std of resids
    binwidth = .1;
    bins = [-1.6:binwidth:-.4];
    v = [];
    sd = [];
    for j = 1:length(bins)
        m(j) = mean(resid(log10(data(:,predidx)) > (bins(j)-binwidth/2) & log10(data(:,predidx)) < (bins(j)+binwidth/2)));
        sd(j) = std(resid(log10(data(:,predidx)) > (bins(j)-binwidth/2) & log10(data(:,predidx)) < (bins(j)+binwidth/2)));
    end
    %plot(bins,sd,'color',colors(i,:),'LineWidth',2)
end
ylabel('Residual','FontSize',12);
set(gca,'Xtick',log10([.01, .03, .1, .3, 1]));
set(gca,'XtickLabel',[.01, .03, .1, .3, 1]);
xlabel('Predicted threshold','FontSize',12);
set(gca,'Ylim',[-2.5 1],'Xlim',[log10(0.01) .2],'TickDir','out');
%%
% Section 14
% Analysis of residuals from final model fits

Lidx = 1;
Midx = 2;
TFidx = 3;
Xidx = 4;
Yidx = 5;
threshidx = 6;
predidx = 7;

load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
subjects = {'A','U','E','G'};
colors = [.9 .4 .2;0 0 1; 1 0 1; 0 0 0];
hfig(1) = figure; axes; hold on;
hfig(2) = figure; axes; hold on;
maxabsresid = 0;
RESIDPLOTSCALE = 60;
allsubjs_thresh_pred = [];
for i = 1:length(subjects)
    s = eval(subjects{i});
    neccs = length(s.raw);
    data = [];
    for j = 1:neccs
        Loog = logical(s.raw{j}(:,4));
        [~,r,f] = tf_fiterr2(s.legacy.mode5models(:,j), s.raw{j});
        data = [data; s.raw{j}(~Loog,[1:3,5:6]) r(~Loog) f(~Loog)];
    end  
    resid = log10(data(:,threshidx))-log10(data(:,predidx));
    allsubjs_thresh_pred = [allsubjs_thresh_pred; repmat(i,size(data,1),1) data(:,threshidx) data(:,predidx)];
    % Residual as a function of predicted threshold
    figure(hfig(1));
    h = plot(log10(data(:,predidx)), resid,'o','MarkerSize',2);
    set(h,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    
    % Cross correlation of residuals as a function of TF and color
    % direction
    figure(hfig(2));
    subplot(2,length(subjects),(i-1)*2+1);
    colordir = cart2pol(data(:,Lidx),data(:,Midx));
    thetabinedges = linspace(0,pi,9);
    thetabinedges = thetabinedges-(thetabinedges(2)-thetabinedges(1))/2;
    thetabinmidpoints = (thetabinedges(2:end)+thetabinedges(1:end-1))/2;
    colordir(colordir>thetabinedges(end)) = colordir(colordir>thetabinedges(end))-pi;
    logTFbinedges = linspace(log10(1),log10(15),11);
    theta_TF_matrix = zeros(length(thetabinedges)-1,length(logTFbinedges)-1);
    % First calculating mean residuals in each bin
    for j = 1:length(thetabinedges)-1
        for k = 1:length(logTFbinedges)-1
            L = colordir>=thetabinedges(j) & colordir<thetabinedges(j+1) & log10(data(:,TFidx))>=logTFbinedges(k) & log10(data(:,TFidx))<logTFbinedges(k+1);
            theta_TF_matrix(j,k) = nanmedian(resid(L));
        end
    end % nans where there are no residuals
    % Now calculating cross-correlation
    tmpmatrix = theta_TF_matrix;
    rs = zeros(length(thetabinedges)-1,length(logTFbinedges)-1);
    for j = 1:length(thetabinedges)-1
        for k = 1:length(logTFbinedges)-1
            tmpmatrix = theta_TF_matrix;
            if j > 1 % circular scrolling
                tmpmatrix = tmpmatrix([j:end,1:(j-1)],:);
            end
            if k > 1
                tmpmatrix = [nan(length(thetabinedges)-1,k-1),tmpmatrix(:,1:end-(k-1))];
            end
            centereddata = [(theta_TF_matrix(:)-nanmean(theta_TF_matrix(:))), (tmpmatrix(:)-nanmean(tmpmatrix(:)))];
            centereddata = centereddata(~any(isnan(centereddata),2),:);
            r = corrcoef(centereddata);
            rs(j,k) = r(1,2);
            %if j == 1 & k == 1
            %    rs(j,k) = nan; % autocorrelation at zero lag is trivially large
            %end
        end
    end
    % Need to bring correlations into the range 1-64 before calling image
    colormap([linspace(0,1,64)',linspace(0,1,64)', linspace(1,0,64)']);
    colormap(parula(64));
    h = image(rs'*32+33);
    
    axis tight
    axis xy
    thetabinwidth = thetabinedges(2)-thetabinedges(1);
    thetabincenters = thetabinedges(1:end-1)+thetabinwidth/2;
    bincenterstmp = interp1(thetabincenters,1:length(thetabincenters),[0, pi/2 ],'linear','extrap');
    set(gca,'Xtick',bincenterstmp,'Xticklabel',{'0','\pi/2'})
    set(gca,'Ytick',[0:2:length(logTFbinedges)],'Yticklabel',round(10.^logTFbinedges(1:2:end)))
    ylabel('TF'); xlabel('theta');
    
    % residual as a function of screen location
    subplot(2,length(subjects),(i-1)*2+2); hold on;
    for j = 1:neccs
        L = all(data(:,[Xidx Yidx]) == repmat(s.eccs(j,:),size(data,1),1),2);
        resids = log10(data(L,threshidx))-log10(data(L,predidx));
        netresid = nanmedian(resids);
        if abs(netresid) > maxabsresid
            maxabsresid = netresid;
        end
        h = plot(s.eccs(j,1)/10,s.eccs(j,2)/10,'o');
        if (netresid > 0)
            set(h,'MarkerSize',RESIDPLOTSCALE*netresid+1,'Color','black','MarkerFaceColor','black');
        else
            set(h,'MarkerSize',-RESIDPLOTSCALE*netresid+1,'Color','red','MarkerFaceColor','red');
        end
    end
    set(gca,'Plotboxaspectratio',[1 2 1])
    set(gca,'Xlim',[0 12],'Ylim',[-12 12]);
    
    xlabel('X');
    ylabel('Y');
end

% Colorbar
figure; 
axes('position',[0.46 0.05 .06 .9]);
image([1:64]');
axis xy
colormap(parula(64));
set(gca,'Ytick',linspace(1,64,5),'Yticklabel',[-1 -.5 0 .5 1],'Xtick',[],'TickDir','out');

% Symbol bar
subplot(2,4,8);% Same size as residual dot plot
set(gca,'Plotboxaspectratio',[1 2 1])
set(gca,'Xlim',[-6 6],'Ylim',[-12 12]);
hold on;
symbolsizes = RESIDPLOTSCALE*maxabsresid*linspace(-1,1,6)+sign(linspace(-1,1,6));
ypositions = 10*linspace(-1,1,6);
for i = 1:length(symbolsizes)
   h = plot(0,ypositions(i),'o','MarkerSize',abs(symbolsizes(i)));
   if symbolsizes(i) > 0
       set(h,'MarkerFaceColor','black','MarkerEdgeColor','black');
   else
       set(h,'MarkerFaceColor','red','MarkerEdgeColor','red');
   end
   text(2,ypositions(i),num2str(10.^(symbolsizes(i)/RESIDPLOTSCALE),3))
end
set(gca,'Xtick',[],'Ytick',[]);

figure(hfig(1));
ylabel('Residual','FontSize',12);
set(gca,'Xtick',log10([.01, .03, .1, .3, 1]));
set(gca,'XtickLabel',[.01, .03, .1, .3, 1]);
xlabel('Predicted threshold','FontSize',12);
set(gca,'Ylim',[-2.5 1],'Xlim',[log10(0.01) .2],'TickDir','out');

% Reviewer 2 asked for an analysis of overall residuals
% columns of allsubjs_thresh_pred: (1) subject ID, (2) threshold, (3)
% predicted threshsold.
% Not including out of gamut points
hist(allsubjs_thresh_pred(:,2)./allsubjs_thresh_pred(:,3))
% Looking t upper and lower percentiles
prctile(allsubjs_thresh_pred(:,2)./allsubjs_thresh_pred(:,3),[5 10 25 50 75 90 95])


%%
% Section 15: Low temporal frequency luminance sensitivity for monkeys and
% humans. Emily sent these data on 2/28/18 (using "senscomp_mvh.m"). 
% Thresholds are for L+M modulations from 1-1.459 Hz.

A = [ 0.5043
    0.6370
    0.4326
    0.4956
    0.4956
    0.6084
    0.6336
    0.4869
    0.4171
    0.5971
    0.3494];
U = [0.6272
    0.5018
    0.4799
    0.7463
    0.4033
    0.4363
    0.6651
    0.4561
    0.4825
    0.5118
    0.5563
    0.5966
    0.3793
    0.6051
    0.4043
    0.5477
    0.2965
    0.4558
    0.5242
    0.4230
    0.4607
    0.4357
    0.3932
    0.3739
    0.0726
    0.3567
    0.1704
    0.4218];

E = [0.6185
    0.5281
    0.6609
    0.7196
    0.6483
    0.6918
    0.5701
    0.8405
    1.0179
    0.5494
    0.4848
    0.5511
    0.5667
    0.6918
    0.8514];
G = [0.6152
    0.8888
    0.5081
    0.5036
    0.5117
    0.5192
    0.5634
    0.5477
    0.5364
    0.6233
    0.5928
    0.6332
    0.7811
    0.6924
    0.6606
    0.5340
    0.5902
    0.5902
    0.8022];

data = [A;U;E;G];

subj={};
[subj{1:length(A),1}] = deal('A');
[subj{length(subj)+1:length(subj)+length(U),1}] = deal('U');
[subj{length(subj)+1:length(subj)+length(E),1}] = deal('E');
[subj{length(subj)+1:length(subj)+length(G),1}] = deal('G');

species = {};
[species{1:length(A)+length(U),1}] = deal('Monkey');
[species{length(species)+1:length(species)+length(G)+length(E),1}] = deal('Human');

[p,tbl] = anovan(data,{species,subj},'model','linear','varnames',{'species','subject'},'nested',[0,0;1 0],'random',2); % subject nested within species, subject random


%%
% Section 16
% Data for Table 2 (containing final model parameters)
% Order: zeta_lum, n1_lum, n2_lum, tau1_lum, tau2_lum, zeta_rg, n1_rg,
% n2_rg, tau1_rg, tau2_rg, theta, b0_lum, b1_lum, b2_lum, b0_rg, b1_rg,
% b2_rg, b3.

load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
subs= {'A','U','E','G'};
for i = 1:length(subs)
    params = eval([char(subs(i)),'.legacy.mode5params'])
    [params(1:2); params(2)+params(3); 10^params(4); 10^(params(4)*params(5)); 
     params(6:7); params(7)+params(8); 10^params(9); 10^(params(9)*params(10)); 
     params(11);
     params(12:14); params(16:18); params(15)]
    sprintf('%2.2f\n',params)
end

%% 
% Section 17
% Figure to the reviewers showing konio cell responses to L+M
% and optionally magno responses to L+M
WNfilenames = {'A060317001.nex'}; % konio
IsoSampfilenames = {'A060317002.nex'}; % konio
WNfilenames = {'A061717010.nex'}; % magno
IsoSampfilenames = {'A061717011.nex'}; % magno

stro = {};
for i = 1:length(WNfilenames)
    stro = nex2stro(findfile(WNfilenames{i}));
end
WN = strocat(stro);
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
% Reconstructing the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = WN.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
STAs = out{1};
nspikes = out{3};

temporal_energy = sum(STAs.^2);
whichframe = find(temporal_energy == max(temporal_energy));
STA = reshape(STAs(:,whichframe),nstixperside,nstixperside,3);
figure; subplot(2,2,1);
image(STA./(2*max(abs(STA(:))))+.5);
axis square;
set(gca,'Xtick',[],'Ytick',[]);

% Now IsoSamp data (GDLH stopped here)
stro = {};
for i = 1:length(IsoSampfilenames)
    stro = nex2stro(findfile(IsoSampfilenames{i}));
end
IsoSamp = strocat(stro);
Lcc = IsoSamp.trial(:,strcmp(IsoSamp.sum.trialFields(1,:),'stim_l'));
Mcc = IsoSamp.trial(:,strcmp(IsoSamp.sum.trialFields(1,:),'stim_m'));
TF = IsoSamp.trial(:,strcmp(IsoSamp.sum.trialFields(1,:),'tf'));
stimon_t = IsoSamp.trial(:,strcmp(IsoSamp.sum.trialFields(1,:),'stimon_t'));
stimoff_t = IsoSamp.trial(:,strcmp(IsoSamp.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
spikeidxs = strncmp(IsoSamp.sum.rasterCells(1,:),'sig',3);
uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
bins = linspace(0,dur,100);
data = zeros(size(uniquestim,1), length(bins));
temporalenvelope = ones(size(bins));
temporalenvelope(1:round((1/4)*length(bins))) = linspace(0,1,round((1/4)*length(bins)));
temporalenvelope(end:-1:round((3/4)*length(bins))+1) = linspace(0,1,round((1/4)*length(bins)));

for i = 1:size(uniquestim,1)
    contrast = sqrt(uniquestim(i,1).^2+uniquestim(i,2).^2);
    data(i,:) = contrast*temporalenvelope.*sin(2*pi*uniquestim(i,3)*bins);
end

for spikeidx = find(spikeidxs)
    spikes = IsoSamp.ras(:,spikeidx);
    Lstimtype = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
    subplot(2,2,3);
    im = 255*(data - min(data(:)))/(max(data(:)) - min(data(:)));
    image(flipud(im(Lstimtype,:)));
    xtick = 0:.2:dur;
    set(gca,'XTick',interp1([bins(1) bins(end)],[1 length(bins)],xtick));
    set(gca,'Xticklabel',xtick);
    
    colormap(jet(255));
    axis square;
    set(gca,'Ytick',[]);
    
    % Plotting rasters
    offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
    subplot(2,2,4); hold on; counter = 0;
    for j = find(Lstimtype)'
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        for i = find(L)'
            tmpspikes = spikes{i}-stimon_t(i);
            tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
            nspikestot = length(tmpspikes);
            plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
            counter = counter + 1;
        end
    end
    set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
    axis square;
end
set(gcf,'name',stro.sum.rasterCells{spikeidx});
set(gcf,'Renderer','painters');
print -dpdf junk2

%%
% Section 18
% See section 2