% Figures for the PBIO retreat, 2016
% A lot of dependences in this code. Don't expect to be able to 
% run one cell from this middle of this thing and have it work.
%
% Section 1) A simulated cone mosaic
% Section 2) Spinning cone thresholds
% Section 3) Monkey behavior on the same axes as the cone thresholds(?)
% Section 4) Monkey neurophysiological data
% Section 5) Example Gabor (refers to iSETBIOFigs.m, section 5)
% Section 6) Cone thresholds and monkey behavioral thresholds on the same
% axes
% Section 7) Cone thresholds with fits (spinning movie)
% Section 8) Monkey behavioral thresholds with fits (spinning movie)


%%
% Section 1) A simulated cone mosaic
addpath('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg/SamplingToolbox/ColMosaic_Toolbox')
addpath('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg/SamplingToolbox/MonMosaic_Toolbox')
addpath('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg/SamplingToolbox/LowLevel_Toolbox')

% cones per mm^2
ecc = sqrt((-5)^2 + (-3.5)^2); % in deg
mmperdeg = .22;
conespermm2 = 150.9*10^3*exp(-1.2*ecc)+35.9*10^3*exp(-.16*ecc)+9.9*10^3*exp(-.03*ecc); % actually cones/mm^2 ~27,000
Sconespermm2 = 2.5*10^3 *exp(-.2*ecc)+1.8*10^3*exp(-.05*ecc);

gaborsddeg = 0.15; % Should be 0.4
gaborsdmm = gaborsddeg*mmperdeg;

% Ok, we're going to make a cone mosaic that extends "n" sds in both
% directions
mm2 = .6; % size of cone array in mm2
nstds_array = sqrt(mm2)/gaborsdmm; % Size of cone array re: Gabor SD
totconesinfield = round(mm2*conespermm2);
Sconesinfield = round(mm2*Sconespermm2);
Mconesinfield = round((totconesinfield-Sconesinfield)/2);

imsize = 2000; ratio = .8; jitter = 0.1; % paramaters taken from Curcio and Sloan paper
[Clist,nCones,nL,nM,nS] = MakeSloanClist(totconesinfield,Mconesinfield,Sconesinfield,imsize,ratio,jitter);
figure; axes; hold on;
L = Clist(:,3) == 1;
plot(Clist(L,1),Clist(L,2),'ro','MarkerFaceColor','red','MarkerSize',2); % 1.5 is good for sd = 0.4°
L = Clist(:,3) == 2;
plot(Clist(L,1),Clist(L,2),'go','MarkerFaceColor','green','MarkerSize',2);
L = Clist(:,3) == 3;
plot(Clist(L,1),Clist(L,2),'bo','MarkerFaceColor','blue','MarkerSize',2);
axis square; axis equal;


%%
% Section 2
% Spinning cone thresholds

% Looking at heat maps of threshold ratios using L:M ratios that are not
% 1:1. And the ideal observer weighting function stuff that he sent on
% 8/16/15.
% Need to load data_new_LtoM_1.mat, etc. Let's navigate to the appropriate
% directory.
MAKEMOVIE = 0;
cd ('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg');
load('data_new_LtoM_1.mat');
threshold_pts = gab.colorDirs .* repmat(cones.alpha_analytic,1,size(gab.colorDirs,2));
tmp = [threshold_pts;-threshold_pts];
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
evals = evals([3 2 1]);  % back in LMS order. Fragile
radii = sqrt(1./evals);
[lcc,mcc,scc] = meshgrid(linspace(-radii(1),radii(1),40),...
    linspace(-radii(2),radii(2),40),...
    linspace(-radii(3),radii(3),40));
fn = (lcc./radii(1)).^2+(mcc./radii(2)).^2+(scc./radii(3)).^2;

figure; axes; hold on;
set(gcf,'Color','none');
plot3(tmp(:,1),tmp(:,2),tmp(:,3),'k.');
p = patch(isosurface(lcc,mcc,scc,fn,1));
set(p,'EdgeColor', 'none', 'FaceAlpha',.8,'FaceColor',[0 .5 0],'Edgealpha',0,'SpecularExponent',.6);
set(p,'edgecolor','none','facecolor',[0 .6 0]);
lighting gouraud;
axis equal; 
set(gca,'View',[17 14]);
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
set(gca,'color','none','visible','off')
axis vis3d;
% Here are the radii for different L:M cone ratios:
% 1:1 .0022 .0023 .0079

moviefilename = 'SpinningConeIdlObsThresholds';
if (MAKEMOVIE)
    originalview = get(gca,'View');
    deltaviews = linspace(0,360,120);
    deltaviews(end) = [];
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[deltaviews(i) 0]);
        %set(gca,'xtick',[-.2 0 .2], 'ytick',[-.2 0 .2]);
        h_light = camlight(20, 20);
        drawnow;
        frame = getframe(gca);
        writeVideo(writerObj,frame);
        h_light.Visible = 'off';
    end
    close(writerObj);
end

%%
% Section 3) Monkey behavioral data
% These are DTNT data so they're not *directly* comparable to the
% physiology, but it's the best I can do for now.
load kali_DTNT_022815.mat
monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);
sfIdx = logical([1;0;0;0])  % which SF should be analyzed? [0.5 1 2 4]
sfs = unique(cat(1,dtnt.sfs{:}))

% make a logical vector to weed out conditions (if need be)
l_questBeta = dtnt.l_beta > 0;
l_bkgndrgb = dtnt.l_bkgndrgb > 0;

% pull out the spatial frequncy condition
l_sfs = ismember(cat(1,dtnt.sfs{:}),  sfs(sfIdx));

% grab the appropriate files
l_valid = [l_sfs & l_questBeta & l_bkgndrgb];
colors = cat(1, dtnt.colorDirs{l_valid});
alphas = cat(1, monkeyAlphas{l_valid});

% Trying to find thresholds for the four colors that Charlie investigated
% while recording from V1 neurons
Charliescolors = mkbasis([1 -1 0; 0 0 1; .1386 -.1386 .9806; -.1386 .1386 .9806]')';
Charliesalphas = [0 0 0 0]; % Kali's data, but Charlie's special color directions
for i = 1:size(Charliescolors,1)
    r = colors*Charliescolors(i,:)';
    L = r > .999;
    sum(L)
    Charliesalphas(i) = geomean(alphas(L));
end

figure; axes; hold on;
for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*Charliesalphas(i), ...
        Charliescolors(i,2)*Charliesalphas(i),...
        Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','black');
    plot3(-Charliescolors(i,1)*Charliesalphas(i), ...
        -Charliescolors(i,2)*Charliesalphas(i),...
        -Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','black');
end
axis equal;
set(gca,'Visible','off');
set(gca,'view',[-45 0]);

% S-cone threshold = 0.035
% L-M threshold = [0.005, -0.005]

%%
% Section 4
% Monkey neurophysiological (V1) data
load('/Volumes/NO BACKUP/NexFiles/Charlie/Batch Data And Text Files/CardVsInt_13-Oct-2012.mat')
n = size(out.dat,2);
intcolordirs = reshape([out.dat.prefInt],3,n)';
cardcolordirs = reshape([out.dat.prefCard],3,n)';
intTRs = [out.dat.intTR]';
cardTRs = [out.dat.cardTR]';
sfs = [];
for i = 1:length(out.dat)
    sfs(i) = out.dat(i).expt.sfs;
end
% Filtering out junk
L = ~isnan(intTRs) & ~isnan(cardTRs);
intcolordirs = intcolordirs(L,:);
cardcolordirs = cardcolordirs(L,:);
intTRs = intTRs(L);
cardTRs = cardTRs(L);
TRs = [cardTRs; intTRs];
colordirs = [cardcolordirs; intcolordirs];
Charliescolors = mkbasis([1 -1 0; 0 0 1; .1386 -.1386 .9806; -.1386 .1386 .9806]')';

geomeanTRs = [];
for i = 1:size(Charliescolors,1)
    r = Charliescolors(i,:)*colordirs';
    geomeanTRs(i) = geomean(TRs(abs(r) > .999 & sfs<=1.1));
end

% Picking up from the previous cell

figure; axes; hold on;
for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*Charliesalphas(i), ...
        Charliescolors(i,2)*Charliesalphas(i),...
        Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','black');
    plot3(-Charliescolors(i,1)*Charliesalphas(i), ...
        -Charliescolors(i,2)*Charliesalphas(i),...
        -Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','black');
end
for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*Charliesalphas(i)*geomeanTRs(i), ...
        Charliescolors(i,2)*Charliesalphas(i)*geomeanTRs(i),...
        Charliescolors(i,3)*Charliesalphas(i)*geomeanTRs(i),'go','MarkerFaceColor','green');
    plot3(-Charliescolors(i,1)*Charliesalphas(i)*geomeanTRs(i), ...
        -Charliescolors(i,2)*Charliesalphas(i)*geomeanTRs(i),...
        -Charliescolors(i,3)*Charliesalphas(i)*geomeanTRs(i),'go','MarkerFaceColor','green');
end

axis equal;
set(gca,'Visible','off');
set(gca,'view',[-45 0]);

%%
% Section 5
% An example Gabor
% (Go to iSETBIOFigs.m section 5

%%
% Section 6
% Plotting behavior, V1 neurons, and cones all on the same axes 
% Special "Charlie's colors" only

cd ('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg');
load('data_new_LtoM_1.mat');
threshold_pts = gab.colorDirs .* repmat(cones.alpha_analytic,1,size(gab.colorDirs,2));
tmp = [threshold_pts;-threshold_pts];
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
evals = evals([3 2 1]);  % back in LMS order. Fragile
radii = sqrt(1./evals);

figure; axes; hold on;
% Cones
for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*radii(1), ...
        Charliescolors(i,2)*radii(2),...
        Charliescolors(i,3)*radii(3),'ko','MarkerFaceColor','magenta');
    plot3(-Charliescolors(i,1)*radii(1), ...
        -Charliescolors(i,2)*radii(2),...
        -Charliescolors(i,3)*radii(3),'ko','MarkerFaceColor','magenta');
end

% Behavior
for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*Charliesalphas(i), ...
        Charliescolors(i,2)*Charliesalphas(i),...
        Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','yellow');
    plot3(-Charliescolors(i,1)*Charliesalphas(i), ...
        -Charliescolors(i,2)*Charliesalphas(i),...
        -Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','yellow');
end

% V1 neurons
for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*Charliesalphas(i)*geomeanTRs(i), ...
        Charliescolors(i,2)*Charliesalphas(i)*geomeanTRs(i),...
        Charliescolors(i,3)*Charliesalphas(i)*geomeanTRs(i),'go','MarkerFaceColor','green');
    plot3(-Charliescolors(i,1)*Charliesalphas(i)*geomeanTRs(i), ...
        -Charliescolors(i,2)*Charliesalphas(i)*geomeanTRs(i),...
        -Charliescolors(i,3)*Charliesalphas(i)*geomeanTRs(i),'go','MarkerFaceColor','green');
end
axis equal;

set(gca,'Visible','off');
set(gca,'view',[-45 0]);

%%
% Section 7
% Cone thresholds for the special stimuli plus fits
MAKEMOVIE = true;
PLOTTHRESHOLDPOINTS = false;
cd ('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg');
load('data_new_LtoM_1.mat');
threshold_pts = gab.colorDirs .* repmat(cones.alpha_analytic,1,size(gab.colorDirs,2));
tmp = [threshold_pts;-threshold_pts];
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
evals = evals([3 2 1]);  % back in LMS order. Fragile
radii = sqrt(1./evals);

[lcc,mcc,scc] = meshgrid(linspace(-radii(1),radii(1),40),...
    linspace(-radii(2),radii(2),40),...
    linspace(-radii(3),radii(3),40));
fn = (lcc./radii(1)).^2+(mcc./radii(2)).^2+(scc./radii(3)).^2;

figure; axes; hold on;
set(gcf,'Color','none');
if PLOTTHRESHOLDPOINTS
    plot3(tmp(:,1),tmp(:,2),tmp(:,3),'k.');
end
p = patch(isosurface(lcc,mcc,scc,fn,1));
set(p,'EdgeColor', 'none', 'FaceAlpha',1,'FaceColor','magenta','Edgealpha',0,'SpecularExponent',.6);
lighting gouraud;
axis equal; 
set(gca,'View',[-45 0]);
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
set(gca,'color','none','visible','off');
axis vis3d;
% Here are the radii for different L:M cone ratios:
% 1:1 .0022 .0023 .0079
if PLOTTHRESHOLDPOINTS
    for i = 1:size(Charliescolors,1)
        plot3(Charliescolors(i,1)*radii(1), ...
            Charliescolors(i,2)*radii(2),...
            Charliescolors(i,3)*radii(3),'mo','MarkerFaceColor','magenta','MarkerSize',10);
        plot3(-Charliescolors(i,1)*radii(1), ...
            -Charliescolors(i,2)*radii(2),...
            -Charliescolors(i,3)*radii(3),'mo','MarkerFaceColor','magenta','MarkerSize',10);
    end
end
%set(gca,'Zlim',[-.0085 .0085])
%set(gca,'XLimMode','manual','YLimMode','manual','ZLimMode','manual')
moviefilename = 'SpinningConeIdlObsThresholdsWithData';
if (MAKEMOVIE)
    originalview = get(gca,'View');
    deltaviews = linspace(0,360,120);
    deltaviews(end) = [];
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[deltaviews(i) 0]);
        h_light = camlight(20, 20);
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        h_light.Visible = 'off';
    end
    close(writerObj);
end

%%
% Section 8
% Behavioral thresholds for the special stimuli plus fits
MAKEMOVIE = 1;
PLOTSURFACE = 0;
% Plotting the surface fits
observers = {'kali'};         % Kali_DTNT_022815.mat or Sedna_DTNT_022815.mat
load modIsoSurf;
load ('kali_DTNT_022815.mat');
sfidx = 1; % 0.5 cpd
clear fpar_monkey
% divide the behavioral data by 100 so that % is b/w 0&1.
monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);

% which SF condition should be analyzed? The behvioral data has many SFs,
% but the model data is only at one SF.
sfs = unique(cat(1,dtnt.sfs{:}));

% make a logical vector to weed out conditions (if need be)
l_questBeta = dtnt.l_beta > 0;
l_bkgndrgb = dtnt.l_bkgndrgb > 0;

% pull out the spatial frequncy condition
l_sfs = cat(1,dtnt.sfs{:}) == sfs(sfidx);

% grab the appropriate files
l_valid = [l_sfs & l_questBeta & l_bkgndrgb];
colors = cat(1, dtnt.colorDirs{l_valid});
alphas = cat(1, monkeyAlphas{l_valid});

% default init params on the monkey
[fpar_monkey_default, fval_monkey_default] = fitDetectionSurface(colors, alphas, 'ellipsoid');

% init params from kali at 0.5 cpd (on the money)
initparams_kali = [1.15; -15.82; 38.29; -14.36; -17.29; 6.62; -8.74; -81.82; 68.40; 9.36];  % optional input argument for threshSurfPlot for 0.5cpd
[fpar_monkey, fval_monkey_kali] = fitDetectionSurface(colors, alphas, initparams_kali);

coordinates = bsxfun(@times, colors, alphas);

% Plotting the fit
figure('Position',[831   175   480   630]); 
axes; hold on;
plotlims = max(abs(coordinates));
[x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),50),linspace(-plotlims(2),plotlims(2),50),linspace(-plotlims(3),plotlims(3),50));
tmp = [x(:) y(:) z(:)];
v = sum(abs(tmp * reshape(fpar_monkey(2:end),3,3)).^fpar_monkey(1),2);
fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
if (PLOTSURFACE)
    p = patch(fv);
    set(p,'EdgeColor', 'none', 'FaceAlpha',.8,'FaceColor',[.5 .5 0],'Edgealpha',0,'SpecularExponent',.6);
end
lighting gouraud;
xlabel('L'); ylabel('M'); zlabel('S');
h_light = camlight(20,20);
plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'ko','MarkerFaceColor',[.4 .4 0],'MarkerEdgeColor','none','markersize', 4);
plot3(-coordinates(:,1),-coordinates(:,2),-coordinates(:,3),'ko','MarkerFaceColor',[.4 .4 0],'MarkerEdgeColor','none','markersize', 4);
set(gcf,'renderer','painters');
set(gca,'Visible','off');% for surface and points
set(gcf,'Color','none');
set(gca,'Color','none');
set(gca,'view',[-45 0]);
axis equal
axis vis3d

for i = 1:size(Charliescolors,1)
    plot3(Charliescolors(i,1)*Charliesalphas(i), ...
        Charliescolors(i,2)*Charliesalphas(i),...
        Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','yellow','MarkerSize',8);
    plot3(-Charliescolors(i,1)*Charliesalphas(i), ...
        -Charliescolors(i,2)*Charliesalphas(i),...
        -Charliescolors(i,3)*Charliesalphas(i),'ko','MarkerFaceColor','yellow','MarkerSize',8);
end

set(gca,'View',[-45 0]);
originalview = get(gca,'View');

moviefilename = 'SpinningPsychoThresholdsWithData';
if (MAKEMOVIE)
    originalview = get(gca,'View');
    deltaviews = linspace(0,360,200);
    deltaviews(end) = [];
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[deltaviews(i) 0]);
        h_light = camlight(20, 20);
        drawnow;
        frame = getframe(gca);
        writeVideo(writerObj,frame);
        h_light.Visible = 'off';
    end
    close(writerObj);
end
%%
% Section 9
% make a map of the ratio b/w Monkey and Retina. Where the ratio is large,
% the monkey does worse than the retina.
cd ('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg');
load modIsoSurf;
fpar_cones = fitDetectionSurface(modIsoSurf.colorDirs, modIsoSurf.thresholds, 'ellipsoid');

[X,Y,Z] = sphere(400);
spherecolors = [X(:), Y(:), Z(:)];
spherecolors = bsxfun(@rdivide, spherecolors, sqrt(sum(spherecolors.^2, 2)));
monkeyThresh = coleThresh(reshape(fpar_monkey(2:end),3,3), fpar_monkey(1), spherecolors);
conesThresh = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), spherecolors);
C = monkeyThresh ./ conesThresh;
C = log10(C);
C = reshape(C, size(X,1), size(X,2)); % for plotting
minVal = min(C(:))
minVec = [X(C(:) == minVal), Y(C(:) == minVal), Z(C(:) == minVal)] .* 1.65;
if size(minVec,1) ==1; minVec = [minVec; -minVec]; end
maxVal = max(C(:))
maxVec = [X(C(:) == maxVal), Y(C(:) == maxVal), Z(C(:) == maxVal)] .* 1.65;
if size(maxVec,1) ==1; maxVec = [maxVec; -maxVec]; end

figure; axes; hold on;
hand = surf(X,Y,Z,C);
set(hand, 'edgealpha', 0,'edgecolor','none');
plot3([-1.4, 1.4], [0 0], [0,0], '-', 'linewidth', 5, 'color', 'r') % the L direction
plot3([0 0], [-1.4, 1.4],  [0,0], '-', 'linewidth', 5, 'color', 'g') % the M direction
plot3([0, 0], [0 0], [-1.4 1.4], '-', 'linewidth', 5, 'color', 'b') % S-iso color direction
%plot3(minVec(:,1), minVec(:,2), minVec(:,3), '--', 'linewidth', 5, 'color', 'm')
xlim([-1.3, 1.3]); ylim([-1.3, 1.3]); zlim([-1.3, 1.3]);
set(gca, 'view', [-45 0],'color','none','visible', 'off')
axis square;
axis vis3d;
maps = {'Edge', 'CubicL', 'CubicYF', 'IsoL'};
colormap(pmkmp(250, maps{2}))
set(gcf,'renderer','zbuffer','color','none');

moviefilename = 'SpinningTRHeatMap';
if (MAKEMOVIE)
    originalview = get(gca,'View');
    deltaviews = linspace(0,1000,350);
    deltaviews(end) = [];
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:length(deltaviews)
        set(gca,'View',originalview+[deltaviews(i) 0]);
        %h_light = camlight(20, 20);
        drawnow;
        frame = getframe(gca);
        writeVideo(writerObj,frame);
        h_light.Visible = 'off';
    end
    close(writerObj);
end
figure;
colormap(pmkmp(250, maps{2}))
colorbar;
