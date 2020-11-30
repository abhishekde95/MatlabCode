% Figures for Johns Hopkins 2017 talk
%
% Contents
% 1) Behavioral thresholds from Patrick's psychophysical data
% 1.1) As above but for more directions
% 2) Rotating principal axes of neurothresh ellipses
% (Single exmaple neurothresh ellipse, see SFN2010.m section 12)
% 3) Isodetection thresholds in 3-D cone space
% 4) Isodetection thresholds in 3-D cone space (with single mechanism
% plane superimposed)
% 5) Drifting Gabor movie
% 6) Still Gabors from example linear cell ('K082509007.nex') on iso.
% surface.
% 7) Bubble plots from Patrick's ephys data
% 8) Stimuli at each of the locations in the LM plane that Patrick tested
% 9) Patrick's white noise stimulus 
% -------------------------------------
%%
% Section 1) behavioral thresholds.
% Data file: N061513002
data = [   -0.0900    0.0000   60.0000   72.0000
   -0.0831   -0.0344   65.0000   76.0000
   -0.0831    0.0344   71.0000   78.0000
   -0.0788    0.0000   72.0000   79.0000
   -0.0728   -0.0301   65.0000   76.0000
   -0.0728    0.0301   74.0000   76.0000
   -0.0675    0.0000   66.0000   76.0000
   -0.0636   -0.0636   58.0000   76.0000
   -0.0636    0.0636   72.0000   78.0000
   -0.0624   -0.0258   66.0000   78.0000
   -0.0624    0.0258   71.0000   77.0000
   -0.0563    0.0000   68.0000   79.0000
   -0.0557   -0.0557   54.0000   77.0000
   -0.0557    0.0557   67.0000   77.0000
   -0.0520   -0.0215   67.0000   79.0000
   -0.0520    0.0215   71.0000   76.0000
   -0.0477   -0.0477   43.0000   78.0000
   -0.0477    0.0477   70.0000   78.0000
   -0.0450    0.0000   68.0000   75.0000
   -0.0416   -0.0172   54.0000   77.0000
   -0.0416    0.0172   67.0000   75.0000
   -0.0398   -0.0398   44.0000   78.0000
   -0.0398    0.0398   70.0000   75.0000
   -0.0344   -0.0831   73.0000   77.0000
   -0.0344    0.0831   74.0000   78.0000
   -0.0338    0.0000   52.0000   74.0000
   -0.0318   -0.0318   43.0000   78.0000
   -0.0318    0.0318   69.0000   78.0000
   -0.0312   -0.0129   52.0000   75.0000
   -0.0312    0.0129   64.0000   73.0000
   -0.0301   -0.0728   69.0000   75.0000
   -0.0301    0.0728   74.0000   78.0000
   -0.0258   -0.0624   67.0000   76.0000
   -0.0258    0.0624   70.0000   78.0000
   -0.0239   -0.0239   35.0000   77.0000
   -0.0239    0.0239   61.0000   77.0000
   -0.0225    0.0000   54.0000   75.0000
   -0.0215   -0.0520   67.0000   78.0000
   -0.0215    0.0520   72.0000   80.0000
   -0.0208   -0.0086   41.0000   76.0000
   -0.0208    0.0086   50.0000   77.0000
   -0.0172   -0.0416   57.0000   79.0000
   -0.0172    0.0416   65.0000   79.0000
   -0.0159   -0.0159   39.0000   74.0000
   -0.0159    0.0159   55.0000   79.0000
   -0.0129   -0.0312   49.0000   77.0000
   -0.0129    0.0312   57.0000   78.0000
   -0.0113    0.0000   45.0000   74.0000
   -0.0104   -0.0043   40.0000   77.0000
   -0.0104    0.0043   42.0000   77.0000
   -0.0086   -0.0208   45.0000   78.0000
   -0.0086    0.0208   53.0000   77.0000
   -0.0080   -0.0080   37.0000   78.0000
   -0.0080    0.0080   50.0000   79.0000
   -0.0043   -0.0104   41.0000   77.0000
   -0.0043    0.0104   42.0000   76.0000
   -0.0000   -0.0900   69.0000   73.0000
   -0.0000   -0.0788   74.0000   77.0000
   -0.0000   -0.0675   73.0000   76.0000
   -0.0000   -0.0563   71.0000   77.0000
   -0.0000   -0.0450   71.0000   78.0000
   -0.0000   -0.0338   67.0000   76.0000
   -0.0000   -0.0225   49.0000   78.0000
   -0.0000   -0.0113   44.0000   77.0000
    0.0000    0.0113   39.0000   76.0000
    0.0000    0.0225   52.0000   77.0000
    0.0000    0.0338   62.0000   76.0000
    0.0000    0.0450   63.0000   75.0000
    0.0000    0.0563   71.0000   75.0000
    0.0000    0.0675   69.0000   76.0000
    0.0000    0.0788   72.0000   76.0000
    0.0000    0.0900   67.0000   78.0000
    0.0043   -0.0104   40.0000   77.0000
    0.0043    0.0104   33.0000   76.0000
    0.0080   -0.0080   43.0000   77.0000
    0.0080    0.0080   41.0000   77.0000
    0.0086   -0.0208   61.0000   76.0000
    0.0086    0.0208   42.0000   77.0000
    0.0104   -0.0043   53.0000   79.0000
    0.0104    0.0043   39.0000   77.0000
    0.0113         0   42.0000   78.0000
    0.0129   -0.0312   69.0000   77.0000
    0.0129    0.0312   42.0000   77.0000
    0.0159   -0.0159   61.0000   76.0000
    0.0159    0.0159   40.0000   77.0000
    0.0172   -0.0416   69.0000   77.0000
    0.0172    0.0416   57.0000   78.0000
    0.0208   -0.0086   63.0000   78.0000
    0.0208    0.0086   45.0000   80.0000
    0.0215   -0.0520   73.0000   77.0000
    0.0215    0.0520   63.0000   77.0000
    0.0225         0   59.0000   78.0000
    0.0239   -0.0239   61.0000   76.0000
    0.0239    0.0239   45.0000   77.0000
    0.0258   -0.0624   67.0000   74.0000
    0.0258    0.0624   64.0000   77.0000
    0.0301   -0.0728   72.0000   78.0000
    0.0301    0.0728   66.0000   79.0000
    0.0312   -0.0129   69.0000   78.0000
    0.0312    0.0129   56.0000   78.0000
    0.0318   -0.0318   67.0000   76.0000
    0.0318    0.0318   46.0000   77.0000
    0.0338         0   65.0000   78.0000
    0.0344   -0.0831   77.0000   79.0000
    0.0344    0.0831   66.0000   77.0000
    0.0398   -0.0398   69.0000   76.0000
    0.0398    0.0398   57.0000   80.0000
    0.0416   -0.0172   70.0000   77.0000
    0.0416    0.0172   68.0000   76.0000
    0.0450         0   69.0000   76.0000
    0.0477   -0.0477   69.0000   76.0000
    0.0477    0.0477   57.0000   76.0000
    0.0520   -0.0215   71.0000   78.0000
    0.0520    0.0215   66.0000   77.0000
    0.0557   -0.0557   76.0000   79.0000
    0.0557    0.0557   63.0000   76.0000
    0.0563         0   72.0000   78.0000
    0.0624   -0.0258   71.0000   77.0000
    0.0624    0.0258   71.0000   77.0000
    0.0636   -0.0636   72.0000   76.0000
    0.0636    0.0636   65.0000   76.0000
    0.0675         0   72.0000   78.0000
    0.0728   -0.0301   70.0000   79.0000
    0.0728    0.0301   72.0000   77.0000
    0.0788         0   69.0000   77.0000
    0.0831   -0.0344   70.0000   76.0000
    0.0831    0.0344   68.0000   76.0000
    0.0900         0   71.0000   76.0000];

figure; axes; hold on;
for i = 1:size(data,1)
    h = plot3(data(i,1),data(i,2),data(i,3)./data(:,4),'ko','Markerfacecolor','black');
end

ccnorms = sqrt(sum(data(:,[1 2]).^2,2));
whichcomparison = [1 1; -1 -1; 1 0; 0 1];
threshold = [];
for i = 1:size(whichcomparison,1)
    unitvector = whichcomparison(i,:)./norm(whichcomparison(i,:));
    proj = data(:,[1 2])*unitvector'
    L = softEq(proj-ccnorms,0);
    [fittedparams, success] = weibullFitGH(ccnorms(L), data(L,[3,4]), 'mle',[0.03    1 nan nan]);
    threshold(i) = fittedparams(1);
end

figure; axes; hold on; 
plot(data(:,1),data(:,2),'.');
for i = 1:size(whichcomparison,1)
    unitvector = whichcomparison(i,:)./norm(whichcomparison(i,:));
    plot(threshold(i)*unitvector(1), threshold(i)*unitvector(2),'r*');
end
axis square

%%
% Section 1.1
% As above, but now getting detection thresholds in every direction
[t, ~] = cart2pol(data(:,1),data(:,2));
[a,b] = hist(t,400)
L = a > 0;
thetas = b(L);

for i = 1:length(thetas)
        unitvector = [cos(thetas(i)); sin(thetas(i))];
    proj = data(:,[1 2])*unitvector;
    L = softEq(proj-ccnorms,0,4);
    [fittedparams, success] = weibullFitGH(ccnorms(L), data(L,[3,4]), 'mle',[0.03    1 nan nan]);
    threshold(i) = fittedparams(1);
end

figure; axes; hold on;
[x,y] = pol2cart(thetas,threshold)
plot(x,y,'o')
axis square
set(gca,'Xlim',[-.08 .08],'ylim',[-.08 .08])


%%
% Section 2
% Principal axes
WHICHPCS = 1; % 1 or [1 2 3]
ELEVVIEWANGLE = 12;
STARTANGLE = -30;
if length(WHICHPCS) == 1
    TILTELEVATION = 0;
    ENDANGLE = STARTANGLE+360;
else
    TILTELEVATION = 1;
    ENDANGLE = 360;    
end
MAKEMOVIE = 1;
PLOTFRAME = 1;

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
eigenvectors = []; eigenvalues = [];
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
    
    out = NTpreprocess(NT,0,Inf);  % Getting the termination points
    scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    Loog = logical(out(:,end));
 
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    A = [quadparams(1) quadparams(4) quadparams(5);...
         quadparams(4) quadparams(2) quadparams(6);...
         quadparams(5) quadparams(6) quadparams(3)];
    % Transforming quadparams into un-whitened space
    [evecs,evals] = eig(xformmat*A*xformmat');
    [evals,i] = sort(diag(evals),1,'ascend');
    evecs = evecs(:,i);
    
    eigenvectors = cat(3,eigenvectors,evecs);
    eigenvalues = cat(2,eigenvalues,evals);
end

% Principal axes for ellipsoidal cells
evalratiothresh = 8;
out = [];
for i = 1:size(eigenvalues,2)
    out = [out; reshape(eigenvectors(:,:,i),1,9)];
    ratio12 = max(abs(sqrt(eigenvalues(1,i))./sqrt(eigenvalues(2,i))),sqrt(abs(eigenvalues(2,i))./sqrt(eigenvalues(1,i))));
    ratio23 = max(abs(sqrt(eigenvalues(2,i))./sqrt(eigenvalues(3,i))),sqrt(abs(eigenvalues(3,i))./sqrt(eigenvalues(2,i))));
    if (ratio12 < evalratiothresh)
        out(end,[1:6]) = nan;
    end
    if (ratio23 < evalratiothresh)
        out(end,[4:9]) = nan;
    end
end
Lellipse = all(eigenvalues>0)';
Ltiedaxes = any(isnan(out),2);
L = Lellipse&~Ltiedaxes;
L = Lellipse;
out = out(L,:);

figure; axes;
set(gcf,'Position',[250   238   990   860]); % Nice and big
set(gca,'Position',[0.05 0.05 .9 .9]); % fill the figure
set(gca,'Color','none');
set(gcf,'Color','none');
set(gca,'Visible','off');
set(gca,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set(gca,'View',[44 47],'Units','pixels') % To make sure everything stays in the frame
axis vis3d;
hold on;
colors = [.5 .5 .5; 0 0 1; 1 1 0];

for j = 1:size(out,1)
    for i = WHICHPCS
        idxs = [1 2 3]+(i-1)*3;
        h = plot3(out(j,idxs(1))*[-1 1], out(j,idxs(2))*[-1 1], out(j,idxs(3))*[-1 1],'-','linewidth',2);
        set(h,'Color',colors(i,:));
    end
end
if (PLOTFRAME)
    plot3([-1 1 1 -1 -1],[-1 -1 1 1 -1],[-1 -1 -1 -1 -1],'w-'); % face 1
    plot3([-1 1 1 -1 -1],[-1 -1 1 1 -1],[1 1 1 1 1],'w-'); % face 2
    plot3([-1 -1],[-1 -1],[-1 1],'w-'); % edge
    plot3([-1 -1],[-1 -1],[-1 1],'w-'); % edge
    plot3([-1 -1],[1 1],[-1 1],'w-'); % edge
    plot3([1 1],[1 1],[-1 1],'w-'); % edge
    plot3([1 1],[-1 -1],[-1 1],'w-'); % edge
end

if (MAKEMOVIE)
    clear M;
    viewangles = [STARTANGLE:3:ENDANGLE];
    viewangles([end+1 end+2]) = viewangles(end)+[3 6]; % having to tack on extra frames for Keynote, not sure why this is neccesary
    if (TILTELEVATION)
        elevangledeltas = linspace(0, 90-ELEVVIEWANGLE, length(viewangles-2));
        elevangledeltas([end+1 end+2]) = 90-ELEVVIEWANGLE;
    else
        elevangledeltas = zeros(length(viewangles),1);
    end
    for i = 1:length(viewangles)
        set(gca,'View',[viewangles(i) ELEVVIEWANGLE+elevangledeltas(i)])
        axespos = get(gca,'position');
        M(i) = getframe(gcf,axespos);% 850 750 works
    end
    repeat = 1;     %default = 1
    pSearch = 1;    %default = 0
    bSearch = 1;    %default = 1
    reference = 1;  %default = 0
    pixRange = 10;  %default = 10
    iFrame = 8;     %default = 8
    pFrame = 10;    %default = 10
    bFrame = 25;    %default = 25
    options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
    mpgwrite(M, gray, 'rotdata3.mpg', options);
end

%%
% Section 3
% Isodetection thresholds in 3-D space
% Based on PBIO2016.m section 8
% Path: /Users/greghorwitz/Documents/Manuscripts/Completed/Charlie's model/For Greg
matfilepath = '/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg';
MAKEMOVIE = 0;
PLOTSURFACE = 0;
PLOTFRAME = 1;
FRAMECOLOR = 'white';
% Plotting the surface fits
load ([matfilepath,filesep,'kali_DTNT_022815.mat']);
sfidx = 3; % 0.5, 1, 2, 4 cpd
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
% fpar_monkey. First parameter is beta, after that, each row of three is a
% mechanism.

% Getting raw data points ------------
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
coordinates = bsxfun(@times, colors, alphas);
% ----------------------------- 

% Plotting the fit
figure('Position',[831   175   480   630]); 
set(gcf,'InvertHardCopy','off');
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
set(gca,'view',[-45 0]); % Getting the initial axes size right
axis equal
axis vis3d

if (PLOTFRAME)
    
    edgedist = .1;
    h = plot3(edgedist*[-1 1 1 -1 -1],edgedist*[-1 -1 1 1 -1],edgedist*[-1 -1 -1 -1 -1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % face 1
    plot3(edgedist*[-1 1 1 -1 -1],edgedist*[-1 -1 1 1 -1],edgedist*[1 1 1 1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % face 2
    plot3(edgedist*[-1 -1],edgedist*[-1 -1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[-1 -1],edgedist*[-1 -1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[-1 -1],edgedist*[1 1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[1 1],edgedist*[1 1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[1 1],edgedist*[-1 -1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
end

set(gca,'View',[-55 22]);
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
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        h_light.Visible = 'off';
    end
    close(writerObj);
end

%%
% Section 4
% Plane showing one of the fitted detection mechanisms
matfilepath = '/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg';
MAKEMOVIE = 0;
PLOTFRAME = 1;
FRAMECOLOR = 'white';
WHICHMECHANISMS = [1 2 3];
% Plotting the surface fits
load ([matfilepath,filesep,'kali_DTNT_022815.mat']);
sfidx = 3; % 0.5, 1, 2, 4 cpd
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
% fpar_monkey. First parameter is beta, after that, each row of three is a
% mechanism.

coordinates = bsxfun(@times, colors, alphas);
[x,y,z] = meshgrid(linspace(-.1,.1,2),linspace(-.1,.1,2),linspace(-.1,.1,2));
mechanisms = reshape(fpar_monkey(2:end),3,3);
mechanism_activations = [x(:), y(:), z(:)]*mechanisms

% Plotting
figure('Position',[831   175   480   630]); 
axes; hold on;
lighting gouraud;
xlabel('L'); ylabel('M'); zlabel('S');
h_light = camlight(20,20);
plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'ko','MarkerFaceColor',[.4 .4 0],'MarkerEdgeColor','none','markersize', 4);
plot3(-coordinates(:,1),-coordinates(:,2),-coordinates(:,3),'ko','MarkerFaceColor',[.4 .4 0],'MarkerEdgeColor','none','markersize', 4);
set(gcf,'renderer','painters');
set(gca,'Visible','off');% for surface and points
set(gcf,'Color','none');
set(gca,'Color','none');
set(gca,'view',[-45 0]); % Getting the initial axes size right
axis equal
axis vis3d

for i = 1:length(WHICHMECHANISMS)
    patch_info = isosurface(x,y,z,reshape(mechanism_activations(:,WHICHMECHANISMS(i)),size(x)),1)
    patch(patch_info,'FaceColor',[0 1 0],'FaceAlpha',.5,'EdgeColor','none')
    patch_info.vertices = -1*patch_info.vertices;
    patch(patch_info,'FaceColor',[0 1 0],'FaceAlpha',.5,'EdgeColor','none')
end
if PLOTFRAME
    edgedist = .1;
    h = plot3(edgedist*[-1 1 1 -1 -1],edgedist*[-1 -1 1 1 -1],edgedist*[-1 -1 -1 -1 -1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % face 1
    plot3(edgedist*[-1 1 1 -1 -1],edgedist*[-1 -1 1 1 -1],edgedist*[1 1 1 1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % face 2
    plot3(edgedist*[-1 -1],edgedist*[-1 -1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[-1 -1],edgedist*[-1 -1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[-1 -1],edgedist*[1 1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[1 1],edgedist*[1 1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
    plot3(edgedist*[1 1],edgedist*[-1 -1],edgedist*[-1 1],'-','LineWidth',1.5,'Color',FRAMECOLOR); % edge
end
h_light = camlight(20,20);


set(gca,'View',[-55 22]);
originalview = get(gca,'View');

moviefilename = 'SpinningPlaneWithData';
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
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        h_light.Visible = 'off';
    end
    close(writerObj);
end

set(gcf,'InvertHardCopy','off','Renderer','OpenGL')
%%
% Section 5
% Drifting Gabor movie
% M matrix taken from K082509007.nex
% Isostim are stimuli that stimulate 'K082509007.nex' at criterion fr.

isostim = [-0.3398 -.06525 0.3947; 
            0.7957 0.8724 0.6506;
            0.9611 0.9611 0;
            0.42821 0.5991 -0.6407;
            0 0 -1.139;
            -0.8014 -0.9864 -0.9209];

M = [   0.0608    0.1219    0.0175
    0.0220    0.1266    0.0257
    0.0019    0.0095    0.0976];

MAKEMOVIE = 1;
MAKESTILLS = 0;
bkgndrgb = [.5 .5 .5];
gaborrgb = inv(M)*[.5 .5 0]';
gaborrgb = gaborrgb./(2*max(gaborrgb));
gaborrgb = gaborrgb
deltaphi = 1;
edgergb = [0 0 0];
theta = 0;
lambda = 2;
sigma = .5;
gamma = 1;
phi = 0;
xoff = 0;
yoff = 0;
etheta = 0;
edisp = 0;
gausslim = .999;
pixperdeg = 20;
nframes = 150;
phis = phi-[1:nframes]*deltaphi;
ramp = linspace(0,1,nframes/4);
temporalenvelope = [ramp ones(1,nframes-2*length(ramp)), fliplr(ramp)];

im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg);
if MAKESTILLS
    figure;
    set(gcf,'Color',[.5 .5 .5])
    subplot(2,1,1);
    axis image
    set(gca,'Visible','off','Xtick',[],'Ytick',[])
    im = DrawGaborEdge(bkgndrgb, temporalenvelope(round(nframes/2))*gaborrgb, edgergb, theta, lambda, sigma, gamma, phis(1), xoff, yoff, etheta, edisp, gausslim, pixperdeg);
    image(im);
    set(gca,'Visible','off')
    axis image
    subplot(2,1,2)
    axis image
    im = DrawGaborEdge(bkgndrgb, temporalenvelope(round(nframes/2))*gaborrgb, edgergb, theta, lambda, sigma, gamma, phis(1)+pi, xoff, yoff, etheta, edisp, gausslim, pixperdeg);
    image(im);
    set(gca,'Visible','off')
    axis image
end

% make a movie
if MAKEMOVIE
    image(im)
    axis image
    set(gca,'Visible','off','Xtick',[],'Ytick',[])
    set(gcf,'Color',[.5 .5 .5])
    
    moviefilename = 'GaborMovie1';
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:nframes
        im = DrawGaborEdge(bkgndrgb, temporalenvelope(i)*gaborrgb, edgergb, theta, lambda, sigma, gamma, phis(i), xoff, yoff, etheta, edisp, gausslim, pixperdeg);
        image(im);
        set(gca,'Visible','off')
        axis image
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    close(writerObj);
end


%%
% Section 6
% Isostim are stimuli that stimulate 'K082509007.nex' at criterion fr.

isostim = [-0.3398 -.06525 0.3947; 
            0.7957 0.8724 0.6506;
            0.9611 0.9611 0;
            0.42821 0.5991 -0.6407;
            0 0 -1.139;
            -0.8014 -0.9864 -0.9209];

M = [   0.0608    0.1219    0.0175
    0.0220    0.1266    0.0257
    0.0019    0.0095    0.0976];
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';

edgergb = [0 0 0];
theta = 0;
lambda = 2;
sigma = .5;
gamma = 1;
phi = 0;
xoff = 0;
yoff = 0;
etheta = 0;
edisp = 0;
gausslim = .999;
pixperdeg = 20;

figure;
for i = 1:size(isostim,1)
    gaborrgb = inv(M)*(bkgndlms.*(1+isostim(i,:))') - bkgndrgb;
    im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg);
    subplot(ceil(sqrt(size(isostim,1))),ceil(sqrt(size(isostim,1))),i)
    image(im);
    set(gca,'Visible','off')
    axis image
end
%%
% Section 7 GLMSubunit data
LUMPROJTHRESH = .7; % threshold for projection onto the L+M axis
PLOTCONTOURS = 0;
patrickpath = '~/Documents/MATLAB/DataFromPatrick';
filename = 'N042615001_pancolor';
filename = 'N011215002_lumHS'; % asymmetric pancolor
%filename = 'N120613002'; % Chrom horseshoe
filename = 'N010515003_mLpM_HS'
%filename = 'N112614001'; % asymmetric pancolor
%filename = 'N111414001'; % bipolar luminance
%filename = 'N020215001' % sparse bipolar luminance
%filename = 'N012115001' % M-L
%filename = 'N070313003' % M-L (circular stim distn)
%filename = 'M080315001'; % L-M
%filename = 'N022715001' % ON luminance
%filename = 'N022315004' % OFF luminance
filename = 'N052915001s1' % ON luminance
load([patrickpath,filesep,filename]);
data = eval(filename);
uniqueLMs = unique(data(:,[1 2]),'rows');
nreps = zeros(size(uniqueLMs,1),1);
for i = 1:length(uniqueLMs)
    nreps(i) = sum((data(:,1) == uniqueLMs(i,1)) & (data(:,2) == uniqueLMs(i,2)));
end

% weeding out non-repeated stimuli
uniqueLMs(nreps<2,:) = [];
% weeding out very high luminance stimuli
proj = uniqueLMs*[1/sqrt(2); 1/sqrt(2)];
uniqueLMs(abs(proj) > LUMPROJTHRESH,:) = [];

figure; axes('Units','inches','Position',[.5 .5 5 5]); hold on;
set(gca,'Color','black');
frbounds = [0 inf]; % (max and min mean firing rate across stimuli)
markersizebounds = [0 0];
SCALEFACTOR = max(data(:,3));
mn = [];
for i = 1:size(uniqueLMs,1)
    L = data(:,1) == uniqueLMs(i,1) &  data(:,2) == uniqueLMs(i,2);
    mn(i) = mean(data(L,3));
    sd = std(data(L,3));
    n = sum(L);
    h = plot(uniqueLMs(i,1), uniqueLMs(i,2),'ko','Markerfacecolor','yellow');
    set(h,'Markersize',50*mn(i)/(SCALEFACTOR)+2);
    if (mn(i) > frbounds(1))
        frbounds(1) = mn(i);
        markersizebounds(1) = get(h,'MarkerSize');
    end
    if (mn(i) < frbounds(2))
        frbounds(2) = mn(i);
        markersizebounds(2) = get(h,'MarkerSize');
    end
end

if (PLOTCONTOURS)
    f = tpaps(uniqueLMs',mn,.95)
    [v,d] = eig(cov(uniqueLMs));
    [x,y] = meshgrid([d(1,1)*[-30:1:30]],[d(2,2)*[-30:1:30]]);
    rotpoints = [x(:),y(:)]*v;
    F = fnval(f,rotpoints');
    [~,h] = contour(reshape(rotpoints(:,1),size(x)),reshape(rotpoints(:,2),size(x)),reshape(F,size(x)),[2 4 6 8])
    set(h,'Linewidth',2,'LineColor','white')
end

% Legend bar
nsymbols = 4;
axes('Units','inches','Position',[6 .5 1 5]); hold on;
set(gca,'Color','black');
markersizes = linspace(markersizebounds(1), markersizebounds(2), nsymbols);
ypos = linspace(.9,.2, nsymbols);
frs = linspace(frbounds(1),frbounds(2),nsymbols);
for i = 1:nsymbols
    plot(.5,ypos(i),'yo','MarkerFaceColor','yellow','MarkerSize',markersizes(i));
    text(2,ypos(i),num2str(frs(i)*5,'%1.0f')); % *5 to get to spikes/sec (200 ms counting window)
end
set(gca,'Ylim',[0 1],'XTick',[],'YTick',[]);
set(gcf,'InvertHardcopy','off')


%%
% Section 8
% Stimuli in the locations Patrick tested for one cell (in their
% appropriate colors).
M = [   0.0608    0.1219    0.0175
    0.0220    0.1266    0.0257
    0.0019    0.0095    0.0976]; % From Section 5, above
LUMPROJTHRESH = .7; % threshold for projection onto the L+M axis
PLOTCONTOURS = 1;
patrickpath = '~/Documents/MATLAB/DataFromPatrick';
filename = 'N042615001_pancolor';
filename = 'N012115001'
%filename = 'N112614001'; % asymmetric pancolor
load([patrickpath,filesep,filename]);
data = eval(filename);
uniqueLMs = unique(data(:,[1 2]),'rows');
nreps = zeros(size(uniqueLMs,1),1);
for i = 1:length(uniqueLMs)
    nreps(i) = sum((data(:,1) == uniqueLMs(i,1)) & (data(:,2) == uniqueLMs(i,2)));
end

% weeding out non-repeated stimuli
uniqueLMs(nreps<2,:) = [];
% weeding out very high luminance stimuliM
proj = uniqueLMs*[1/sqrt(2); 1/sqrt(2)];
uniqueLMs(abs(proj) > LUMPROJTHRESH,:) = [];

figure; axes('Units','inches','Position',[.5 .5 5 5]); hold on;
set(gca,'Color',[.5 .5 .5],'Xtick',[],'Ytick',[]);
rgbs = [];
bkgndlms = M*[.5 .5 .5]';
for i = 1:size(uniqueLMs,1)
   cone_exc = bkgndlms.*(1+[uniqueLMs(i,:) 0]');
   rgbs(i,:) = inv(M)*cone_exc
end
   
for i = 1:size(uniqueLMs,1)
    h = plot(uniqueLMs(i,1),uniqueLMs(i,2),'ks','MarkerSize',10);
    set(h,'MarkerFaceColor',rgbs(i,:));
end

%%
% Section 9 Patrick's white noise stimulus

uniqueLMs = [.7 .7; .09  -.09; -.7 -.7; -.09 .09];

M = [   0.0608    0.1219    0.0175
    0.0220    0.1266    0.0257
    0.0019    0.0095    0.0976]; % From Section 5, above
rgbs = [];
bkgndlms = M*[.5 .5 .5]';
for i = 1:size(uniqueLMs,1)
    cone_exc = bkgndlms.*(1+[uniqueLMs(i,:) 0]');
    rgbs(i,:) = inv(M)*cone_exc
end


moviefilename = 'Patrick_stimulus';
figure; axes('units','normalized','position',[0 0 .7 1]);
nframes = 100;
colormap(rgbs);

mat = unidrnd(size(rgbs,1),10, 10);
image(mat)
set(gca,'Xtick',[],'Ytick',[],'Ztick',[],'Box','off');
drawnow;
set(gca,'Projection','Perspective')
set(gca,'View',[-50 40],'Color',[.5 .5 .5])
if (MAKEMOVIE)
    writerObj = VideoWriter(moviefilename);
    open(writerObj);
    for i = 1:nframes
        mat = unidrnd(size(rgbs,1),10, 10);
        image(mat)
        set(gca,'View',[-76.8000   68.8000],'Visible','off')
        set(gca,'Color',[.5 .5 .5])
        set(gcf,'Color',[.5 .5 .5])
        set(gca,'Xtick',[],'Ytick',[],'Ztick',[],'Box','off');
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    close(writerObj);
end




