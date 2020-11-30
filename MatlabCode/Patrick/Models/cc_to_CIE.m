% This script transforms Cone Contrast into CIE-xy
%global GLMP % this is just to load in the background

% Load variables from rawdata file
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    library = 'C:\Documents and Settings\JPatrickWeller\My Documents\Dropbox\Patrick\GLMS Data\nex files\';
end
datafile = 'N060613003.nex'; % +L-M
rawdata = nex2stro([char(library) datafile]);

% Load some variables from datafile
bkgndrgb = rawdata.sum.exptParams.bkgndrgb;
funds = reshape(rawdata.sum.exptParams.fundamentals,numel(rawdata.sum.exptParams.fundamentals)/3,3); % have to adjust the spacing in wls
fundswls = linspace(380,780,size(funds,1))';
monspds = reshape(rawdata.sum.exptParams.monspd,numel(rawdata.sum.exptParams.monspd)/3,3);
wls = linspace(380,780,size(monspds,1))';
fundamentals = SplineSpd(fundswls,funds,wls);
RGB2LMS = monspds' * fundamentals;

% Establish some coordinates in cone contrast space
LMweights = [10 -10];
lumcc = .9;
colcc = .09;
nRnds = 3;
thetaspace = pi/4;
rhospace = .5;
nconts = 20;
thetas = shiftdim(0:(thetaspace/nRnds):(2*pi-(thetaspace/nRnds)));
rhos = shiftdim(rhospace/nRnds:rhospace/nRnds:1);
FFIdx = fullfact([numel(thetas) numel(rhos)]);
LMlist = repmat([thetas(FFIdx(:,1)) rhos(FFIdx(:,2))],1,1);

% Transform from polar to cartesian coordinates
[Lcc,Mcc] = pol2cart(LMlist(:,1),LMlist(:,2));
    
% Scale Cone Contrast Units
scale = lumcc*colcc./sqrt((colcc.*cos(LMlist(:,1)-pi/4)).^2 ...
    +(lumcc.*sin(LMlist(:,1)-pi/4)).^2);
Lcc = Lcc .* scale;
Mcc = Mcc .* scale;
Scc = zeros(size(Lcc));
stimCC = cat(2,Lcc,Mcc,Scc);

% Generate mesh for rgb
x = linspace(0,1,100);
y = x;
[meshR,meshG] = meshgrid(x,y);
meshB = ones(size(meshR)) .* bkgndrgb(3);
meshrgb = cat(2,meshR(:),meshG(:),meshB(:));

% Convert mesh to cone exciation
bkgndCE = bkgndrgb' * RGB2LMS;
meshCE = meshrgb * RGB2LMS;
meshLCE = reshape(meshCE(:,1),size(meshR));
meshMCE = reshape(meshCE(:,2),size(meshR));

% To cone contrast
meshCD = meshCE - repmat(bkgndCE,size(meshCE,1),1);
meshCC = meshCD./repmat(bkgndCE,size(meshCE,1),1);
meshLCC = reshape(meshCC(:,1),size(meshR));
meshMCC = reshape(meshCC(:,2),size(meshR));
meshresps = max([meshLCC(:) meshMCC(:)] * LMweights',0);
meshresps = reshape(meshresps,size(meshR));
stimresps = max(([Lcc Mcc] * LMweights'),0);

% plot
figure(1); clf; hold on; grid on; box on; %axis equal;
plot3(Lcc,Mcc,stimresps,'ko')
contour(meshLCC,meshMCC,meshresps,nconts)
title('Cone Contrast Space')
xlabel('Lcc')
ylabel('Mcc')
zlabel('Response')


%% Convert from cone contrast space to cone excitation space

% Convert stimuli
bkgndCE =  bkgndrgb' * RGB2LMS;
stimCD = repmat(bkgndCE,numel(Lcc),1) .* stimCC;
stimCE = stimCD + repmat(bkgndCE,numel(Lcc),1);

% Plot
% figure(2); clf; hold on; grid on; %axis equal;
% plot3(stimCE(:,1),stimCE(:,2),stimresps,'go')
% contour(meshLCE,meshMCE,meshresps,nconts)
% title('Cone Excitation Space')
% xlabel('Le')
% ylabel('Me')
% zlabel('Response')


%% Converts from RGB space to CIE_xy space (using built-in function)

% Convert stim from cc to rgb to cie
stimrgb = stimCE * inv(RGB2LMS);
stimXYZ = rgb2xyz(stimrgb);
stimX = stimXYZ(:,1);
stimY = stimXYZ(:,2);
stimZ = stimXYZ(:,3);

% Convert mesh from cc to rgb to cie
meshXYZ = rgb2xyz(meshrgb);
meshX = meshXYZ(:,1);
meshY = meshXYZ(:,2);
meshZ = meshXYZ(:,3);

% Convert stim from CIE-XYZ to CIE-xy
stimx = stimX ./ (stimX + stimY + stimZ);
stimy = stimY ./ (stimX + stimY + stimZ);

% Convert mesh from CIE-XYZ to CIE-xy
meshx = reshape(meshX./sum(meshXYZ,2), size(meshR));
meshy = reshape(meshY./sum(meshXYZ,2), size(meshR));

% plot
figure(3); clf; hold on; grid on; box on; 
plot3(stimx,stimy,stimresps,'mo')
contour(meshx,meshy,meshresps,nconts)
title('RGB to xy (built in function)')
xlabel('x')
ylabel('y')
zlabel('Response')


%% Transform from LMS to CIE_XYZ

% Load in xyz color matching functions
load T_xyz1931
xyzvals = T_xyz1931;
xyzwls = S_xyz1931(1):S_xyz1931(2):(S_xyz1931(1) + S_xyz1931(2) * (S_xyz1931(3)-1));
xyzCMF = SplineSpd(xyzwls',xyzvals',wls);

% Transformation matrix from CE to XYZ
%LMS2XYZ = fundamentals' * xyzCMF;
LMS2XYZ = fundamentals\xyzCMF;

% Transform stim from ce to XYZ
stimXYZ = stimCE * LMS2XYZ;

% Transform mesh from ce to XYZ
meshXYZ = meshCE * LMS2XYZ;
meshX = meshXYZ(:,1);
meshY = meshXYZ(:,2);

% Transform stimuli from XYZ to xy space
stimx = stimXYZ(:,1) ./ sum(stimXYZ,2);
stimy = stimXYZ(:,2) ./ sum(stimXYZ,2);

% Transform mesh from XYZ to xy space
meshx = reshape(meshX./sum(meshXYZ,2), size(meshR));
meshy = reshape(meshY./sum(meshXYZ,2), size(meshR));

% Plot xy space
figure(4); clf; hold on; grid on; box on; 
plot3(stimx,stimy,stimresps,'mo')
contour(meshx,meshy,meshresps,nconts)
title('LMS to xy')
xlabel('x')
ylabel('y')
zlabel('Response')

%% Transform from RGB to XYZ

% Load in xyz color matching functions
load T_xyz1931
xyzvals = T_xyz1931;
xyzwls = S_xyz1931(1):S_xyz1931(2):(S_xyz1931(1) + S_xyz1931(2) * (S_xyz1931(3)-1));
xyzCMF = SplineSpd(xyzwls',xyzvals',wls);

% Transformation matrix from RGB to XYZ
RGB2XYZ = monspds' * xyzCMF;

% Transform stim from RGB to XYZ
stimXYZ = stimrgb * RGB2XYZ;

% Transform mesh from ce to XYZ
meshXYZ = meshrgb * RGB2XYZ;
meshX = meshXYZ(:,1);
meshY = meshXYZ(:,2);

% Transform stimuli from XYZ to xy space
stimx = stimXYZ(:,1) ./ sum(stimXYZ,2);
stimy = stimXYZ(:,2) ./ sum(stimXYZ,2);

% Transform mesh from XYZ to xy space
meshx = reshape(meshX./sum(meshXYZ,2), size(meshR));
meshy = reshape(meshY./sum(meshXYZ,2), size(meshR));

% Plot xy space
figure(5); clf; hold on; grid on; box on; 
plot3(stimx,stimy,stimresps,'mo')
contour(meshx,meshy,meshresps,nconts)
title('RGB to xy')
xlabel('x')
ylabel('y')
zlabel('Response')

%% LMS to XYZ by multiplying inv(RGB2LMS) * RGB2XYZ

% Transformation matrix (the long way)
LMS2XYZ = inv(monspds' * fundamentals) * (monspds' * xyzCMF)

% Transform stim from ce to XYZ
stimXYZ = stimCE * LMS2XYZ;

% Transform mesh from ce to XYZ
meshXYZ = meshCE * LMS2XYZ;
meshX = meshXYZ(:,1);
meshY = meshXYZ(:,2);

% Transform stimuli from XYZ to xy space
stimx = stimXYZ(:,1) ./ sum(stimXYZ,2);
stimy = stimXYZ(:,2) ./ sum(stimXYZ,2);

% Transform mesh from XYZ to xy space
meshx = reshape(meshX./sum(meshXYZ,2), size(meshR));
meshy = reshape(meshY./sum(meshXYZ,2), size(meshR));

% Plot xy space
figure(6); clf; hold on; grid on; box on; 
plot3(stimx,stimy,stimresps,'mo')
contour(meshx,meshy,meshresps,nconts)
title('rgb to xyz')
xlabel('x')
ylabel('y')
zlabel('Response')






