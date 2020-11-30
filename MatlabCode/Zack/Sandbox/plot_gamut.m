function plot_gamut(calfilename)
if ischar(calfilename) && length(calfilename) > 3 && ~strcmp(calfilename(end-3:end), '.mat')
    calfilename = [calfilename '.mat'];
end
if ~exist(calfilename, 'file')
    error('Cannot find the calibration file you specified');
end

% define the color directions.... just putting points on a sphere for now
nColors = 2000;
tmp = ceil(sqrt(nColors));
az = linspace(0, 2*pi, tmp);
el = linspace(0, pi/2, tmp);
inds = fullfact([tmp tmp]);
dirs_pol = [az(inds(:,1))', el(inds(:,2))'];
[x,y,z] = sph2cart(dirs_pol(:,1), dirs_pol(:,2), ones(size(dirs_pol,1),1));

colorDirs = [x y z];
norms = sqrt(sum(colorDirs.^2,2));
colorDirs = bsxfun(@rdivide, colorDirs, norms);
colorDirs(abs(colorDirs)<1e-14) = 0;
colorDirs = [colorDirs; 1 -1 0; -1 1 0; 0 0 -1; 0 0 1];
colorDirs = unique(colorDirs, 'rows'); %remove the duplicates

lvsm_idxs = abs(colorDirs(:,1)) == 1 & abs(colorDirs(:,2)) == 1 & colorDirs(:,3) == 0;
siso_idxs = colorDirs(:,1) == 0 & colorDirs(:,2) == 0 & abs(colorDirs(:,3)) == 1;

fprintf('loading calibration data from %s\n', which(calfilename));
calData = load(calfilename);
s = load('T_cones_smj10.mat');
fns = fieldnames(s);
fundamentals = s.(fns{1});
fundWavelengthSpacing = s.(fns{2});
calData = calData.cals{end};
bkgndRGB = round(255 .* calData.bgColor);
bkgndrgb = [calData.gammaTable(bkgndRGB(1)+1, 1), calData.gammaTable(bkgndRGB(2)+1, 2), ...
    calData.gammaTable(bkgndRGB(3)+1, 3)];

monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
M = fundamentals * monSpd;

[~,maxContrast] = gamutCheck(colorDirs',bkgndrgb,M,'both');

gamut_in_cc = bsxfun(@times,colorDirs,maxContrast');

fprintf('%s:\n', calData.describe.monitor);
fprintf('\tL-M gamut length: %g\n', sqrt(sum(diff(gamut_in_cc(lvsm_idxs,:)).^2)));
fprintf('\tS-iso gamut length: %g\n', sqrt(sum(diff(gamut_in_cc(siso_idxs,:)).^2)));

figure; hold on;
ourcolors = 'br';
for i = 0:1
    plot3((2*i-1)*gamut_in_cc(:,1),(2*i-1)*gamut_in_cc(:,2),(2*i-1)*gamut_in_cc(:,3),[ourcolors(i+1) '.']);
    plot3((2*i-1)*gamut_in_cc(:,1),(2*i-1)*gamut_in_cc(:,2),(2*i-1)*gamut_in_cc(:,3),[ourcolors(i+1) '.']);
end
title(sprintf('%s''s gamut in cc', calData.describe.monitor));
xlabel('L'); ylabel('M'); zlabel('S');
