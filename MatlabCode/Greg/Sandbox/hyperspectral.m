% Playing around with hyperspectral images
% Need to load one of these files to get the variable "reflectances" in the
% workspace.
%
% Contents
%
% 0: just loading a file
% 1: rendering the scene under different illuminants
% 2: Gui for selecting points (and seeing how they change across
% illuminants
% 3: Comparing a "realistically" renders RGB image and one that constrains
% the cone exitation ratios to be fixed (?!))
%% Section 0
%load /Users/greghorwitz/Documents/MATLAB/scene4/ref_cyflower1bb_reg1.mat
load /Users/greghorwitz/Documents/MATLAB/scene7/ref_ribeira1bbb_reg1.mat

%% Section 1
% 10 nm bins fom 400 to 720 nm
% Reference: Foster, D.H., Nascimento, S.M.C., & Amano, K. (2004), 
% Information limits on neural identification of coloured surfaces in 
% natural scenes, Visual Neuroscience, 21, 331-336

wls = 400:10:720;
load T_xyz1964.mat
load T_cones_smj.mat
load B_cieday

% splining the fundamentals to match the hyperspectral data
fundamentals = SplineSpd([380:5:780]',T_cones_smj',wls')';

% Setting up illuminants
daylight_basisvectors = SplineSpd([380:5:780]',B_cieday,wls');
% Coefficients for the daylight spectra
% (from Judd et al. 1964 table 3)
daylight_coefficients = [-1.14 .677; 
    -0.784 -0.195;
    -0.293 -0.698;
    0.145 -0.752;
    1.005 -0.378];
nillums = 5;
coeffs = linspace(daylight_coefficients(1,1),daylight_coefficients(end,1),nillums);
daylight_coefficients_interp = [coeffs', interp1(daylight_coefficients(:,1),...
    daylight_coefficients(:,2),...
    coeffs,'spline')'];
illuminants = repmat(daylight_basisvectors(:,1),1,nillums)+daylight_basisvectors(:,[2 3])*daylight_coefficients_interp';

figure;
data =[];
for whichillum = 1:nillums
% plot(ref_n7(:,1),ref_n7(:,2));
illuminant = illuminants(:,whichillum);
%illuminant = flipud(illuminant); % Just for fun
%illuminant = illuminant(randperm(length(illuminant))); % Just for fun
%illuminant = zeros(length(illuminant),1); illuminant(round(whichillum/nillums*length(illuminant)))=1; % monochromatic lights
illuminant = normpdf(linspace(0,1,length(wls)),whichillum/nillums,.5)'; % narrow spectrum lights

lights = reflectances.*permute(repmat(illuminant,[1, size(reflectances,1),size(reflectances,2)]),[2 3 1]);
% lights is by far the largest variable in the workspace
lms = fundamentals*reshape(lights,size(lights,1)*size(lights,2),size(lights,3))';
lms = permute(reshape(lms,3,size(lights,1),size(lights,2)),[2 3 1]);
data = cat(4,data,lms);

% % Plotting cone excitations
% figure;
% for i = 1:3
%     subplot(2,2,i);
%     imagesc(lms(:,:,i));
%     colormap(gray)
% end
% subplot(2,2,4)
% imagesc(lms(:,:,1)./mean(mean(lms(:,:,1)))-lms(:,:,2)./mean(mean(lms(:,:,2))));

% Using RGBs from a random monitor for rendering
load('Dell4BitsCal')
rgbspect = SplineSpd([380:4:780]',cals{end}.P_device,wls');
M = fundamentals*rgbspect;

% Have to invert the gamma functions
gamma = cals{end}.gammaTable;
invgamma = InvertGamma(cals{end},0);

% Sythesizing an RGB image from the LMS image
RGBimage = zeros(size(lms));
n_by_3_lms = reshape(permute(lms,[3 1 2]),3,size(lms,1)*size(lms,2))';
three_by_n_rgb = inv(M)*n_by_3_lms';
three_by_n_rgb(three_by_n_rgb<0) = 0; % can't have values < 0
%three_by_n_rgb(three_by_n_rgb>prctile(three_by_n_rgb(:),99.99)) = 1; % HACK. Artificially clipping the highest intensities to brighten the image
three_by_n_rgb = three_by_n_rgb./max(three_by_n_rgb(:)); % can't have values > 1 to fill up the gamut
for i = 1:3
    three_by_n_rgb(i,:) = interp1(linspace(0,1,size(invgamma,1)),invgamma(:,i),three_by_n_rgb(i,:));
end

rgb_im = permute(reshape(three_by_n_rgb,3,size(lms,1),size(lms,2)),[2 3 1]);
%figure; 
%subplot(1,2,1);
%plot(rgb_im(:));
subplot(1,nillums,whichillum);
image(rgb_im);
axis square;
set(gca,'visible','off');
end

% -------------------
% Making an L-M cone contrast image
% normalizing to the overall cone excitations of the scene
% figure;
% for whichillum = 1:size(data,4)
%     Lcc = data(:,:,1,whichillum)./mean(mean(data(:,:,1,whichillum)));
%     Mcc = data(:,:,2,whichillum)./mean(mean(data(:,:,2,whichillum)));
%     Scc = data(:,:,3,whichillum)./mean(mean(data(:,:,3,whichillum)));
%     subplot(3,3,whichillum);
%     imagesc(Lcc-Mcc);
%    % imagesc(Lcc+Mcc);
%    % imagesc(Scc-(Lcc+Mcc)/2);
%    set(gca,'Visible','off');
%    axis square;
% end
% colormap(jet);

opp_cc = zeros(size(data,1),size(data,2),size(data,4));
for whichillum = 1:size(data,4)
    Lcc = (data(:,:,1,whichillum)-mean(mean(data(:,:,1,whichillum))))./mean(mean(data(:,:,1,whichillum)));
    Mcc = (data(:,:,2,whichillum)-mean(mean(data(:,:,2,whichillum))))./mean(mean(data(:,:,2,whichillum)));
    Scc = (data(:,:,3,whichillum)-mean(mean(data(:,:,3,whichillum))))./mean(mean(data(:,:,3,whichillum)));
    LM_opp_cc(:,:,whichillum) = Lcc-Mcc;
end

%Places in the scene where Lcc-Mcc changes a lot over the course of the day
figure;
tmp = LM_opp_cc(:,:,1)-LM_opp_cc(:,:,end); % difference of two zero mean vars is also zero mean.
image(255*(tmp.*1/(2*max(max(abs(tmp))))+.5));
colormap(jet(255))
title('changes in L-M over the course of the day');
% red means L-M is higher in the AM. 
% blue means L-M is higher in the PM
% (Remember, cone contrasts are normalized to the whole scene, so
% everything is getting bluer but something things are changing faster than
% others)

% Looking at trajectories in the Lcc-Mcc direction over the "course of the
% day"
figure; subplot(2,1,1); hold on; npix = 300;
tmp = [unidrnd(size(LM_opp_cc,1),npix,1),unidrnd(size(LM_opp_cc,2),npix,1)];
for i = 1:npix
    plot(LM_opp_cc(tmp(i,1),tmp(i,2),1),LM_opp_cc(tmp(i,1),tmp(i,2),end),'o');
end
xlabel('L-M under illuminant 1')
ylabel(['L-M under illuminant ',num2str(nillums)])
plot(0,0,'m*');
plot([-.2 .4],[-.2 .4],'k-')

% The changes are pretty modest

% Trying to make a DO filter
npix = 30;
kernel1 = normpdf(linspace(-3,3,npix))'*normpdf(linspace(-3,3,npix)).*repmat(sin(linspace(0,4*pi,npix)),npix,1);
kernel2 = normpdf(linspace(-3,3,npix))'*normpdf(linspace(-3,3,npix)).*repmat(sin(linspace(0,4*pi,npix))',1,npix);
figure;
subplot(2,2,1);
imagesc(conv2(LM_opp_cc(:,:,1),kernel1));
subplot(2,2,2);
imagesc(conv2(LM_opp_cc(:,:,1),kernel2));
subplot(2,2,3);
imagesc(conv2(LM_opp_cc(:,:,1),kernel1)-conv2(LM_opp_cc(:,:,end),kernel1));
title('Places where idealized DO response change the most between illuminants')
subplot(2,2,4);
imagesc(conv2(LM_opp_cc(:,:,1),kernel2)-conv2(LM_opp_cc(:,:,end),kernel2));


%%
% Making a gui for clicking on pixels and seeing how they evolve 
% in color space.
% Remember, "data" is in cone excitations

xyz = SplineSpd([380:5:780]',T_xyz1964',wls')';
M_lms_to_xyz = xyz*pinv(fundamentals);

lms = squeeze(data(:,:,:,4));
n_by_3_lms = reshape(permute(lms,[3 1 2]),3,size(lms,1)*size(lms,2))';
three_by_n_rgb = inv(M)*n_by_3_lms';
three_by_n_rgb(three_by_n_rgb<0) = 0; % can't have values < 0
three_by_n_rgb = three_by_n_rgb./max(three_by_n_rgb(:)); % can't have values > 1 to fill up the gamut
for i = 1:3
    three_by_n_rgb(i,:) = interp1(linspace(0,1,size(invgamma,1)),invgamma(:,i),three_by_n_rgb(i,:));
end

rgb_im = permute(reshape(three_by_n_rgb,3,size(lms,1),size(lms,2)),[2 3 1]);
h1 = figure;
image(rgb_im);
set(gca,'Visible','off'); axis square; hold on;
[x, y] = ginput;

x = round(x); y = round(y);
colors = jet(length(x));
figure(h1);
for i = 1:length(x)
    h = plot(x(i),y(i),'*');
    set(h,'color',colors(i,:))
end
% Now analyzing the selected points
mean_cone_excs = squeeze(mean(mean(data))); % mean across x and y (mean across whole image, for each cone type, for each illuminant)

h2 = figure;
subplot(2,2,1); hold on;
for i = 1:length(x) 
    lms = squeeze(data(y(i),x(i),:,:)); % intentionally flipping x and y here
    lms_ccs = (lms-mean_cone_excs)./mean_cone_excs;
    subplot(2,2,1); hold on;
    h = plot(lms_ccs(1,:)-lms_ccs(2,:),lms_ccs(3,:),'.-');
    set(h,'color',colors(i,:));
    subplot(2,2,2); hold on;
    XYZ = M_lms_to_xyz*lms;
    x_chrom = XYZ(1,:)./sum(XYZ);
    y_chrom = XYZ(2,:)./sum(XYZ);
    h = plot(x_chrom,y_chrom,'.-');
    set(h,'color',colors(i,:));
end
subplot(2,2,1); xlabel('L-M'); ylabel('S');
set(gca,'color',[.5 .5 .5])
subplot(2,2,2); xlabel('x'); ylabel('y');
set(gca,'color',[.5 .5 .5])

weights = [1 -1 0; .5 .5 -1; .5 -1 .5; 1 -.5 -.5];
if rem(length(x), 2) == 0 % Are there an even number of points? If so we assume each pair spans an edge
    figure;
    for i = 1:length(x)/2 % guaranteed to be an integer
        lms1 = squeeze(data(y(1+(i-1)*2),x(1+(i-1)*2),:,:));
        lms2 = squeeze(data(y(2+(i-1)*2),x(2+(i-1)*2),:,:));

        lms_ccs_1 = (lms1-mean_cone_excs)./mean_cone_excs; % This definition of cone contrast might have a big influence
        lms_ccs_2 = (lms2-mean_cone_excs)./mean_cone_excs;
        deltalms = lms1-lms2;
        for j = 1:size(weights,1)
            lincomb1 = weights(j,:)*lms_ccs_1; % weighted sum of cone signals (e.g. L-M) in pixel 1 as a function of illuminant
            lincomb2 = weights(j,:)*lms_ccs_2; % weighted sum of cone signals (e.g. L-M) in pixel 2 as a function of illuminant
            subplot(2,2,j); hold on;
            h = plot(lincomb1,lincomb2,'.-');
            set(h,'color',mean(colors([1 2]+(i-1)*2,:)))
            xlabel([num2str(weights(j,:)),' (cc) in pixel 1']);
            ylabel([num2str(weights(j,:)),' (cc) in pixel 2']);
            axis square;
            daspect([1 1 1]);
        end
    end
end

% Now looking at raw cone excitations
figure; axes; hold on;
%plot(lms1(1,:),lms2(1,:),'ro-'); % L-cone excitations
%plot(lms1(2,:),lms2(2,:),'go-'); % M-cone excitations

for i = 1:size(lms1,2) % looping over illuminants
    l = [lms1(1,i) lms2(1,i)];
    m = [lms1(2,i) lms2(2,i)];
    local_contrast_l = [l(1)/sum(l) l(2)/sum(l)];
    local_contrast_m = [m(1)/sum(m) m(2)/sum(m)];
    lum = local_contrast_l+local_contrast_m
    rg = local_contrast_l-local_contrast_m
    plot(lum(1),-lum(2),'ko');
    plot(rg(1),-rg(2),'ro');
end
xlabel('activation of subunit 1')
ylabel('activation of subunit 2')
% Under this model, whatever the L-M signal is in one subunit, it's always
% exactly the opposite in the other subunit. 

%%
% What does a straight line in delta-rgb space look like in (x,y)
deltargb = [.1 -.2 -.3]; % Just made this up
bkgndrgb = [.5 .5 .5];
load 
load T_xyz1964.mat
load('Dell4BitsCal')
rgbspect = SplineSpd([380:4:780]',cals{end}.P_device,[380:5:780]'); % to match XYZ
intensities = linspace(0,1,20);
rgbs = [];
XYZs = [];
for i = 1:length(intensities)
    rgbs(i,:) = intensities(i)*deltargb+bkgndrgb;
    XYZs(i,:) = T_xyz1964*rgbspect*rgbs(i,:)'
end

figure; axes; hold on;
plot(XYZs(:,1)./sum(XYZs,2),XYZs(:,2)./sum(XYZs,2));

% Now trying it with two equal-and-opposite subunits
rgbs = [];
XYZs = [];
for i = 1:length(intensities)
    rgbs(i,[1:3]) = intensities(i)*deltargb+bkgndrgb;
    rgbs(i,[4:6]) = intensities(end-i+1)*-deltargb+bkgndrgb;

    XYZs(i,[1:3]) = T_xyz1964*rgbspect*rgbs(i,[1:3])'
    XYZs(i,[4:6]) = T_xyz1964*rgbspect*rgbs(i,[4:6])' 
end
figure; axes; hold on; % RGB space
plot3(rgbs(:,1),rgbs(:,2),rgbs(:,3),'r-');
plot3(rgbs(:,4),rgbs(:,5),rgbs(:,6),'k-');

figure; axes; hold on; % xy space
plot(XYZs(:,1)./sum(XYZs(:,[1:3]),2),XYZs(:,2)./sum(XYZs(:,[1:3]),2),'ro');
plot(XYZs(:,4)./sum(XYZs(:,[4:6]),2),XYZs(:,5)./sum(XYZs(:,[4:6]),2),'ko');

%%
% Section 3
% Rendering a natural image under D65 illumination
% Need to have "reflectances" in the workspace already

wls = 400:10:720;
load T_xyz1964.mat
load T_cones_smj.mat
load spd_D65.mat
illuminant = SplineSpd([380:5:780]',spd_D65,wls');

lights = reflectances.*permute(repmat(illuminant,[1, size(reflectances,1),size(reflectances,2)]),[2 3 1]);
% lights is by far the largest variable in the workspace
lms = fundamentals*reshape(lights,size(lights,1)*size(lights,2),size(lights,3))';
lms = permute(reshape(lms,3,size(lights,1),size(lights,2)),[2 3 1]);
data = cat(4,data,lms);

% Settig up to make RGB image
load('Dell4BitsCal')
rgbspect = SplineSpd([380:4:780]',cals{end}.P_device,wls');
M = fundamentals*rgbspect;
gamma = cals{end}.gammaTable;
invgamma = InvertGamma(cals{end},0);

% Synthesizing an RGB image from the LMS image
RGBimage = zeros(size(lms));
n_by_3_lms = reshape(permute(lms,[3 1 2]),3,size(lms,1)*size(lms,2))';
three_by_n_rgb = inv(M)*n_by_3_lms';
three_by_n_rgb(three_by_n_rgb<0) = 0; % can't have values < 0
%three_by_n_rgb(three_by_n_rgb>prctile(three_by_n_rgb(:),99.99)) = 1; % HACK. Artificially clipping the highest intensities to brighten the image
three_by_n_rgb = three_by_n_rgb./max(three_by_n_rgb(:)); % can't have values > 1 to fill up the gamut
for i = 1:3
    three_by_n_rgb(i,:) = interp1(linspace(0,1,size(invgamma,1)),invgamma(:,i),three_by_n_rgb(i,:));
end

rgb_im = permute(reshape(three_by_n_rgb,3,size(lms,1),size(lms,2)),[2 3 1]);
