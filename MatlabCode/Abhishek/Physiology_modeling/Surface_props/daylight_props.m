% Writing the new script for understanding how the perception of an edge
% changes with changes in the daylight spectra.
% Author - Abhishek De, 1/17
% Daylight spectra - sky.asc

% Code begins here
% Loading the illumination spectra
clearvars; close all;
load sky.asc
wave_il = 390:4:1070; % specific to this file
illuminant = [];
hi_wave_ind = find([390:5:1070]==780);
for ii = 1:size(sky,1)
    tmp = sky(ii,:);
    splinefit = spline(wave_il, tmp, [390:5:1070]);
    illuminant = [illuminant; splinefit(1:hi_wave_ind)];
end
% The illuminant has the illumination spectras ranging from 390 to 780 in steps of 5nm.
  
% Loading the reflectance spectra of two surfaces and place them adjacent to each other
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
L = numel(380:1:780); % for the munsell reflectance spectras - ftp://ftp.cs.joensuu.fi/pub/color/spectra/mspec/README.txt
ind = 1:5:L;
ind = ind(3:end);
rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
reflectance_spectra1 = munsell(ind,rand_idxs(1));
reflectance_spectra1 = reflectance_spectra1';
reflectance_spectra2 = munsell(ind,rand_idxs(2));
reflectance_spectra2 = reflectance_spectra2';

% Loading the cone action spectra and monitor spectral distributions
load fundamentals.mat
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
fundamentals = fundamentals(3:end,:);
mon_spd = mon_spd(:,3:end);
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations

hi = size(illuminant,1);
% Declaring the image plane
im_plane = ones(50,25);
L_cone_exc_ratios = zeros(1,hi);
M_cone_exc_ratios = zeros(1,hi);
S_cone_exc_ratios = zeros(1,hi);
RGBsurf1 = [];
RGBsurf2 = [];

wave = 390:5:780;
for ii = 1:1:hi
    % Calculating the net light energy entering the eye of an edge
    net_spectra1 = illuminant(ii,:).*reflectance_spectra1;
    net_spectra2 = illuminant(ii,:).*reflectance_spectra2;

    % Calculating the L,M,S excitations/quantal catches                                                   from the incident light
    L_exc1 = net_spectra1 * fundamentals(:,1);
    M_exc1 = net_spectra1 * fundamentals(:,2);
    S_exc1 = net_spectra1 * fundamentals(:,3);
    L_exc2 = net_spectra2 * fundamentals(:,1);
    M_exc2 = net_spectra2 * fundamentals(:,2);
    S_exc2 = net_spectra2 * fundamentals(:,3);
    
    % Calculating the RGB intensities from the L,M,S excitations after
    % normalizng the responses, 0 - don't normalize, 1 - do normalize
    [L_norm1,M_norm1,S_norm1,L_norm2,M_norm2,S_norm2] = normalize_resp(L_exc1,M_exc1,S_exc1,L_exc2,M_exc2,S_exc2,0);

    RGB1 = inv(M)*([L_norm1; M_norm1; S_norm1]);
    RGB_norm1 = (RGB1)/norm(RGB1);
    im1 = cat(3,im_plane*RGB_norm1(1),im_plane*RGB_norm1(2),im_plane*RGB_norm1(3));
    RGB2 = inv(M)*([L_norm2; M_norm2; S_norm2]);
    RGB_norm2 = (RGB2)/norm(RGB2);
    im2 = cat(3,im_plane*RGB_norm2(1),im_plane*RGB_norm2(2),im_plane*RGB_norm2(3));
    im = cat(2,im1,im2);
    RGBsurf1 = [RGBsurf1; RGB_norm1'];
    RGBsurf2 = [RGBsurf2; RGB_norm2'];
     
    % Plotting the illumination spectra 
    figure(1),subplot(hi/2,6,3*(ii-1)+1),plot(wave,illuminant(ii,:),'r','Linewidth',2); set(gca,'XTick',[],'YTick',[]);
    
    % Plotting the net incident energy reaching the eye -> illuminant spectra * reflectance spectra
    figure(1),subplot(hi/2,6,3*(ii-1)+2),plot(wave,net_spectra1,'k','Linewidth',2); hold on;
    plot(wave,net_spectra2,'g','Linewidth',2); set(gca,'XTick',[],'YTick',[]); hold off;
    
    % Plotting the color of the net spectra with respect to a gray background
    figure(1),subplot(hi/2,6,3*ii),image(im); set(gca,'XTick',[],'YTick',[],'Box','on');

    L_cone_exc_ratios(ii) = L_norm1/L_norm2;
    M_cone_exc_ratios(ii) = M_norm1/M_norm2;
    S_cone_exc_ratios(ii) = S_norm1/S_norm2;
end

figure(2);subplot(221),plot(wave,reflectance_spectra1,'g','Linewidth',2); hold on;
plot(wave,reflectance_spectra2,'m','Linewidth',2); title('Reflectance');
xlabel('wavelength'); ylabel('%'); hold off;

subplot(222),plot(L_cone_exc_ratios, M_cone_exc_ratios,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); 
xlabel('L cone ratios'); ylabel('M cone ratios');

subplot(223),plot(M_cone_exc_ratios, S_cone_exc_ratios,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); 
xlabel('M cone ratios'); ylabel('S cone ratios');

subplot(224),plot(L_cone_exc_ratios, S_cone_exc_ratios,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); 
xlabel('L cone ratios'); ylabel('S cone ratios');

% figure(3),plot3(L_cone_exc_ratios, M_cone_exc_ratios, S_cone_exc_ratios,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); 
% xlabel('L cone ratios'); ylabel('M cone ratios'); ylabel('S cone ratios');

%% Next step is to project the edge onto a the RF of an edge sensitive cell and see how the 
% what the isoresponse contour looks like 

% Create an STA of a hypothetical neuron, RF size L x W pixels
L = 50; W = 50;
mask = zeros(L,W);
mask(:,1:round(2*W/4)) = 1; % Subunit 1
mask(:,round(2*W/4)+1:end) = 2; % Subunit 2, rest all is background
bg = find(mask == 0);
S1 = find(mask == 1);
S2 = find(mask == 2);
RF = zeros(L*W,3);
fig_idx = 4;
for jj = fig_idx:fig_idx+2
    figure(jj)
    num_neurons = 5;
    for ii = 1:num_neurons
%         RGBS1  = [0.5+0.25*rand(1); 0.5+0.25*rand(1); 0.5+0.25*rand(1)]; RGBS2  = [0.5-0.25*rand(1); 0.5-0.25*rand(1); 0.5-0.25*rand(1)];
        RGBS1  = rand(3,1); RGBS2  = rand(3,1);

        RF(S1,1) = RGBS1(1); RF(S1,2) = RGBS1(2); RF(S1,3) = RGBS1(3);
        RF(S2,1) = RGBS2(1); RF(S2,2) = RGBS2(2); RF(S2,3) = RGBS2(3);
        Subunit_act1 = RGBsurf1 * RGBS1;
        Subunit_act2 = RGBsurf2 * RGBS2;
        subplot(num_neurons,2,2*(ii-1)+1); image(reshape(RF,[50 50 3])); set(gca,'XTick',[],'YTick',[]);
        subplot(num_neurons,2,2*ii), plot(Subunit_act1, Subunit_act2,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]);
        set(gca,'XTick',[],'YTick',[]); xlabel('S1'); ylabel('S2');
    end
end

