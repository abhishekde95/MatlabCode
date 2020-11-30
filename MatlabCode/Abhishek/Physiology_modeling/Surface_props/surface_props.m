% Script to understand how different illuminant spectra change the net
% light entering the eye and what computations are needed by neuron in
% order to achieve an invariant representation or a surface or an edge
% Author - Abhishek De, 1/17

% Some spectra from the Psychtoolbox 
% spd_apparatusrel.mat spd_CIEA.mat spd_CIEC.mat spd_D65.mat spd_flourescent.mat
% spd_incanCC.mat spd_phillybright.mat spd_xenonArc.mat spd_xenonFlash.mat

% ASCII files, Load them as - load filename.asc 
% Reference white - baso4.asc
% Daylight spectra - aky.asc , illumination spectra, I guess
% Tree spectra - tree.asc
% Munsell reflectance spectras - munsell380_800_1.mat

% Code begins here
% Loading the illumination spectra
clearvars; close all;
load spd_CIEA.mat;
load spd_D65.mat
load spd_incanCC.mat
load spd_flourescent.mat
wave = 380:5:780;
illuminant = [spd_CIEA spd_D65 spd_incanCC spd_flourescent];

% Loading the reflectance spectra
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
L = numel(380:1:780); % for the munsell reflectance spectras - ftp://ftp.cs.joensuu.fi/pub/color/spectra/mspec/README.txt
ind = 1:5:L;
rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
reflectance_spectra1 = munsell(ind,rand_idxs(1));
reflectance_spectra2 = munsell(ind,rand_idxs(2));

% Loading the cone action spectra and monito spectral distributions
load fundamentals.mat
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations

% Declaring the image plane
im_plane = ones(50,25);
L_cone_exc_ratios = zeros(1,4);
M_cone_exc_ratios = zeros(1,4);
S_cone_exc_ratios = zeros(1,4);

for ii = 1:4
    % Calculating the net light energy entering the eye of an edge
    net_spectra1 = illuminant(:,ii).*reflectance_spectra1;
    net_spectra2 = illuminant(:,ii).*reflectance_spectra2;

    % Calculating the L,M,S excitations from the incident light
    L_exc1 = net_spectra1' * fundamentals(:,1);
    M_exc1 = net_spectra1' * fundamentals(:,2);
    S_exc1 = net_spectra1' * fundamentals(:,3);
    L_exc2 = net_spectra2' * fundamentals(:,1);
    M_exc2 = net_spectra2' * fundamentals(:,2);
    S_exc2 = net_spectra2' * fundamentals(:,3);
    
    % Calculating the RGB intensities  from the L,M,S excitations
    L_norm1 = L_exc1/norm([L_exc1;M_exc1;S_exc1]);
    M_norm1 = M_exc1/norm([L_exc1;M_exc1;S_exc1]);
    S_norm1 = S_exc1/norm([L_exc1;M_exc1;S_exc1]);
    L_norm2 = L_exc2/norm([L_exc2;M_exc2;S_exc2]);
    M_norm2 = M_exc2/norm([L_exc2;M_exc2;S_exc2]);
    S_norm2 = S_exc2/norm([L_exc2;M_exc2;S_exc2]);
   
    RGB1 = inv(M)*([L_norm1; M_norm1; S_norm1]);
    RGB_norm1 = (0.5*RGB1)/norm(RGB1);
    im1 = cat(3,im_plane*RGB_norm1(1),im_plane*RGB_norm1(2),im_plane*RGB_norm1(3));
    RGB2 = inv(M)*([L_norm2; M_norm2; S_norm2]);
    RGB_norm2 = (0.5*RGB2)/norm(RGB2);
    im2 = cat(3,im_plane*RGB_norm2(1),im_plane*RGB_norm2(2),im_plane*RGB_norm2(3));
    im = cat(2,im1,im2);
     
    % Plotting the illumination spectra 
    figure(1),subplot(4,4,4*(ii-1)+1),plot(wave,illuminant(:,ii),'r','Linewidth',2); title('Illumination');
    xlabel('wavelength'); ylabel('Energy');
    
    % Plotting the reflectance spectra
    figure(1),subplot(4,4,4*(ii-1)+2),plot(wave,reflectance_spectra1,'g','Linewidth',2); hold on;
    plot(wave,reflectance_spectra2,'m','Linewidth',2); title('Reflectance');
    xlabel('wavelength'); ylabel('%'); hold off;
    
    % Plotting the net incident energy reaching the eye -> illuminant spectra * reflectance spectra
    figure(1),subplot(4,4,4*(ii-1)+3),plot(wave,net_spectra1,'k','Linewidth',2); hold on;
    plot(wave,net_spectra2,'g','Linewidth',2); title('Net spectra');
    xlabel('wavelength'), ylabel('Incident Energy'); hold off;
    
    % Plotting the color of the net spectra with respect to a gray background
    figure(1),subplot(4,4,4*ii),image(im); title('Image Color'); set(gca,'XTick',[],'YTick',[],'Box','on');

    L_cone_exc_ratios(ii) = L_norm1/L_norm2;
    M_cone_exc_ratios(ii) = M_norm1/M_norm2;
    S_cone_exc_ratios(ii) = S_norm1/S_norm2;
end

%% Trying to replicate Foster and Nacimento 1996 study about relational color constancy based on preserved cone excitation ratios

% Selecting just two illuminants
for jj = 1:4
    ntrials = 1000;
    size_illum = size(illuminant,2);
    illuminant1 = illuminant(:,randi(size_illum));
    illuminant2 = illuminant(:,randi(size_illum));
    size_munsell = size(munsell,2);
    exc_Lratios_ill1 = zeros(1,ntrials); exc_Lratios_ill2 = zeros(1,ntrials);
    exc_Mratios_ill1 = zeros(1,ntrials); exc_Mratios_ill2 = zeros(1,ntrials);
    exc_Sratios_ill1 = zeros(1,ntrials); exc_Sratios_ill2 = zeros(1,ntrials);
    card_RGratios_ill1 = zeros(1,ntrials); card_RGratios_ill2 = zeros(1,ntrials);
    card_BGratios_ill1 = zeros(1,ntrials); card_BGratios_ill2 = zeros(1,ntrials);
    card_lumratios_ill1 = zeros(1,ntrials); card_lumratios_ill2 = zeros(1,ntrials);
    for ii = 1:ntrials
        reflectance_spec1 = munsell(ind,randi(size_munsell));
        reflectance_spec2 = munsell(ind,randi(size_munsell));
        % Calculating the net light energy entering the eye of an edge
        surf1_ill1 = illuminant1'.*reflectance_spec1';
        surf2_ill1 = illuminant1'.*reflectance_spec2';
        surf1_ill2 = illuminant2'.*reflectance_spec1';
        surf2_ill2 = illuminant2'.*reflectance_spec2';
        
        % Calculating the L,M,S excitations/quantal catches from 2 different illuminants
        L_exc1_ill1 = surf1_ill1 * fundamentals(:,1); L_exc2_ill1 = surf2_ill1 * fundamentals(:,1);
        M_exc1_ill1 = surf1_ill1 * fundamentals(:,2); M_exc2_ill1 = surf2_ill1 * fundamentals(:,2);
        S_exc1_ill1 = surf1_ill1 * fundamentals(:,3); S_exc2_ill1 = surf2_ill1 * fundamentals(:,3);
        
        L_exc1_ill2 = surf1_ill2 * fundamentals(:,1); L_exc2_ill2 = surf2_ill2 * fundamentals(:,1);
        M_exc1_ill2 = surf1_ill2 * fundamentals(:,2); M_exc2_ill2 = surf2_ill2 * fundamentals(:,2);
        S_exc1_ill2 = surf1_ill2 * fundamentals(:,3); S_exc2_ill2 = surf2_ill2 * fundamentals(:,3);
        
        % Normalizing the responses
        [L_norm1_ill1,M_norm1_ill1,S_norm1_ill1,L_norm2_ill1,M_norm2_ill1,S_norm2_ill1] = normalize_resp(L_exc1_ill1,M_exc1_ill1,S_exc1_ill1,L_exc2_ill1,M_exc2_ill1,S_exc2_ill1,0);
        [L_norm1_ill2,M_norm1_ill2,S_norm1_ill2,L_norm2_ill2,M_norm2_ill2,S_norm2_ill2] = normalize_resp(L_exc1_ill2,M_exc1_ill2,S_exc1_ill2,L_exc2_ill2,M_exc2_ill2,S_exc2_ill2,0);
        
        % Calculating the cone excitation ratios 
        exc_Lratios_ill1(ii) = L_norm1_ill1/L_norm2_ill1;
        exc_Mratios_ill1(ii) = M_norm1_ill1/M_norm2_ill1;
        exc_Sratios_ill1(ii) = S_norm1_ill1/S_norm2_ill1;
        exc_Lratios_ill2(ii) = L_norm1_ill2/L_norm2_ill2;
        exc_Mratios_ill2(ii) = M_norm1_ill2/M_norm2_ill2;
        exc_Sratios_ill2(ii) = S_norm1_ill2/S_norm2_ill2;
        
        % Calculating the cardinal mechanism ratios 
        card_RGratios_ill1(ii) = (L_norm1_ill1-M_norm1_ill1)/(L_norm2_ill1-M_norm2_ill1);
        card_BYratios_ill1(ii) = (S_norm1_ill1-L_norm1_ill1-M_norm1_ill1)/(S_norm2_ill1-L_norm2_ill1-M_norm2_ill1);
        card_lumratios_ill1(ii) = (L_norm1_ill1+M_norm1_ill1)/(L_norm2_ill1+M_norm2_ill1);
        card_RGratios_ill2(ii) = (L_norm1_ill2-M_norm1_ill2)/(L_norm2_ill2-M_norm2_ill2);
        card_BYratios_ill2(ii) = (S_norm1_ill2-L_norm1_ill2-M_norm1_ill2)/(S_norm2_ill2-L_norm2_ill2-M_norm2_ill2);
        card_lumratios_ill2(ii) = (L_norm1_ill2+M_norm1_ill2)/(L_norm2_ill2+M_norm2_ill2);
    end
    
    figure(1+jj),subplot(231), plot(exc_Lratios_ill1,exc_Lratios_ill2 ,'o'); xlabel('ill1'), ylabel('ill2'); title('L');
    subplot(232), plot(exc_Mratios_ill1,exc_Mratios_ill2 ,'o'); xlabel('ill1'), ylabel('ill2'); title('M');
    subplot(233), plot(exc_Sratios_ill1,exc_Sratios_ill2 ,'o'); xlabel('ill1'), ylabel('ill2'); title('S');
    subplot(234), plot(card_RGratios_ill1, card_RGratios_ill2,'o'); xlabel('ill1'), ylabel('ill2'); title('L-M');
    subplot(235), plot(card_lumratios_ill1, card_lumratios_ill2,'o'); xlabel('ill1'), ylabel('ill2'); title('lum');
    subplot(236), plot(card_BYratios_ill1, card_BYratios_ill2,'o'); xlabel('ill1'), ylabel('ill2'); title('S/L-M');
end
