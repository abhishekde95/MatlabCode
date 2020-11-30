% Some of the sanity checks suggested by Greg,% The goal of this exercise is to keep two surfaces the same but have a fraction of incident light energy falling onto one of the
% surfaces. This light then falls onto the RF of a double opponent cell. We
% want to understand how this light excites the two subunits of a double
% opponent RF.
% Author - Abhishek 4/17
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

% Loading the cone action spectra and monitor spectral distributions
load fundamentals.mat
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
fundamentals = fundamentals(3:end,:); % Starting the fundamentals from 390 nm
mon_spd = mon_spd(:,3:end);
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations

% Loading the reflectance spectra of two surfaces and place them adjacent to each other
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
L = numel(380:1:780); % for the munsell reflectance spectras - ftp://ftp.cs.joensuu.fi/pub/color/spectra/mspec/README.txt
ind = 1:5:L;
ind = ind(3:end);
wave = 390:5:780;
frac  = 0.5;
hi = size(illuminant,1);
RGBS1 = [0.3; -0.3; 0.0]; RGBS2 = [0.3; -0.3; 0];
color = [0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0];
figure(1),subplot(121);
avg_Subunit_act1 = [];
avg_Subunit_act2 = [];
for jj = 1:5
    % Declaring the image plane
    Lnorm1 = zeros(1,hi);
    Mnorm1 = zeros(1,hi);
    Snorm1 = zeros(1,hi);
    Lnorm2 = zeros(1,hi);
    Mnorm2 = zeros(1,hi);
    Snorm2 = zeros(1,hi);
    RGBsurf1 = [];
    RGBsurf2 = [];
    rand_idxs = randi(size(munsell,2),1); % choose 1 random number from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs);
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = reflectance_spectra1;
    
    for ii = 1:1:hi
        % Calculating the net light energy entering the eye of an edge
        net_spectra1 = illuminant(ii,:).*reflectance_spectra1;
        net_spectra2 = frac*(illuminant(ii,:).*reflectance_spectra2);
               
        % Calculating the L,M,S excitations/quantal catches from the incident light
        L_exc1 = net_spectra1 * fundamentals(:,1);
        M_exc1 = net_spectra1 * fundamentals(:,2);
        S_exc1 = net_spectra1 * fundamentals(:,3);
        L_exc2 = net_spectra2 * fundamentals(:,1);
        M_exc2 = net_spectra2 * fundamentals(:,2);
        S_exc2 = net_spectra2 * fundamentals(:,3);
        
        % Calculating the RGB intensities from the L,M,S excitations after
        % normalizng the responses, 0 - don't normalize, 1 - do normalize
        [L1,M1,S1,L2,M2,S2] = normalize_resp(L_exc1,M_exc1,S_exc1,L_exc2,M_exc2,S_exc2,0);
        RGB1 = inv(M)*([L1; M1; S1]);
        RGB_norm1 = (RGB1);
        RGB2 = inv(M)*([L2; M2; S2]);
        RGB_norm2 = (RGB2);
        Lnorm1(ii) = L1;
        Mnorm1(ii) = M1;
        Snorm1(ii) = S1;
        Lnorm2(ii) = L2;
        Mnorm2(ii) = M2;
        Snorm2(ii) = S2;
        RGBsurf1 = [RGBsurf1; RGB_norm1'];
        RGBsurf2 = [RGBsurf2; RGB_norm2'];
    end
    Subunit_act1 = RGBsurf1 * RGBS1;
    Subunit_act2 = RGBsurf2 * RGBS2;
    avg_Subunit_act1 = [avg_Subunit_act1; mean(Subunit_act1)];
    avg_Subunit_act2 = [avg_Subunit_act2; mean(Subunit_act2)];
    plot([Subunit_act1;0], [Subunit_act2; 0],'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',color(jj,:));
    xlabel('S1'), ylabel('S2'); grid on; hold on;
end
figure(1),subplot(121),title('Changing reflectances'); hold off;

% Now keeping the reflectances the same but changing the fraction of incident light energy falling onto one of the surfaces
rand_idxs = randi(size(munsell,2),1); % choose 2 random numbers from 1269 possible reflectance surfaces
reflectance_spectra1 = munsell(ind,rand_idxs);
reflectance_spectra1 = reflectance_spectra1';
reflectance_spectra2 = reflectance_spectra1;
RGBS1 = [0.3; -0.5; 0]; RGBS2 = [-0.3; 0.5; 0];
figure(1); subplot(122);
for jj = 1:6
    % Declaring the image plane
    Lnorm1 = zeros(1,hi);
    Mnorm1 = zeros(1,hi);
    Snorm1 = zeros(1,hi);
    Lnorm2 = zeros(1,hi);
    Mnorm2 = zeros(1,hi);
    Snorm2 = zeros(1,hi);
    RGBsurf1 = [];
    RGBsurf2 = [];
    frac  = 0.1*jj;
        
    for ii = 1:1:hi
        % Calculating the net light energy entering the eye of an edge
        net_spectra1 = illuminant(ii,:).*reflectance_spectra1;
        net_spectra2 = frac*(illuminant(ii,:).*reflectance_spectra2);
        
        % Calculating the cone excitatios pretending as if the surface is a
        % perfect reflector
        illL1 = illuminant(ii,:) * fundamentals(:,1);
        illM1 = illuminant(ii,:) * fundamentals(:,2);
        illS1 = illuminant(ii,:) * fundamentals(:,3);
        
        % Calculating the L,M,S excitations/quantal catches from the incident light
        L_exc1 = net_spectra1 * fundamentals(:,1);
        M_exc1 = net_spectra1 * fundamentals(:,2);
        S_exc1 = net_spectra1 * fundamentals(:,3);
        L_exc2 = net_spectra2 * fundamentals(:,1);
        M_exc2 = net_spectra2 * fundamentals(:,2);
        S_exc2 = net_spectra2 * fundamentals(:,3);
        
        % Calculating the RGB intensities from the L,M,S excitations after
        % normalizng the responses, 0 - don't normalize, 1 - do normalize
        [L1,M1,S1,L2,M2,S2] = normalize_resp(L_exc1,M_exc1,S_exc1,L_exc2,M_exc2,S_exc2,0);
        RGB1 = inv(M)*([L1; M1; S1]);
        RGB_norm1 = RGB1; % Without normalizing it
        RGB2 = inv(M)*([L2; M2; S2]);
        RGB_norm2 = RGB2;
        Lnorm1(ii) = L1;
        Mnorm1(ii) = M1;
        Snorm1(ii) = S1;
        Lnorm2(ii) = L2;
        Mnorm2(ii) = M2;
        Snorm2(ii) = S2;
        RGBsurf1 = [RGBsurf1; RGB_norm1'];
        RGBsurf2 = [RGBsurf2; RGB_norm2'];
    end
    Subunit_act1 = RGBsurf1 * RGBS1;
    Subunit_act2 = RGBsurf2 * RGBS2;
    plot([Subunit_act1;0], [Subunit_act2; 0],'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',color(jj,:));
    xlabel('S1'), ylabel('S2'); grid on; hold on;
end
figure(1), subplot(122); title('Changing frac');hold off;

%% Calculating where does the illumination (daylight) spectra lie in the LMS space
illL = zeros(1,hi);
illM = zeros(1,hi);
illS = zeros(1,hi);
for ii = 1:1:hi
    % Calculating the cone excitations pretending as if the surface is a perfect reflector
    illL(ii) = illuminant(ii,:) * fundamentals(:,1);
    illM(ii) = illuminant(ii,:) * fundamentals(:,2);
    illS(ii) = illuminant(ii,:) * fundamentals(:,3);
end
figure(2), subplot(221), plot3(illL,illM,illS,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); grid on;
xlabel('L'), ylabel('M'), zlabel('S'); title('Daylight illumination spectra');

% Calculating where does the reflectance spectras of munsell chips lie in the LMS space 
L = size(munsell,2); 
reflL = zeros(1,L);
reflM = zeros(1,L);
reflS = zeros(1,L);
for ii = 1:L
    % Calculating the cone excitatios of the reflectance spectras
    reflL(ii) = munsell(ind,ii)' * fundamentals(:,1);
    reflM(ii) = munsell(ind,ii)' * fundamentals(:,2);
    reflS(ii) = munsell(ind,ii)' * fundamentals(:,3);
end
figure(2), subplot(222), plot3(reflL,reflM,reflS,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); grid on;
xlabel('L'), ylabel('M'), zlabel('S'); title('Reflectance spectra');

% Calculating where the does reflectance spectras of tree leaves lie in the LMS space
load tree.asc
T = size(tree,1); 
tree_spec = [];
for ii = 1:T
    tmp = tree(ii,:);
    splinefit = spline(wave_il, tmp, [390:5:1070]);
    tree_spec = [tree_spec; splinefit(1:hi_wave_ind)];
end
treeL = zeros(1,T);
treeM = zeros(1,T);
treeS = zeros(1,T);
for ii = 1:T
    % Calculating the cone excitatios of the reflectance spectras
    treeL(ii) = tree_spec(ii,:) * fundamentals(:,1);
    treeM(ii) = tree_spec(ii,:) * fundamentals(:,2);
    treeS(ii) = tree_spec(ii,:) * fundamentals(:,3);
end
figure(2), subplot(223), plot3(treeL,treeM,treeS,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); grid on;
xlabel('L'), ylabel('M'), zlabel('S'); title('Tree spectra');
%
load birch.mat
B = size(birch,2);
birch_spec = []; 
for ii = 1:B
    tmp = birch(1:hi_wave_ind,ii)';
    birch_spec = [birch_spec; tmp];
end
birchL = zeros(1,B);
birchM = zeros(1,B);
birchS = zeros(1,B);
for ii = 1:T
    % Calculating the cone excitations of the reflectance spectras
    birchL(ii) = birch_spec(ii,:) * fundamentals(:,1);
    birchM(ii) = birch_spec(ii,:) * fundamentals(:,2);
    birchS(ii) = birch_spec(ii,:) * fundamentals(:,3);
end
figure(2), subplot(224), plot3(birchL,birchM,birchS,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 1 0]); grid on;
xlabel('L'), ylabel('M'), zlabel('S'); title('Birch spectra');

% Now I need to plot all the spectras
figure(3), subplot(221),plot(wave,illuminant'), title('Illumination Spectra');
subplot(222),plot(wave,munsell(ind,:)'), title('Reflectance Spectra');
subplot(223),plot(wave,tree_spec'), title('Tree Spectra');
subplot(224),plot(wave,birch_spec'), title('Birch Spectra');

% Now I need to plot the data in DKL space
figure(4), subplot(221), plot3(illL+illM+illS, illL-illM, illS-(illL+illM),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); grid on;
xlabel('L+M+S'), ylabel('L-M'), zlabel('S-(L+M)'); title('Illumination spectra');
subplot(222), plot3(reflL+reflM+reflS, reflL-reflM, reflS-(reflL+reflM),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); grid on;
xlabel('L+M+S'), ylabel('L-M'), zlabel('S-(L+M)'); title('Reflectance spectra');
subplot(223), plot3(treeL+treeM+treeS, treeL-treeM, treeS-(treeL+treeM),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); grid on;
xlabel('L+M+S'), ylabel('L-M'), zlabel('S-(L+M)'); title('Tree spectra');
subplot(224), plot3(birchL+birchM+birchS, birchL-birchM, birchS-(birchL+birchM),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 1 0]); grid on;
xlabel('L+M+S'), ylabel('L-M'), zlabel('S-(L+M)'); title('Birch spectra');


%% Now I want to see how the cone excitations change when the illumination is scaled keeping its relative power spectral density constant
fracs = linspace(0,3,1001);
RGBS1 = [0.5; -0.5; 0];
RGBS2 = [-0.5; 0.5; 0];
rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
reflectance_spectra1 = munsell(ind,rand_idxs(1));
reflectance_spectra1 = reflectance_spectra1';
reflectance_spectra2 = munsell(ind,rand_idxs(2));
reflectance_spectra2 = reflectance_spectra2';
Subunit1_act = [];
Subunit2_act = [];
for ii = 1:numel(fracs)
    net_spectra1 = fracs(ii)*(illuminant(1,:).*reflectance_spectra1);
    net_spectra2 = fracs(ii)*(illuminant(1,:).*reflectance_spectra2);
    
    % Calculating the L,M,S excitations/quantal catches from the incident light
    L1 = net_spectra1 * fundamentals(:,1);
    M1 = net_spectra1 * fundamentals(:,2);
    S1 = net_spectra1 * fundamentals(:,3);
    L2 = net_spectra2 * fundamentals(:,1);
    M2 = net_spectra2 * fundamentals(:,2);
    S2 = net_spectra2 * fundamentals(:,3);
    
    RGB1 = inv(M)*([L1; M1; S1]);
    RGB2 = inv(M)*([L2; M2; S2]);
    Subunit1_act = [Subunit1_act; RGB1'*RGBS1];
    Subunit2_act = [Subunit2_act; RGB2'*RGBS2];
end
figure(5), subplot(121),plot(Subunit1_act, Subunit2_act,'o','MarkerSize',5,'Linewidth',0.5,'MarkerFaceColor',[0 0 1]); grid on;
xlabel('S1'), ylabel('S2'); 

% Next I want to see if changes in illumination spectrum and changes in illuminant amplitude yield different trajectories in LMS space
Subunit1_act1 = [];
Subunit2_act1 = [];
for ii = 1:hi
    net_spectra1 = (illuminant(ii,:).*reflectance_spectra1);
    net_spectra2 = (illuminant(ii,:).*reflectance_spectra2);
    
    % Calculating the L,M,S excitations/quantal catches from the incident light
    L1 = net_spectra1 * fundamentals(:,1);
    M1 = net_spectra1 * fundamentals(:,2);
    S1 = net_spectra1 * fundamentals(:,3);
    L2 = net_spectra2 * fundamentals(:,1);
    M2 = net_spectra2 * fundamentals(:,2);
    S2 = net_spectra2 * fundamentals(:,3);
    
    RGB1 = inv(M)*([L1; M1; S1]);
    RGB2 = inv(M)*([L2; M2; S2]);
    Subunit1_act1 = [Subunit1_act1; RGB1'*RGBS1];
    Subunit2_act1 = [Subunit2_act1; RGB2'*RGBS2];
end

figure(5), subplot(121), hold on;
plot(Subunit1_act1, Subunit2_act1,'o','MarkerSize',5,'Linewidth',0.5,'MarkerFaceColor',[1 0 0]); hold off;
figure(5), subplot(122), plot(illuminant'); xlabel('wavelength'),ylabel('Power'), title('SPD illuminant');

%% Loading illuminant SPD from Psychtoolbox
% This cell did not prove to useful as the units of the SPD from
% psychtoolbox and the daylight database are different
load spd_CIEA.mat
load spd_CIEC.mat
load spd_D65.mat
load spd_flourescent.mat
load spd_incanCC.mat
load spd_phillybright.mat
load spd_xenonArc.mat
load spd_xenonFlash.mat
spd_phillybright = spline([380:4:780], spd_phillybright, [380:5:780]); % fitting a cubic spline
% The spectrum of these illuminants range from 380-780 nm which is
% different from the range of the daylight spectrum which is 390-780nm,
% therefore in order to maintain consistency, I am only including values
% between 390-790 nm
art_illuminant = [spd_CIEA spd_CIEC spd_D65 spd_flourescent spd_incanCC spd_phillybright' spd_xenonFlash];
art_illuminant = art_illuminant(3:end,:)';
figure(6),subplot(121),plot(wave,illuminant'),
subplot(122),plot(wave,art_illuminant');

% So now basically try the same exercise as before 
all_illuminants = cat(1,1000*illuminant,art_illuminant); % Using a multiplicative factor of 1000 for 
L = size(all_illuminants,1);
Subunit1_act2 = [];
Subunit2_act2 = [];
for ii = 1:L
    net_spectra1 = (all_illuminants(ii,:).*reflectance_spectra1);
    net_spectra2 = (all_illuminants(ii,:).*reflectance_spectra2);
    
    % Calculating the L,M,S excitations/quantal catches from the incident light
    L1 = net_spectra1 * fundamentals(:,1);
    M1 = net_spectra1 * fundamentals(:,2);
    S1 = net_spectra1 * fundamentals(:,3);
    L2 = net_spectra2 * fundamentals(:,1);
    M2 = net_spectra2 * fundamentals(:,2);
    S2 = net_spectra2 * fundamentals(:,3);
    
    RGB1 = inv(M)*([L1; M1; S1]);
    RGB2 = inv(M)*([L2; M2; S2]);
    Subunit1_act2 = [Subunit1_act2; RGB1'*RGBS1];
    Subunit2_act2 = [Subunit2_act2; RGB2'*RGBS2];
    
end

figure(7), plot(1000*Subunit1_act, 1000*Subunit2_act,'o','MarkerSize',5,'Linewidth',0.5,'MarkerFaceColor',[0 0 1]); hold on;
plot(Subunit1_act2, Subunit2_act2,'o','MarkerSize',5,'Linewidth',0.5,'MarkerFaceColor',[1 0 0]); grid on;
xlabel('S1'), ylabel('S2'); 

% Projecting all the illumination spectras into LMS space, both the previous spectras and the illuminants from Psychtoolbox
all_illL = zeros(1,L);
all_illM = zeros(1,L);
all_illS = zeros(1,L);
for ii = 1:1:L
    % Calculating the cone excitations pretending as if the surface is a perfect reflector
    all_illL(ii) = all_illuminants(ii,:) * fundamentals(:,1);
    all_illM(ii) = all_illuminants(ii,:) * fundamentals(:,2);
    all_illS(ii) = all_illuminants(ii,:) * fundamentals(:,3);
end
figure(8), subplot(121),plot3(all_illL,all_illM,all_illS,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); grid on;
xlabel('L'), ylabel('M'), zlabel('S'); title('Daylight illumination spectra');
subplot(122),plot(wave,all_illuminants'), xlabel('Wavelength'), ylabel('Energy');

%% Trying a new code where I want to see how a single STA (neuron) responds to the different edges and see when it could use 
% linear and nonlinear operation 

% Choosing the STA of the cell 
RGBS1_DO = [0.3; -0.3;0]; % Double opponent cell
RGBS2_DO = [-0.3; 0.3; 0];
RGBS1 = [0.3; -0.1; 0.2]; % Some random cell, not a double opponent cell
RGBS2 = [0.15; 0.5; 0.1];
RGBS1_lum = [0.3; 0.2; 0.1]; % A luminance cell 
RGBS2_lum = [-0.3; -0.2; -0.1];
RGBS1_unm = [0.3; -0.5; 0.0]; % DO but unmatched  
RGBS2_unm = [-0.3; 0.2; 0.0];
hi = size(illuminant,1);
N = 25;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs(1));
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = munsell(ind,rand_idxs(2));
    reflectance_spectra2 = reflectance_spectra2';
    Subunit1_act_DO = [];
    Subunit2_act_DO = [];
    Subunit1_act = [];
    Subunit2_act = [];
    Subunit1_act_lum = [];
    Subunit2_act_lum = [];
    Subunit1_act_unm = [];
    Subunit2_act_unm = [];
    for ii = 1:hi
        net_spectra1 = illuminant(ii,:).*reflectance_spectra1;
        net_spectra2 = (illuminant(ii,:).*reflectance_spectra2);
        
        L1 = net_spectra1 * fundamentals(:,1);
        M1 = net_spectra1 * fundamentals(:,2);
        S1 = net_spectra1 * fundamentals(:,3);
        L2 = net_spectra2 * fundamentals(:,1);
        M2 = net_spectra2 * fundamentals(:,2);
        S2 = net_spectra2 * fundamentals(:,3);
        
        RGB1 = (inv(M)*log10(([L1; M1; S1])));
        RGB2 = (inv(M)*log10(([L2; M2; S2])));
        Subunit1_act_DO = [Subunit1_act_DO; RGB1'*RGBS1_DO];
        Subunit2_act_DO = [Subunit2_act_DO; RGB2'*RGBS2_DO];
        Subunit1_act = [Subunit1_act; RGB1'*RGBS1];
        Subunit2_act = [Subunit2_act; RGB2'*RGBS2];
        Subunit1_act_lum = [Subunit1_act_lum; RGB1'*RGBS1_lum];
        Subunit2_act_lum = [Subunit2_act_lum; RGB2'*RGBS2_lum];
        Subunit1_act_unm = [Subunit1_act_unm; RGB1'*RGBS1_unm];
        Subunit2_act_unm = [Subunit2_act_unm; RGB2'*RGBS2_unm];
    end
    figure(9), subplot(n,n,jj), plot(Subunit1_act_DO, Subunit2_act_DO,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal, drawnow;
%     figure(10), subplot(n,n,jj), plot(Subunit1_act, Subunit2_act,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); drawnow;
%     figure(11), subplot(n,n,jj), plot(Subunit1_act_lum, Subunit2_act_lum,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); drawnow;
%     figure(12), subplot(n,n,jj), plot(Subunit1_act_unm, Subunit2_act_unm,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 1]); drawnow;
end

%% Tying out few more codes based on meeting with Greg 4/14
% Trying out things in LMS space, both lights and STA in LMS

% Choosing the STA of the cell 
LMSS1_DO = [0.5; -0.5;0]; % Double opponent cell
LMSS2_DO = [-0.5; 0.5; 0];

hi = size(illuminant,1);
N = 9;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs(1));
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = munsell(ind,rand_idxs(2));
    reflectance_spectra2 = reflectance_spectra2';
    Subunit1_act_DO = [];
    Subunit2_act_DO = [];
    Subunit1log = [];
    Subunit2log = [];
    Subunit1logrect = [];
    Subunit2logrect = [];
    for ii = 1:hi
        net_spectra1 = illuminant(ii,:).*reflectance_spectra1;
        net_spectra2 = (illuminant(ii,:).*reflectance_spectra2);
        
        L1 = net_spectra1 * fundamentals(:,1);
        M1 = net_spectra1 * fundamentals(:,2);
        S1 = net_spectra1 * fundamentals(:,3);
        L2 = net_spectra2 * fundamentals(:,1);
        M2 = net_spectra2 * fundamentals(:,2);
        S2 = net_spectra2 * fundamentals(:,3);
        
        LMS1 = [L1; M1; S1];
        LMS2 = [L2; M2; S2];
        Subunit1_act_DO = [Subunit1_act_DO; LMS1'*LMSS1_DO];
        Subunit2_act_DO = [Subunit2_act_DO; LMS2'*LMSS2_DO];
        resp1 = log10(LMS1')*LMSS1_DO;
        resp2 = log10(LMS2')*LMSS2_DO;
        Subunit1logrect = [Subunit1logrect; resp1*(resp1>0)];
        Subunit2logrect = [Subunit2logrect; resp2*(resp2>0)];
        
        Subunit1log = [Subunit1log; resp1];
        Subunit2log = [Subunit2log; resp2];
        
    end
    figure(13), subplot(n,n,jj), plot(Subunit1_act_DO, Subunit2_act_DO,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal, drawnow;
    figure(14), subplot(n,n,jj), plot(Subunit1log, Subunit2log,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); axis equal, drawnow;
    figure(15), subplot(n,n,jj), plot(Subunit1logrect, Subunit2logrect,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); axis equal, drawnow;


end
