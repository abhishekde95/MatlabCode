% Writing a new script daylight_props_4.m
% Post second thesis committee meeting
% Trying out some control cases

%% Case I: using ISETBIO illuminant spectras

close all; clearvars;
wave = 390:5:780;
dayBasis = ieReadSpectra('cieDaylightBasis',wave); % Daylight spectra basis functions from isetbio
num_spectras = 100;
coeffs = cat(2,ones(num_spectras,1),2*rand(num_spectras,2)-1); % Limiting the coefficients between 0 and 1
illuminants = coeffs * dayBasis';
figure(1), subplot(121), plot(wave,dayBasis'),xlabel('Wavelength'), ylabel('Energy'); set(gca,'Xlim',[390 780]);
subplot(122),plot(wave,illuminants'), xlabel('Wavelength'), ylabel('Energy'); set(gca,'Xlim',[390 780]);

% Loading the cone action spectra and monitor spectral distributions
load fundamentals.mat
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
fundamentals = fundamentals(3:end,:); % Starting the fundamentals from 390 nm
mon_spd = mon_spd(:,3:end);
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations

% One can also use the munsell chips spectra from psychtoolbox stored as
% sur_nickerson.mat. An important thing to note would be that the SPDs are
% measured from 380 to 780 nm in steps of 5nm gap
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
L = numel(380:1:780); % for the munsell reflectance spectras - ftp://ftp.cs.joensuu.fi/pub/color/spectra/mspec/README.txt
ind = 1:5:L;
ind = ind(3:end);

%% RF of a double opponent cell
% The subunit vector represents cone weights
RGDO_S1 = [0.3; -0.6; 0.0]; RGDO_S2 = [-0.3; 0.6; 0];
BYDO_S1 = [0.15;0.15;-0.6]; BYDO_S2 = [-0.15;-0.15;0.6];
hi = size(illuminants,1);
N = 16;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs(1));
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = munsell(ind,rand_idxs(2));
    reflectance_spectra2 = reflectance_spectra2';
    Subunit1_act_RGDO = [];
    Subunit2_act_RGDO = [];
    Subunit1_act_BYDO = [];
    Subunit2_act_BYDO = [];
    for ii = 1:hi
        net_spectra1 = illuminants(ii,:).*reflectance_spectra1;
        net_spectra2 = (illuminants(ii,:).*reflectance_spectra2);
        
        L1 = net_spectra1 * fundamentals(:,1);
        M1 = net_spectra1 * fundamentals(:,2);
        S1 = net_spectra1 * fundamentals(:,3);
        L2 = net_spectra2 * fundamentals(:,1);
        M2 = net_spectra2 * fundamentals(:,2);
        S2 = net_spectra2 * fundamentals(:,3);
        
        LMS1 = [L1; M1; S1];
        LMS2 = [L2; M2; S2];
        Subunit1_act_RGDO = [Subunit1_act_RGDO; LMS1'*RGDO_S1];
        Subunit2_act_RGDO = [Subunit2_act_RGDO; LMS2'*RGDO_S2];
        Subunit1_act_BYDO = [Subunit1_act_BYDO; LMS1'*BYDO_S1];
        Subunit2_act_BYDO = [Subunit2_act_BYDO; LMS2'*BYDO_S2];
    end
    figure(2), subplot(n,n,jj), plot(Subunit1_act_RGDO, Subunit2_act_RGDO,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); axis equal; drawnow;
    figure(3), subplot(n,n,jj), plot(Subunit1_act_BYDO, Subunit2_act_BYDO,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal; drawnow;
    
end
%% Calculating the mean cone excitations so that in future I can deal with cone contrasts.
% Taking the first basis function as the light that reaches the cones after
% being reflected from a mirror(uniform broad spectral power disribution)
mean_cone_exc = fundamentals'*dayBasis(:,1); % Mean L,M,S cone excitations
% RF of a double opponent cell
% The subunit vector represents cone weights
RGDO_S1 = [0.3; -0.3; 0.0]; RGDO_S2 = [-0.3; 0.3; 0];
BYDO_S1 = [0.15;0.15;-0.3]; BYDO_S2 = [-0.15;-0.15;0.3];
hi = size(illuminants,1);
N = 9;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs(1));
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = munsell(ind,rand_idxs(2));
    reflectance_spectra2 = reflectance_spectra2';
    Subunit1_act_RGDO = [];
    Subunit2_act_RGDO = [];
    Subunit1_act_BYDO = [];
    Subunit2_act_BYDO = [];
    for ii = 1:hi
        net_spectra1 = illuminants(ii,:).*reflectance_spectra1;
        net_spectra2 = (illuminants(ii,:).*reflectance_spectra2);
        
        L1 = net_spectra1 * fundamentals(:,1);
        M1 = net_spectra1 * fundamentals(:,2);
        S1 = net_spectra1 * fundamentals(:,3);
        L2 = net_spectra2 * fundamentals(:,1);
        M2 = net_spectra2 * fundamentals(:,2);
        S2 = net_spectra2 * fundamentals(:,3);
        
        LMS1 = [L1; M1; S1];
        LMS2 = [L2; M2; S2];
        
        % Converting cone excitations into cone contrasts (Weber's
        % contrast)
        
        LMS1 = (LMS1 - mean_cone_exc)./mean_cone_exc; % cone contrasts
        LMS2 = (LMS2 - mean_cone_exc)./mean_cone_exc; % cone contrasts
        
        Subunit1_act_RGDO = [Subunit1_act_RGDO; LMS1'*RGDO_S1];
        Subunit2_act_RGDO = [Subunit2_act_RGDO; LMS2'*RGDO_S2];
        Subunit1_act_BYDO = [Subunit1_act_BYDO; LMS1'*BYDO_S1];
        Subunit2_act_BYDO = [Subunit2_act_BYDO; LMS2'*BYDO_S2];
    end
    figure(4), subplot(n,n,jj), plot(Subunit1_act_RGDO, Subunit2_act_RGDO,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal, drawnow;
    figure(5), subplot(n,n,jj), plot(Subunit1_act_BYDO, Subunit2_act_BYDO,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal, drawnow;
    
end

%% Seems like an interesting direction to head towards. 
L = 10; W = 10; % Creating image and RF of size 10x10 pixels
RGDO_S1 = [0.3; -0.3; 0.0]; RGDO_S2 = [-0.3; 0.3; 0];
S1 = cat(2,ones(10,5),zeros(10,5));
S1 = cat(3,RGDO_S1(1)*S1,RGDO_S1(2)*S1,RGDO_S1(3)*S1);
S2 = cat(2,zeros(10,5),ones(10,5));
S2 = cat(3,RGDO_S2(1)*S2,RGDO_S2(2)*S2,RGDO_S2(3)*S2);
hi = size(illuminants,1);
len = size(munsell,2);
mean_cone_exc = fundamentals'*dayBasis(:,1);
mean_cone_exc = cat(3,mean_cone_exc(1)*ones(L,W),mean_cone_exc(2)*ones(L,W),mean_cone_exc(3)*ones(L,W));
N = 9; n = sqrt(N);
for mm = 1:N
    LMS = zeros(L,W,3);
    surf_ind = randi(len,[L W]);
    S1_act = []; S2_act = [];
    for ii = 1:hi
        for jj = 1:L
            for kk = 1:W
                reflectance_spectra = munsell(ind,surf_ind(jj,kk));
                net_spectra = illuminants(ii,:).*reflectance_spectra';
                L1 = net_spectra * fundamentals(:,1);
                M1 = net_spectra * fundamentals(:,2);
                S1 = net_spectra * fundamentals(:,3);
                LMS(jj,kk,1) = L1;
                LMS(jj,kk,2) = M1;
                LMS(jj,kk,3) = S1;
            end
        end
        LMS = (LMS - mean_cone_exc)./mean_cone_exc;
        tmp = S1.*LMS; S1_act = [S1_act; sum(tmp(:))];
        tmp = S2.*LMS; S2_act = [S2_act; sum(tmp(:))];
    end
    figure(6),subplot(n,n,mm); plot(S1_act,S2_act,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis square;
end

%% Code that is similar to the previous code but I am changing sizes of the blocks of the squares to see where
% does the nonlinearity convert to linearity
% ( Working on it currently, not finished it yet)
% I am not sure if this is worth is pursuing, investing time in a more
% interesting simulation


%% Case II: Trying to create hypothetical gaussian shaped illumination SPD
wave = 390:5:780;
gauss_means = 400:10:800;
gauss_spectras_bwl = zeros(numel(gauss_means),numel(wave));
gauss_spectras_bwh = zeros(numel(gauss_means),numel(wave));
gauss_spectras_scaled = zeros(numel(gauss_means),numel(wave));
for ii = 1:numel(gauss_means)
    gauss_spectras_bwl(ii,:) = 10 * normpdf(wave,gauss_means(ii),40); % low bandwidth
    gauss_spectras_bwh(ii,:) = 10 * normpdf(wave,gauss_means(ii),100); % high bandwidth
    gauss_spectras_scaled(ii,:) = ii* normpdf(wave,gauss_means(20),40); % the illumination spectras are scaled version of each other
end
figure(7),subplot(221),plot(wave, gauss_spectras_bwl','Linewidth',1); set(gca,'Xlim',[wave(1) wave(end)]);
subplot(222),plot(wave, gauss_spectras_bwh','Linewidth',1); set(gca,'Xlim',[wave(1) wave(end)]);
subplot(223),plot(wave, gauss_spectras_scaled','Linewidth',1); set(gca,'Xlim',[wave(1) wave(end)]);

% A similar simulation as done before but with new artificial illumination spectras
% The subunit vector represents cone weights
RGDO_S1 = [0.3; -0.3; 0.2]; RGDO_S2 = [-0.3; 0.3; -0.2];
hi = numel(gauss_means);
N = 16;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs(1));
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = munsell(ind,rand_idxs(2));
    reflectance_spectra2 = reflectance_spectra2';
    Subunit1_act_RGDO_bwl = []; Subunit2_act_RGDO_bwl = [];
    Subunit1_act_RGDO_bwh = []; Subunit2_act_RGDO_bwh = [];
    Subunit1_act_RGDO_scaled = []; Subunit2_act_RGDO_scaled = [];
    for ii = 1:hi
        % low bandwidth
        net_spectra1 = gauss_spectras_bwl(ii,:).*reflectance_spectra1;
        net_spectra2 = gauss_spectras_bwl(ii,:).*reflectance_spectra2;
        [LMS1, LMS2] = getLMSvals(net_spectra1, net_spectra2, fundamentals);
        Subunit1_act_RGDO_bwl = [Subunit1_act_RGDO_bwl; LMS1'*RGDO_S1];
        Subunit2_act_RGDO_bwl = [Subunit2_act_RGDO_bwl; LMS2'*RGDO_S2];
        
        % high bandwidth
        net_spectra1 = gauss_spectras_bwh(ii,:).*reflectance_spectra1;
        net_spectra2 = gauss_spectras_bwh(ii,:).*reflectance_spectra2;
        [LMS1, LMS2] = getLMSvals(net_spectra1, net_spectra2, fundamentals);
        Subunit1_act_RGDO_bwh = [Subunit1_act_RGDO_bwh; LMS1'*RGDO_S1];
        Subunit2_act_RGDO_bwh = [Subunit2_act_RGDO_bwh; LMS2'*RGDO_S2];
        
        % scaled illumination spectras
        net_spectra1 = gauss_spectras_scaled(ii,:).*reflectance_spectra1;
        net_spectra2 = gauss_spectras_scaled(ii,:).*reflectance_spectra2;
        [LMS1, LMS2] = getLMSvals(net_spectra1, net_spectra2, fundamentals);
        Subunit1_act_RGDO_scaled = [Subunit1_act_RGDO_scaled; LMS1'*RGDO_S1];
        Subunit2_act_RGDO_scaled = [Subunit2_act_RGDO_scaled; LMS2'*RGDO_S2];
    end
    figure(8), subplot(n,n,jj), plot(Subunit1_act_RGDO_bwl, Subunit2_act_RGDO_bwl,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); axis equal, drawnow; set(gca,'XTick',[],'YTick',[]);
    figure(9), subplot(n,n,jj), plot(Subunit1_act_RGDO_bwh, Subunit2_act_RGDO_bwh,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); axis equal, drawnow; set(gca,'XTick',[],'YTick',[]);
    figure(10), subplot(n,n,jj), plot(Subunit1_act_RGDO_scaled, Subunit2_act_RGDO_scaled,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal, drawnow; set(gca,'XTick',[],'YTick',[]);
    figure(11),subplot(n,n,jj), plot([reflectance_spectra1; reflectance_spectra2]', 'Linewidth',2); drawnow; set(gca,'XTick',[],'YTick',[]);
end

%%  Case III: Trying to see if the bandwidth of the reflectance spectra has something to do with the spatial computation required
wave = 390:5:780;
gauss_means = 400:10:800;
art_reflectance_spectras_bwl = zeros(numel(gauss_means),numel(wave));
art_reflectance_spectras_bwh = zeros(numel(gauss_means),numel(wave));
for ii = 1:numel(gauss_means)
    art_reflectance_spectras_bwl(ii,:) = 10 * normpdf(wave,gauss_means(ii),10); % low bandwidth
    art_reflectance_spectras_bwh(ii,:) = 80 * normpdf(wave,gauss_means(ii),100); % low bandwidth
end

% A similar simulation as done before but with new artificial illumination spectras
% The subunit vector represents cone weights
RGDO_S1 = [0.3; -0.3; 0]; RGDO_S2 = [-0.3; 0.3; 0];
hi = numel(gauss_means);
N = 16;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(numel(gauss_means),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = art_reflectance_spectras_bwl(rand_idxs(1),:);
    reflectance_spectra2 = art_reflectance_spectras_bwl(rand_idxs(2),:);
    Subunit1_act_RGDO = []; Subunit2_act_RGDO = [];
    for ii = 1:hi
        % low bandwidth
        net_spectra1 = gauss_spectras_bwh(ii,:).*reflectance_spectra1;
        net_spectra2 = gauss_spectras_bwh(ii,:).*reflectance_spectra2;
        [LMS1, LMS2] = getLMSvals(net_spectra1, net_spectra2, fundamentals);
        Subunit1_act_RGDO = [Subunit1_act_RGDO; LMS1'*RGDO_S1];
        Subunit2_act_RGDO = [Subunit2_act_RGDO; LMS2'*RGDO_S2];
    end
    figure(12), subplot(n,n,jj), plot(Subunit1_act_RGDO, Subunit2_act_RGDO,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0]); axis equal, drawnow; set(gca,'XTick',[],'YTick',[]);
    figure(13), subplot(n,n,jj), plot([reflectance_spectra1; reflectance_spectra2]', 'Linewidth',2); drawnow; set(gca,'XTick',[],'YTick',[]);
end

%% Case IV: Trying out a new follow up exercise to see how the chromaticity of a surfaces changes over different times of the day
bandwidths = 10:40:200;
mean_cone_exc = fundamentals'*dayBasis(:,1);
color = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0];
for jj = 1: numel(bandwidths)
    surf_rs= 10 * normpdf(wave,600,bandwidths(jj)); % low bandwidth
    chromaticity1 = []; % for storing chromaticity obtained from cone excitations
    chromaticity2 = []; % for storing chromaticity obtained from cone contrasts 
    for ii = 1:size(illuminants,1)
        net_spectra_bwl = illuminants(ii,:).*surf_rs;
        [LMS1, ~] = getLMSvals(net_spectra_bwl, net_spectra_bwl, fundamentals);
        chromaticity1 = [chromaticity1; LMSToMacBoyn(LMS1)']; % Converting cone excitations to Macleod-Boynton chromaticity
        chromaticity2 = [chromaticity2; LMSToMacBoyn((LMS1-mean_cone_exc)./mean_cone_exc)']; % Converting cone contrasts to Macleod-Boynton chromaticity
    end    
    figure(14), subplot(121), plot(chromaticity1(:,1), chromaticity1(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',color(jj,:)); hold on;
    subplot(122), plot(chromaticity2(:,1), chromaticity2(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',color(jj,:)); hold on;
end
figure(14),subplot(121);title('Macleod-Boynton Space - exc'), xlabel('r'), ylabel('b'); hold off;
subplot(122);title('Macleod-Boynton Space - contrast'), xlabel('r'), ylabel('b'); hold off;
% Big message - greater the bandwidth, greater is the spread in
% chromaticity over different times of the day

% Now if I will try illuminants which are not the daylight spectra 
gauss_means = 400:10:800;
gauss_spectras_bwl = zeros(numel(gauss_means),numel(wave));
gauss_spectras_bwh = zeros(numel(gauss_means),numel(wave));
gauss_spectras_scaled = zeros(numel(gauss_means),numel(wave));
for ii = 1:numel(gauss_means)
    gauss_spectras_bwl(ii,:) = 10 * normpdf(wave,gauss_means(ii),40); % low bandwidth
    gauss_spectras_bwh(ii,:) = 10 * normpdf(wave,gauss_means(ii),100); % high bandwidth
    gauss_spectras_scaled(ii,:) = ii* normpdf(wave,gauss_means(20),40); % the illumination spectras are scaled version of each other
end
mean_cone_exc = fundamentals'*mean(gauss_spectras_bwl)';
for jj = 1: numel(bandwidths)
    surf_rs= 10 * normpdf(wave,600,bandwidths(jj)); % low bandwidth
    chromaticity1 = []; % for storing chromaticity obtained from cone excitations
    chromaticity2 = []; % for storing chromaticity obtained from cone contrasts 
    for ii = 1:numel(gauss_means)
        net_spectra = gauss_spectras_bwl(ii,:).*surf_rs;
        [LMS1, ~] = getLMSvals(net_spectra, net_spectra, fundamentals);
        chromaticity1 = [chromaticity1; LMSToMacBoyn(LMS1)']; % Converting cone excitations to Macleod-Boynton chromaticity
        chromaticity2 = [chromaticity2; LMSToMacBoyn((LMS1-mean_cone_exc)./mean_cone_exc)']; % Converting cone contrasts to Macleod-Boynton chromaticity
    end
    figure(15),subplot(121), plot(chromaticity1(:,1), chromaticity1(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',color(jj,:)); hold on;
    subplot(122), plot(chromaticity2(:,1), chromaticity2(:,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',color(jj,:)); hold on;
end
figure(15),subplot(121),title('Macleod-Boynton Space - exc'), xlabel('r'), ylabel('b'); hold off;
subplot(122),title('Macleod-Boynton Space - contrast'), xlabel('r'), ylabel('b'); hold off;

%% Trying a new way of creating daylight spectra which I had overlooked before in Judd et al., 1964 paper
wave = 390:5:780;
dayBasis = ieReadSpectra('cieDaylightBasis',wave); % Daylight spectra basis functions from isetbio
num_spectras = 100;
x = linspace(0.25,0.40,num_spectras);
y = 2.870*x - 3.000*(x.*x) - 0.275;
coeff1 = (-1.3515-1.7703*x+5.9114*y)./(0.0241+0.2562*x-0.7341*y);
coeff2 = (0.0300-31.4424*x+30.0717*y)./(0.0241+0.2562*x-0.7341*y);
coeffs = cat(2,ones(num_spectras,1),coeff1',coeff2'); % Limiting the coefficients between 0 and 1
illuminants = coeffs * dayBasis';
figure(16), subplot(131), plot(x,y,'Linewidth',2),xlabel('x'), ylabel('y'); title('CIE daylight'); set(gca,'Xlim',[0.25 0.4]);0
subplot(132),plot(wave,illuminants'), xlabel('Wavelength'), ylabel('Energy'); set(gca,'Xlim',[390 780]); title('Daylight spectra');
subplot(133),plot(coeff1,coeff2,'Linewidth',2), xlabel('M1'),ylabel('M2'); title('Coefficients');

% RF of a double opponent cell
% The subunit vector represents cone weights
RGDO_S1 = [0.3; 0.3; 0.3]; RGDO_S2 = [-0.3; -0.3; -0.3];
hi = size(illuminants,1);
N = 16;
n = sqrt(N);
for jj = 1:N
    rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
    reflectance_spectra1 = munsell(ind,rand_idxs(1));
    reflectance_spectra1 = reflectance_spectra1';
    reflectance_spectra2 = munsell(ind,rand_idxs(2));
    reflectance_spectra2 = reflectance_spectra2';
    Subunit1_act_RGDO_cc = []; Subunit2_act_RGDO_cc = []; % cone contrasts
    Subunit1_act_RGDO_ce = []; Subunit2_act_RGDO_ce = []; % cone excitations
    for ii = 1:hi
        net_spectra1 = illuminants(ii,:).*reflectance_spectra1;
        net_spectra2 = (illuminants(ii,:).*reflectance_spectra2);
        
        L1 = net_spectra1 * fundamentals(:,1);
        M1 = net_spectra1 * fundamentals(:,2);
        S1 = net_spectra1 * fundamentals(:,3);
        L2 = net_spectra2 * fundamentals(:,1);
        M2 = net_spectra2 * fundamentals(:,2);
        S2 = net_spectra2 * fundamentals(:,3);
        
        LMS1 = [L1; M1; S1];
        LMS2 = [L2; M2; S2];
        
        Subunit1_act_RGDO_ce = [Subunit1_act_RGDO_ce; LMS1'*RGDO_S1];
        Subunit2_act_RGDO_ce = [Subunit2_act_RGDO_ce; LMS2'*RGDO_S2];
        % Converting cone excitations into cone contrasts (Weber's
        % contrast)
        mean_cone_exc = fundamentals'*illuminants(ii,:)'; % Mean L,M,S cone excitations
        LMS1 = (LMS1 - mean_cone_exc)./mean_cone_exc; % cone contrasts
        LMS2 = (LMS2 - mean_cone_exc)./mean_cone_exc; % cone contrasts
        
        Subunit1_act_RGDO_cc = [Subunit1_act_RGDO_cc; LMS1'*RGDO_S1];
        Subunit2_act_RGDO_cc = [Subunit2_act_RGDO_cc; LMS2'*RGDO_S2];
        
    end
    figure(17), subplot(n,n,jj), plot(Subunit1_act_RGDO_ce, Subunit2_act_RGDO_ce,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis equal; set(gca,'XTick',[],'YTick',[]); 
    figure(18), subplot(n,n,jj), plot(Subunit1_act_RGDO_cc, Subunit2_act_RGDO_cc,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); axis equal; set(gca,'XTick',[],'YTick',[]); 
    figure(19), subplot(n,n,jj), plot([reflectance_spectra1; reflectance_spectra2]','Linewidth',2); set(gca,'XTick',[],'YTick',[]); 

end