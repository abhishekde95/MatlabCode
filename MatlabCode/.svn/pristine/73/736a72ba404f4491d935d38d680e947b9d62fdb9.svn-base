% Author - Abhishek De
% Writing a new script to see what kind of neurons are suitable to discount
% the effect of illumination

% Creating different daylight spectra from the CIE daylight functions
close all; clearvars;
plot_counter = 1;
wave = 390:5:780;
dayBasis = ieReadSpectra('cieDaylightBasis',wave); % Daylight spectra basis functions from isetbio
num_spectras = 100;
x = linspace(0.25,0.40,num_spectras);
y = 2.870*x - 3.000*(x.*x) - 0.275;
coeff1 = (-1.3515-1.7703*x+5.9114*y)./(0.0241+0.2562*x-0.7341*y);
coeff2 = (0.0300-31.4424*x+30.0717*y)./(0.0241+0.2562*x-0.7341*y);
coeffs = cat(2,ones(num_spectras,1),coeff1',coeff2'); % Limiting the coefficients between 0 and 1
illuminants = coeffs * dayBasis';
figure(plot_counter), subplot(131), plot(x,y,'Linewidth',2),xlabel('x'), ylabel('y'); title('CIE daylight'); set(gca,'Xlim',[0.25 0.4]);
subplot(132),plot(wave,illuminants'), xlabel('Wavelength'), ylabel('Energy'); set(gca,'Xlim',[390 780]); title('Daylight spectra');
subplot(133),plot(coeff1,coeff2,'Linewidth',2), xlabel('M1'),ylabel('M2'); title('Coefficients');
plot_counter = plot_counter + 1;

% Loading the cone action spectra and monitor spectral distributions
load fundamentals.mat
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
fundamentals = fundamentals(3:end,:); % Starting the fundamentals from 390 nm
mon_spd = mon_spd(:,3:end);
% M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations

% One can also use the munsell chips spectra from psychtoolbox stored as
% sur_nickerson.mat. An important thing to note would be that the SPDs are
% measured from 380 to 780 nm in steps of 5nm gap
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
len = numel(380:1:780); % for the munsell reflectance spectras - ftp://ftp.cs.joensuu.fi/pub/color/spectra/mspec/README.txt
ind = 1:5:len;
ind = ind(3:end);

% This should be a optimisation process where I wanna determine what cone
% inputs to a double opponent cell does the best job of discounting the
% effect of illuminant

%% In this part of the simulation, the cone weights to the 2 subunits are constrained to be opposite to mediate spatial opponency,
% In short, spatial opponency is constrained but chromatic opponency within a subunit is not
trials = 150;
figure(plot_counter);
numconditions = 2;
mean_cone_exc = fundamentals'*dayBasis(:,1);
cmap = [1 0 0; 0 1 0; 0 0 1];
figuretitle = ['min var'; 'max var'];
for jj = 1:numconditions
    cone_weights_V1_subunit = [];
    OC_cell = 0; % Orange-Cyan
    LM_cell = 0; % Lime-Magenta
    Lum_cell = 0; % Luminance cell
    BY_cell = 0; % Blue- Yellow Double opponent cell
    for ii = 1:trials
        disp(ii);
        rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
%         rand_idxs = [2 1];
        reflectance_spectra1 = munsell(ind,rand_idxs(1));
        reflectance_spectra1 = reflectance_spectra1';
        reflectance_spectra2 = munsell(ind,rand_idxs(2));
        reflectance_spectra2 = reflectance_spectra2';
        [model,fval,success] = V1cellfit(illuminants,reflectance_spectra1,reflectance_spectra2,fundamentals,mean_cone_exc,jj);
        model = model./norm(model);
        cone_weights_V1_subunit = [cone_weights_V1_subunit; model'];
        L = cone_weights_V1_subunit(end,1);
        M = cone_weights_V1_subunit(end,2);
        S = cone_weights_V1_subunit(end,3);
        if sign(L) == -1*sign(M)
            if sign(S) == sign(L)
                LM_cell = LM_cell + 1;
            elseif sign(S) == sign(M)
                OC_cell = OC_cell + 1;
            end
        elseif sign(L) == sign(M)
            if sign(L) == sign(S)
                Lum_cell = Lum_cell + 1;
            elseif sign(L) == -1*sign(S)
                BY_cell = BY_cell + 1;
            end
        end
    end
    
    subplot(2,2,2*(jj-1)+1),plot3(cone_weights_V1_subunit(:,1),cone_weights_V1_subunit(:,2),cone_weights_V1_subunit(:,3),'o','MarkerSize',5,'LineWidth',0.1,'MarkerFaceColor',cmap(jj,:));
    set(gca,'Xlim',[-1 1],'YLim',[-1 1],'Zlim',[-1 1]); xlabel('L'), ylabel('M'), zlabel('S'); title('Cone Inputs to V1 subunit');
    line([-1 1],[0 0],[0 0],'LineWidth',4); line([0 0],[-1 1],[0 0],'LineWidth',4); line([0 0],[0 0],[-1 1],'LineWidth',4); grid on; hold off;
    subplot(2,2,2*(jj-1)+2),bar([LM_cell, OC_cell, Lum_cell, BY_cell]);set(gca,'XTick',[1 2 3 4],'XTickLabel',{'LM','OC','Lum','BY'}); title(figuretitle(jj,:));
end

plot_counter = plot_counter + 1;
%% Trying a grid search on L-M plane to see the nature of the cost function - based on Greg's suggestion

rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
reflectance_spectra1 = munsell(ind,rand_idxs(1));
reflectance_spectra1 = reflectance_spectra1';
reflectance_spectra2 = munsell(ind,rand_idxs(2));
reflectance_spectra2 = reflectance_spectra2';
L_cone_weight = linspace(-1,1,21);
M_cone_weight = linspace(-1,1,21);
S_cone_weight = linspace(-1,1,21);
resp_variance = zeros(numel(L_cone_weight),numel(M_cone_weight),numel(S_cone_weight));
hi = size(illuminants,1);
min_var_ind = zeros(3,1);
min_variance = 10000;
for ii = 1:numel(L_cone_weight)
    for jj = 1:numel(M_cone_weight)
        for mm = 1:numel(S_cone_weight)
            RGDO_S1 = [L_cone_weight(ii); M_cone_weight(jj); S_cone_weight(mm)];
            RGDO_S2 = -1*RGDO_S1;
            Subunit1_act_cc_constancy = []; Subunit2_act_cc_constancy = [];
            for kk = 1:hi
                net_spectra1 = illuminants(kk,:).*reflectance_spectra1;
                net_spectra2 = illuminants(kk,:).*reflectance_spectra2;
                
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
                mean_cone_exc = fundamentals'*illuminants(ii,:)'; % Mean L,M,S cone excitations
                LMS1 = (LMS1 - mean_cone_exc)./mean_cone_exc; % cone contrasts
                LMS2 = (LMS2 - mean_cone_exc)./mean_cone_exc; % cone contrasts
                
                Subunit1_act_cc_constancy = [Subunit1_act_cc_constancy; LMS1'*RGDO_S1];
                Subunit2_act_cc_constancy = [Subunit2_act_cc_constancy; LMS2'*RGDO_S2];
            end
            Subunit_add_constancy = Subunit1_act_cc_constancy + Subunit2_act_cc_constancy ;
            resp_variance(ii,jj,mm) = var(Subunit_add_constancy) + 2*(1-norm(RGDO_S1)) ;
            if resp_variance(ii,jj,mm) < min_variance
                min_var_ind = [ii;jj;mm];
                min_variance = resp_variance(ii,jj,mm);
            end
        end
    end
end

figure(2),subplot(221),surfl(L_cone_weight,M_cone_weight,squeeze(sum(resp_variance,3))), xlabel('L'), ylabel('M'), zlabel('Var'), title('Cost function');
subplot(222), contour(L_cone_weight,M_cone_weight,squeeze(sum(resp_variance,3)),20), xlabel('L'),ylabel('M'); title('Contour plot');
subplot(223),surfl(L_cone_weight,S_cone_weight,squeeze(sum(resp_variance,2))), xlabel('L'), ylabel('S'), zlabel('Var'), title('Cost function');
subplot(224), contour(L_cone_weight,S_cone_weight,squeeze(sum(resp_variance,2)),20), xlabel('L'),ylabel('S'); title('Contour plot');
% Conclusion: 0,0,0 yield a zero variance

%% Need to come up with a new cost function  such that the solution is not 0,0,0. An ideal double opponent cell should not only be invariant to color edges but
% should also respond differently if the color contrast between the edges is different