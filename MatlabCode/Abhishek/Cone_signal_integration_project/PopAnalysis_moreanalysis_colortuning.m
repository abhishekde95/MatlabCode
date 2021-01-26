% A more detailed analysis of the color-tuning
% Author - Abhishek De, 12/20

close all;
clearvars;
plot_counter = plot_counter + 1;

% Section 1: Analysis of subunit and checkberboard color tuning in degrees of angle 

% Section 2: Testing the dependence of color tuning on the number of spikes collected

% Section 3: Checking the relationship between the angle difference and Isoresponse NLI, whitenoise NLI and cone-signal integration NLI

% Section 4: A model simulation for testing the dependence of color tuning on the spatial scale of stimulus 
%            used for probing neuronal tuning

% Section 5: Same as section 4 but for a non-linear cell

%% Section 1:  Analysis of subunit and checkberboard color tuning in degrees of angle 

if ~exist('plot_counter')
    plot_counter = 1;
end

load RGBsubunits.mat
load RGBcheck.mat
load anglediffWNchecksubunit.mat

anglebRFsubfields_Phase1 = [];
anglebRFsubfields_Phase2 = [];
for ii = 1:numel(RGBsubunits)
    vec1 = RGBsubunits{ii}(1:3);
    vec2 = RGBsubunits{ii}(4:end);
    
    vec3 = RGBcheck{ii}(1:3);
    vec4 = RGBcheck{ii}(4:end);
 
    anglebRFsubfields_Phase2 = [anglebRFsubfields_Phase2; 180*acos(dot(vec1,vec2)/(norm(vec1)*norm(vec2)))/pi];
    
    anglebRFsubfields_Phase1 = [anglebRFsubfields_Phase1; 180*acos(dot(vec3,vec4)/(norm(vec3)*norm(vec4)))/pi];

end

bins = 0:10:180;
figure(plot_counter)
subplot(221); histogram(anglediffWNchecksubunit, 0:5:90, 'FaceColor',[0 0 0], 'EdgeColor', [1 1 1]);
xlabel('Angle difference between check and subunit'); ylabel('# cells'); axis square;
set(gca,'Tickdir','out', 'Xlim',[0 90])

subplot(222); histogram(anglebRFsubfields_Phase2, bins, 'FaceColor', [0 0 0], 'EdgeColor', [1 1 1]);
xlabel('Angle between the subunit RGB'); ylabel('# cells'); axis square; set(gca,'Tickdir','out');

subplot(223); histogram(anglebRFsubfields_Phase1, bins, 'FaceColor', [0 0 0], 'EdgeColor', [1 1 1]);
xlabel('Angle between the check RGB'); ylabel('# cells'); axis square; set(gca,'Tickdir','out');

subplot(224); plot(anglebRFsubfields_Phase2, anglediffWNchecksubunit, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1]); hold on; lsline;
xlabel('Angle between the subunit RGB'); ylabel('Angle between the 6-D check and subunit'); axis square; set(gca,'Tickdir','out','Xlim',[0 180],'Ylim',[0 90]);

plot_counter = plot_counter + 1;

% Some stats
[r,p] = corr(anglebRFsubfields_Phase2, anglediffWNchecksubunit,'type','Spearman');

%% Section 2: Testing the dependence of color tuning on the number of spikes collected

if ~exist('plot_counter')
    plot_counter = 1;
end
load nspikes_check.mat
load nspikes_subunit.mat

bins = linspace(min([nspikes_check; nspikes_subunit])-100, max([nspikes_check; nspikes_subunit])+100, 51);

figure(plot_counter);
subplot(121), histogram(nspikes_check,bins,'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth', 2); hold on;
histogram(nspikes_subunit,bins,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth', 2);
xlabel('# spikes'); ylabel('# counts'); axis square; set(gca,'Tickdir','out'); legend('Check','Subunit');  hold off;

subplot(122);
plot(nspikes_check, anglediffWNchecksubunit,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(nspikes_subunit, anglediffWNchecksubunit,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
xlabel('# spikes'); ylabel('Angle difference'); axis square; set(gca,'Tickdir','out','XScale','log'); legend('Check','Subunit'); hold off;


% Some stats on the correlation between the spikes (subunit, checkerboard) and angle diff between check and subunit
[r1, p1] = corr(nspikes_check, anglediffWNchecksubunit,'type','Spearman');
[r2, p2] = corr(nspikes_subunit, anglediffWNchecksubunit, 'type','Spearman');
[r3, p3] = corr(nspikes_check, anglebRFsubfields_Phase1, 'type', 'Spearman');
[r4, p4] = corr(nspikes_subunit, anglebRFsubfields_Phase2, 'type', 'Spearman');

% There doesn't seem to be any dependence the number of spikes.
% This means that the angle difference between the check and subunit is
% real and not an artifact of differences in the spikes collected.

%% Section 3: Checking the relationship between the angle difference and Isoresponse NLI, whitenoise NLI and cone-signal integration NLI

if ~exist('plot_counter')
    plot_counter = 1;
end

% These data are from spatial integration project
load conewts_svd.mat
load vals.mat
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Load the integration across the subunit: white noise analysis data
load AUROClinsubunits_CV.mat
load AUROCquadsubunits_CV.mat

% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Loading the within subunit integration data 
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% For storing the Isoresponse NLI
Isoresponse_NLI = [];

% For storing the white noise NLI
Whitenoise_NLI = [];

% For storing the cone signal integration data
Coneint_NLI = [];

for ii = 1:numel(AUROClinsubunits) 
    
    
    % White noise NLI 
    Error_quad = 1-(AUROCquadsubunits{ii});
    Error_lin = 1-(AUROClinsubunits{ii});
    Whitenoise_NLI = [Whitenoise_NLI; log10(median(Error_lin./Error_quad))];
    
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
    
    % Cone signal integration data 
    Error_linS1 = 1-(AUROClin1{ii}); Error_linS2 = 1-(AUROClin2{ii});
    Error_quadS1 = 1-(AUROCquad1{ii}); Error_quadS2 = 1-(AUROCquad2{ii});
    Coneint_NLI = [Coneint_NLI; mean([log10(median(Error_linS1./Error_quadS1)) log10(median(Error_linS2./Error_quadS2))])];
    
end

% Plotting the results 
figure(plot_counter);
subplot(221); plot(Isoresponse_NLI, anglediffWNchecksubunit, 'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline;
ylabel('Angle deviation'); xlabel('Isoresponse NLI'); set(gca,'Tickdir','out','Ylim',[0 90]); axis square; hold off;

subplot(222); plot(Whitenoise_NLI, anglediffWNchecksubunit, 'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline
ylabel('Angle deviation'); xlabel('White noise NLI'); set(gca,'Tickdir','out','Ylim',[0 90]); axis square; hold off;

subplot(223); plot(Coneint_NLI, anglediffWNchecksubunit, 'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline
ylabel('Angle deviation'); xlabel('Cone signal NLI'); set(gca,'Tickdir','out','Ylim',[0 90]); axis square; hold off;

subplot(224); plot(abs(conewts_svd(3,:)), anglediffWNchecksubunit, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1]); hold on; lsline;
ylabel('Angle deviation'); xlabel('S-cone weight'); set(gca,'Tickdir','out','Ylim',[0 90]); axis square; hold off;

plot_counter = plot_counter + 1;

% Some stats to check the correlation
[r1, p1] = corr(Whitenoise_NLI, anglediffWNchecksubunit, 'type','Spearman');
[r2, p2] = corr(Isoresponse_NLI, anglediffWNchecksubunit, 'type','Spearman'); 
[r3, p3] = corr(Coneint_NLI, anglediffWNchecksubunit, 'type','Spearman');
[r4, p4] = corr(abs(conewts_svd(3,:))', anglediffWNchecksubunit, 'type','Spearman');

% Angle deviation is correlated with all the measures of integration
% Conclusion: Difference in color tuning could be related to non-linearity
% Correlation with S-cone weight almost reached significance

%% Section 4: A model simulation for testing the dependence of color tuning on the spatial scale of stimulus used for probing neuronal tuning

% Step 1: Make a 2-D Gabor filter
% Step 2: Run a white noise RGB noise and collect the response for a linear filter for checkerboard and subunit white noise.
% Step 3: Analysis of color responses obtained using the subunit stim: S1RGB, S2RGB and Allresp

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end

% Step 1: Creating a Gabor filter
theta = pi/2;
Amp = 1;
gamma = 1.5;
sigma = 1.5;
lambda = 1.2;
phi = pi/2;
[X, Y] = meshgrid(-4.5:4.5);
xprime = X.*cos(theta)+Y.*sin(theta);
yprime = -X.*sin(theta)+Y.*cos(theta);
gabor = Amp*exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos((2.*pi.*yprime./lambda)-phi);
gabor = gabor/norm(gabor);

% Creating a mask
thresh = 0.05;
Subunit_mask = gabor;
Subunit_mask(Subunit_mask>thresh) = 1;
Subunit_mask(Subunit_mask<-thresh) = 2;
Subunit_mask(Subunit_mask>=-thresh & Subunit_mask<=thresh) = 0;
S1mask = Subunit_mask == 1;
S2mask = Subunit_mask == 2;

figure(plot_counter);
subplot(321); imagesc(gabor); set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray'); title('Gabor filter')
subplot(322); imagesc(Subunit_mask); set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray'); title('Subunit mask');
subplot(323); histogram(gabor,-0.5:0.05:0.5,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); axis square; set(gca,'Tickdir','out');
xlabel('Pixel intensity'); ylabel('# pixels');
subplot(324); imagesc(gabor.*Subunit_mask); set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray'); title('Masked Gabor');

% Making some changes-eliminating the pixels in the surround for WNGabor and adding RGB 
WNgabor = gabor.*Subunit_mask;
gunweights = [-1 2 0.5];
gunweights = gunweights/sum(abs(gunweights));
WNgabor = cat(3,gunweights(1)*WNgabor,gunweights(2)*WNgabor,gunweights(3)*WNgabor);

% Plotting the WNgabor
im = 0.5 + (0.5*WNgabor/(eps+max(abs(WNgabor(:)))));
subplot(325); bar(gunweights); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Gun Weights');
subplot(326), image(im); axis square; set(gca,'XTick',[],'YTick',[]); title('Colored Gabor');
plot_counter = plot_counter + 1;

% Step 2: Run a white noise RGB noise and collect the response for a linear filter for checkerboard and subunit white noise 
epoch = 100;
stimframes_perepoch = 500;
stimrows = numel(WNgabor(:));

% stim -> WN check 
STCOV_st('init', {numel(gabor(:)) 3 1});
for ii = 1:epoch 
    rgbs = randn(stimrows,stimframes_perepoch);
    response = []; 
    for jj = 1:stimframes_perepoch
        drive = dot(WNgabor(:),rgbs(:,jj));
        response = [response; poissrnd(max([0 drive]).^2)];
    end
    STCOV_st(rgbs(:),response);
end
out = STCOV_st('return');
STA_check = out{1}/out{3};
clear STCOV_st out;

% Calculating the estimated gun weights using the SVD on the WN check
[~,~,v] = svd(reshape(STA_check,[100 3]));
estimated_gunweights = v(:,1)/sum(abs(v(:,1)));
estimated_gunweights = estimated_gunweights *(sign(estimated_gunweights(1))*sign(gunweights(1)));

% stim -> WN subunit 
STCOV_st('init', {numel(gabor(:)) 3 1});
S1RGB = []; S2RGB = []; Allresp = [];
for ii = 1:epoch
    rgbs = [];
    response = []; 
    for jj = 1:stimframes_perepoch
        tmprgb = randn(6,1);
        S1RGB = [S1RGB tmprgb(1:3)]; S2RGB = [S2RGB tmprgb(4:end)];
        tmprgbs = cat(1,S1mask(:)*tmprgb(1),S1mask(:)*tmprgb(2),S1mask(:)*tmprgb(3)) + cat(1,S2mask(:)*tmprgb(4),S2mask(:)*tmprgb(5),S2mask(:)*tmprgb(6));
        rgbs = [rgbs tmprgbs];
        drive = dot(WNgabor(:),tmprgbs);
        response = [response; poissrnd(max([0 drive]).^2)];
    end
    Allresp = [Allresp; response];
    STCOV_st(rgbs(:),response);
end
out = STCOV_st('return');
STA_subunit = out{1}/out{3};
clear STCOV_st out;

% Calculating the estimated gun weights using the SVD on the WN subunit
[~,~,v] = svd(reshape(STA_subunit,[100 3]));
estimated_gunweights_subunit = v(:,1)/sum(abs(v(:,1)));
estimated_gunweights_subunit = estimated_gunweights_subunit *(sign(estimated_gunweights_subunit(1))*sign(gunweights(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLotting the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the STA check 
im1 = 0.5 + (0.5*STA_check/(eps+max(abs(STA_check(:)))));
im1 = reshape(im1, [10 10 3]);

% Plotting the STA subunit 
im2 = 0.5 + (0.5*STA_subunit/(eps+max(abs(STA_subunit(:)))));
im2 = reshape(im2, [10 10 3]);

figure(plot_counter);
subplot(231), image(im); axis square; set(gca,'XTick',[],'YTick',[]); title('Original filter');
subplot(232), image(im1); axis square;  set(gca,'XTick',[],'YTick',[]); title('STA check'); 
subplot(233), image(im2); axis square;  set(gca,'XTick',[],'YTick',[]); title('STA subunit');
subplot(234); bar(gunweights); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Actual Gun Weights');
subplot(235); bar(estimated_gunweights); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Gun Weights check');
subplot(236); bar(estimated_gunweights_subunit); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Gun Weights subunit');
plot_counter = plot_counter + 1;

% Step 3: Analysis of color responses obtained using the subunit stim: S1RGB, S2RGB and Allresp
idx1 = find(S1mask,1); idx1 = [idx1; idx1+100; idx1+200];
idx2 = find(S2mask,1); idx2 = [idx2; idx2+100; idx2+200];
STAS1_RGB = STA_subunit(idx1);
STAS2_RGB = STA_subunit(idx2);

figure(plot_counter);
subplot(321); plot3(S1RGB(1,Allresp>0),S1RGB(2,Allresp>0),S1RGB(3,Allresp>0),'r.'); hold on;
plot3(S1RGB(1,Allresp==0),S1RGB(2,Allresp==0),S1RGB(3,Allresp==0),'k.');
plot3(15*[-STAS1_RGB(1) STAS1_RGB(1)],15*[-STAS1_RGB(2) STAS1_RGB(2)],15*[-STAS1_RGB(3) STAS1_RGB(3)],'g','Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[-5 5],'Ylim',[-5 5],'Zlim',[-5 5]); xlabel('R'); ylabel('G'); zlabel('B'); title('Subunit 1');

subplot(322); plot3(S2RGB(1,Allresp>0),S2RGB(2,Allresp>0),S2RGB(3,Allresp>0),'r.'); hold on;
plot3(S2RGB(1,Allresp==0),S2RGB(2,Allresp==0),S2RGB(3,Allresp==0),'k.'); 
plot3(15*[-STAS2_RGB(1) STAS2_RGB(1)],15*[-STAS2_RGB(2) STAS2_RGB(2)],15*[-STAS2_RGB(3) STAS2_RGB(3)],'g','Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[-5 5],'Ylim',[-5 5],'Zlim',[-5 5]); xlabel('R'); ylabel('G'); zlabel('B'); title('Subunit 2');


% Analysis of color isosurfaces for WN subunit stim  
x = linspace(-3,3,51);
[X,Y,Z] = meshgrid(x,x,x);
covS1spike = cov(S1RGB(:,Allresp>0)');
covS1 = cov(S1RGB');
p1 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS1); p1 = reshape(p1,size(X));
p2 = mvnpdf([X(:) Y(:) Z(:)],STAS1_RGB',covS1spike); p2 = reshape(p2,size(X));
ratio1 = p2./p1;
val1 = prctile(ratio1(:),[25 50 75]);

covS2spike = cov(S2RGB(:,Allresp>0)');
covS2 = cov(S2RGB');
p3 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS2); p3 = reshape(p3,size(X));
p4 = mvnpdf([X(:) Y(:) Z(:)],STAS2_RGB',covS2spike); p4 = reshape(p4,size(X));
ratio2 = p4./p3;
val2 = prctile(ratio2(:),[25 50 75]);

figure(plot_counter);
for jj = 1:3
    fv1 = isosurface(X,Y,Z,ratio1,val1(jj));
    subplot(323); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;

    fv2 = isosurface(X,Y,Z,ratio2,val2(jj));
    subplot(324); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;

end   
subplot(323); plot3(15*[-STAS1_RGB(1) STAS1_RGB(1)],15*[-STAS1_RGB(2) STAS1_RGB(2)],15*[-STAS1_RGB(3) STAS1_RGB(3)],'g','Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[-5 5],'Ylim',[-5 5],'Zlim',[-5 5]); xlabel('R'); ylabel('G'); zlabel('B'); title('Subunit 1'); hold off;

subplot(324); plot3(15*[-STAS2_RGB(1) STAS2_RGB(1)],15*[-STAS2_RGB(2) STAS2_RGB(2)],15*[-STAS2_RGB(3) STAS2_RGB(3)],'g','Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[-5 5],'Ylim',[-5 5],'Zlim',[-5 5]); xlabel('R'); ylabel('G'); zlabel('B'); title('Subunit 2'); hold off;


% Doing some cross-validated GLM/GQM analysis 
responseS1 = Allresp; responseS2 = Allresp;
responseS1(responseS1>0) = 1; responseS2(responseS2>0) = 1; 
responseS1 = logical(responseS1); responseS2 = logical(responseS2);
C1 = cvpartition(responseS1,'KFold',10);
C2 = cvpartition(responseS2,'KFold',10);
AUROClin1 = []; AUROCquad1 = [];
AUROClin2 = []; AUROCquad2 = [];

% Comparing GLM and GQM fits
for jj = 1:C1.NumTestSets
    mdl1lin =  fitglm(S1RGB(:,C1.training(jj))',responseS1(C1.training(jj))','linear','Distribution','binomial','Link','logit');
    mdl1quad =  fitglm(S1RGB(:,C1.training(jj))',responseS1(C1.training(jj))','quadratic','Distribution','binomial','Link','logit');
    mdl2lin =  fitglm(S2RGB(:,C2.training(jj))',responseS2(C2.training(jj))','linear','Distribution','binomial','Link','logit');
    mdl2quad =  fitglm(S2RGB(:,C2.training(jj))',responseS2(C2.training(jj))','quadratic','Distribution','binomial','Link','logit');
    
    % Performing some additional analyses on the model fits (GLM and GQM)
    predlin1 = predict(mdl1lin,S1RGB(:,C1.test(jj))'); % perdiction from GLM subunit 1
    predquad1 = predict(mdl1quad,S1RGB(:,C1.test(jj))'); % perdiction from GQM subunit 1
    predlin2 = predict(mdl2lin,S2RGB(:,C2.test(jj))'); % perdiction from GLM subunit 2
    predquad2 = predict(mdl2quad,S2RGB(:,C2.test(jj))'); % perdiction from GQM subunit 2
    
    % Quantifying accuracy using AUROC
    AUROClin1 = [AUROClin1; rocN(predlin1(responseS1(C1.test(jj))),predlin1(~responseS1(C1.test(jj))))];
    AUROCquad1 = [AUROCquad1; rocN(predquad1(responseS1(C1.test(jj))),predquad1(~responseS1(C1.test(jj))))];
    AUROClin2 = [AUROClin2; rocN(predlin2(responseS2(C2.test(jj))),predlin2(~responseS2(C2.test(jj))))];
    AUROCquad2 = [AUROCquad2; rocN(predquad2(responseS2(C2.test(jj))),predquad2(~responseS2(C2.test(jj))))];
end


subplot(325); plot(AUROClin1, AUROCquad1,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0.5 1],[0.5 1],'k'); axis square; set(gca, 'Tickdir','out', 'Xlim',[0.5 1],'Ylim',[0.5 1]);
xlabel('GLM pred'); ylabel('GQM'); title('Subunit 1'); hold off;

subplot(326); plot(AUROClin2, AUROCquad2,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0.5 1],[0.5 1],'k'); axis square; set(gca, 'Tickdir','out', 'Xlim',[0.5 1],'Ylim',[0.5 1]); 
xlabel('GLM pred'); ylabel('GQM'); title('Subunit 2'); hold off;
plot_counter = plot_counter + 1;

% Checking whether the color tuning estimated using the checkerboard and subunit stim is same 
Color_tuning_checkS1 = [mean(STA_check(find(S1mask(:)==1))); mean(STA_check(find(S1mask(:)==1)+100)); mean(STA_check(find(S1mask(:)==1)+200))];
Color_tuning_checkS2 = [mean(STA_check(find(S2mask(:)==1))); mean(STA_check(find(S2mask(:)==1)+100)); mean(STA_check(find(S2mask(:)==1)+200))];

Color_tuning_check = [Color_tuning_checkS1; Color_tuning_checkS2];
Color_tuning_subunit = [STAS1_RGB; STAS2_RGB];
Anglediffchecksubunit = (180/pi)* acos(dot(Color_tuning_check,Color_tuning_subunit)/(norm(Color_tuning_check)*norm(Color_tuning_subunit)));
AnglediffcheckRGB = (180/pi)* acos(dot(Color_tuning_checkS1,Color_tuning_checkS2)/(norm(Color_tuning_checkS1)*norm(Color_tuning_checkS2)));
AnglediffsubunitRGB = (180/pi)* acos(dot(STAS1_RGB,STAS2_RGB)/(norm(STAS1_RGB)*norm(STAS2_RGB)));


%% % Section 5: Same as section 4 but for a non-linear cell

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end

% Step 1: Creating a Gabor filter, same as the previous section
theta = pi/2;
Amp = 1;
gamma = 1.5;
sigma = 1.5;
lambda = 1.2;
phi = pi/2;
[X, Y] = meshgrid(-4.5:4.5);
xprime = X.*cos(theta)+Y.*sin(theta);
yprime = -X.*sin(theta)+Y.*cos(theta);
gabor = Amp*exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos((2.*pi.*yprime./lambda)-phi);
gabor = gabor/norm(gabor);

% Creating a mask
thresh = 0.05;
Subunit_mask = gabor;
Subunit_mask(Subunit_mask>thresh) = 1;
Subunit_mask(Subunit_mask<-thresh) = 2;
Subunit_mask(Subunit_mask>=-thresh & Subunit_mask<=thresh) = 0;
S1mask = Subunit_mask == 1;
S2mask = Subunit_mask == 2;

% Making some changes-eliminating the pixels in the surround for WNGabor and adding RGB 
WNgabor = gabor.*Subunit_mask;
gunweights = [-1 2 0.5];
gunweights = gunweights/sum(abs(gunweights));
WNgabor = cat(3,gunweights(1)*WNgabor,gunweights(2)*WNgabor,gunweights(3)*WNgabor);

% Step 2: Run a white noise RGB noise and collect the response for a non-linear filter
epoch = 300;
stimframes_perepoch = 500;
stimrows = numel(WNgabor(:));

% stim -> WN check 
STCOV_st('init', {numel(gabor(:)) 3 1});
for ii = 1:epoch 
    rgbs = randn(stimrows,stimframes_perepoch);
    response = []; 
    for jj = 1:stimframes_perepoch
        driveS1 = max(0,dot(WNgabor(1:100)'.*S1mask(:),rgbs(1:100,jj))).^2 + max(0,dot(WNgabor(101:200)'.*S1mask(:),rgbs(101:200,jj))).^2 + max(0,dot(WNgabor(201:300)'.*S1mask(:),rgbs(201:300,jj))).^2;
        driveS2 = max(0,dot(WNgabor(1:100)'.*S2mask(:),rgbs(1:100,jj))).^2 + max(0,dot(WNgabor(101:200)'.*S2mask(:),rgbs(101:200,jj))).^2 + max(0,dot(WNgabor(201:300)'.*S2mask(:),rgbs(201:300,jj))).^2;
        drive = driveS1 + driveS2;
        response = [response; poissrnd(max([0 drive]).^2)];
    end
    STCOV_st(rgbs(:),response);
end
out = STCOV_st('return');
STS = out{1};
nspikes = out{3};
STA_check = out{1}/out{3};
clear STCOV_st out;

% Calculating the estimated gun weights using the SVD on the WN check
[~,~,v] = svd(reshape(STA_check,[100 3]));
estimated_gunweights = v(:,1)/sum(abs(v(:,1)));
estimated_gunweights = estimated_gunweights *(sign(estimated_gunweights(1))*sign(gunweights(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLotting the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the STA check 
im1 = 0.5 + (0.5*STA_check/(eps+max(abs(STA_check(:)))));
im1 = reshape(im1, [10 10 3]);

figure(plot_counter);
im = 0.5 + (0.5*WNgabor/(eps+max(abs(WNgabor(:)))));
subplot(231), image(im); axis square; set(gca,'XTick',[],'YTick',[]); title('Original filter');
subplot(232), image(im1); axis square;  set(gca,'XTick',[],'YTick',[]); title('STA check'); 
subplot(234); bar(gunweights); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Actual Gun Weights');
subplot(235); bar(estimated_gunweights); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Gun Weights check');

% Checking whether the color tuning estimated using the checkerboard and subunit stim is same 
Color_tuning_checkS1 = [mean(STA_check(find(S1mask(:)==1))); mean(STA_check(find(S1mask(:)==1)+100)); mean(STA_check(find(S1mask(:)==1)+200))];
Color_tuning_checkS2 = [mean(STA_check(find(S2mask(:)==1))); mean(STA_check(find(S2mask(:)==1)+100)); mean(STA_check(find(S2mask(:)==1)+200))];

Color_tuning_check = [Color_tuning_checkS1; Color_tuning_checkS2];
AnglediffcheckRGB = (180/pi)* acos(dot(Color_tuning_checkS1,Color_tuning_checkS2)/(norm(Color_tuning_checkS1)*norm(Color_tuning_checkS2)));

% stim -> WN subunit
 
STCOV_st('init', {numel(gabor(:)) 3 1});
S1RGB = []; S2RGB = []; Allresp = [];
for ii = 1:epoch
    rgbs = [];
    response = []; 
    for jj = 1:stimframes_perepoch
        tmprgb = randn(6,1);
        S1RGB = [S1RGB tmprgb(1:3)]; S2RGB = [S2RGB tmprgb(4:end)];
        tmprgbs = cat(1,S1mask(:)*tmprgb(1),S1mask(:)*tmprgb(2),S1mask(:)*tmprgb(3)) + cat(1,S2mask(:)*tmprgb(4),S2mask(:)*tmprgb(5),S2mask(:)*tmprgb(6));
        rgbs = [rgbs tmprgbs];
        
        driveS1 = max(0,dot(WNgabor(1:100)'.*S1mask(:),tmprgbs(1:100))).^2 + max(0,dot(WNgabor(101:200)'.*S1mask(:),tmprgbs(101:200))).^2 + max(0,dot(WNgabor(201:300)'.*S1mask(:),tmprgbs(201:300))).^2;
        driveS2 = max(0,dot(WNgabor(1:100)'.*S2mask(:),tmprgbs(1:100))).^2 + max(0,dot(WNgabor(101:200)'.*S2mask(:),tmprgbs(101:200))).^2 + max(0,dot(WNgabor(201:300)'.*S2mask(:),tmprgbs(201:300))).^2;
        drive = driveS1 + driveS2;

        response = [response; poissrnd(max([0 drive]).^2)];
    end
    Allresp = [Allresp; response];
    STCOV_st(rgbs(:),response);
end
out = STCOV_st('return');
STS = out{1};
STCross = out{2};
nspikes = out{3};
STA_subunit = out{1}/out{3};
clear STCOV_st out;

% Calculating the estimated gun weights using the SVD on the WN subunit
[~,~,v] = svd(reshape(STA_subunit,[100 3]));
estimated_gunweights_subunit = v(:,1)/sum(abs(v(:,1)));
estimated_gunweights_subunit = estimated_gunweights_subunit *(sign(estimated_gunweights_subunit(1))*sign(gunweights(1)));

% Plotting the STA subunit 
im3 = 0.5 + (0.5*STA_subunit/(eps+max(abs(STA_subunit(:)))));
im3 = reshape(im3, [10 10 3]);

figure(plot_counter);
subplot(233), image(im3); axis square;  set(gca,'XTick',[],'YTick',[]); title('STA subunit');
subplot(236); bar(estimated_gunweights_subunit); axis square; set(gca,'Tickdir','out','XTick',[1 2 3],'XTickLabel',{'R','G','B'}); title('Gun Weights subunit');
plot_counter = plot_counter + 1;

% Checking whether the color tuning estimated using the checkerboard and subunit stim is same
idx1 = find(S1mask,1); idx1 = [idx1; idx1+100; idx1+200];
idx2 = find(S2mask,1); idx2 = [idx2; idx2+100; idx2+200];
STAS1_RGB = STA_subunit(idx1);
STAS2_RGB = STA_subunit(idx2);
Color_tuning_checkS1 = [mean(STA_check(find(S1mask(:)==1))); mean(STA_check(find(S1mask(:)==1)+100)); mean(STA_check(find(S1mask(:)==1)+200))];
Color_tuning_checkS2 = [mean(STA_check(find(S2mask(:)==1))); mean(STA_check(find(S2mask(:)==1)+100)); mean(STA_check(find(S2mask(:)==1)+200))];

Color_tuning_check = [Color_tuning_checkS1; Color_tuning_checkS2];
Color_tuning_subunit = [STAS1_RGB; STAS2_RGB];
Anglediffchecksubunit = (180/pi)* acos(dot(Color_tuning_check,Color_tuning_subunit)/(norm(Color_tuning_check)*norm(Color_tuning_subunit)));
AnglediffcheckRGB = (180/pi)* acos(dot(Color_tuning_checkS1,Color_tuning_checkS2)/(norm(Color_tuning_checkS1)*norm(Color_tuning_checkS2)));
AnglediffsubunitRGB = (180/pi)* acos(dot(STAS1_RGB,STAS2_RGB)/(norm(STAS1_RGB)*norm(STAS2_RGB)));

% Analysis of color isosurfaces for WN subunit stim  
x = linspace(-3,3,51);
[X,Y,Z] = meshgrid(x,x,x);
covS1spike = cov(S1RGB(:,Allresp>0)');
covS1 = cov(S1RGB');
p1 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS1); p1 = reshape(p1,size(X));
p2 = mvnpdf([X(:) Y(:) Z(:)],STAS1_RGB',covS1spike); p2 = reshape(p2,size(X));
ratio1 = p2./p1;
val1 = prctile(ratio1(:),[25 50 75]);

covS2spike = cov(S2RGB(:,Allresp>0)');
covS2 = cov(S2RGB');
p3 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS2); p3 = reshape(p3,size(X));
p4 = mvnpdf([X(:) Y(:) Z(:)],STAS2_RGB',covS2spike); p4 = reshape(p4,size(X));
ratio2 = p4./p3;
val2 = prctile(ratio2(:),[25 50 75]);

figure(plot_counter);
for jj = 1:3
    fv1 = isosurface(X,Y,Z,ratio1,val1(jj));
    subplot(221); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;

    fv2 = isosurface(X,Y,Z,ratio2,val2(jj));
    subplot(222); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;

end   
subplot(221); plot3(15*[-STAS1_RGB(1) STAS1_RGB(1)],15*[-STAS1_RGB(2) STAS1_RGB(2)],15*[-STAS1_RGB(3) STAS1_RGB(3)],'g','Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[-5 5],'Ylim',[-5 5],'Zlim',[-5 5]); xlabel('R'); ylabel('G'); zlabel('B'); title('Subunit 1'); hold off;

subplot(222); plot3(15*[-STAS2_RGB(1) STAS2_RGB(1)],15*[-STAS2_RGB(2) STAS2_RGB(2)],15*[-STAS2_RGB(3) STAS2_RGB(3)],'g','Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[-5 5],'Ylim',[-5 5],'Zlim',[-5 5]); xlabel('R'); ylabel('G'); zlabel('B'); title('Subunit 2'); hold off;

% Doing some cross-validated GLM/GQM analysis 
responseS1 = Allresp; responseS2 = Allresp;
responseS1(responseS1>0) = 1; responseS2(responseS2>0) = 1; 
responseS1 = logical(responseS1); responseS2 = logical(responseS2);
C1 = cvpartition(responseS1,'KFold',10);
C2 = cvpartition(responseS2,'KFold',10);
AUROClin1 = []; AUROCquad1 = [];
AUROClin2 = []; AUROCquad2 = [];

% Comparing GLM and GQM fits
for jj = 1:C1.NumTestSets
    mdl1lin =  fitglm(S1RGB(:,C1.training(jj))',responseS1(C1.training(jj))','linear','Distribution','binomial','Link','logit');
    mdl1quad =  fitglm(S1RGB(:,C1.training(jj))',responseS1(C1.training(jj))','quadratic','Distribution','binomial','Link','logit');
    mdl2lin =  fitglm(S2RGB(:,C2.training(jj))',responseS2(C2.training(jj))','linear','Distribution','binomial','Link','logit');
    mdl2quad =  fitglm(S2RGB(:,C2.training(jj))',responseS2(C2.training(jj))','quadratic','Distribution','binomial','Link','logit');
    
    % Performing some additional analyses on the model fits (GLM and GQM)
    predlin1 = predict(mdl1lin,S1RGB(:,C1.test(jj))'); % perdiction from GLM subunit 1
    predquad1 = predict(mdl1quad,S1RGB(:,C1.test(jj))'); % perdiction from GQM subunit 1
    predlin2 = predict(mdl2lin,S2RGB(:,C2.test(jj))'); % perdiction from GLM subunit 2
    predquad2 = predict(mdl2quad,S2RGB(:,C2.test(jj))'); % perdiction from GQM subunit 2
    
    % Quantifying accuracy using AUROC
    AUROClin1 = [AUROClin1; rocN(predlin1(responseS1(C1.test(jj))),predlin1(~responseS1(C1.test(jj))))];
    AUROCquad1 = [AUROCquad1; rocN(predquad1(responseS1(C1.test(jj))),predquad1(~responseS1(C1.test(jj))))];
    AUROClin2 = [AUROClin2; rocN(predlin2(responseS2(C2.test(jj))),predlin2(~responseS2(C2.test(jj))))];
    AUROCquad2 = [AUROCquad2; rocN(predquad2(responseS2(C2.test(jj))),predquad2(~responseS2(C2.test(jj))))];
end


subplot(223); plot(AUROClin1, AUROCquad1,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0.5 1],[0.5 1],'k'); axis square; set(gca, 'Tickdir','out', 'Xlim',[0.5 1],'Ylim',[0.5 1]);
xlabel('GLM pred'); ylabel('GQM'); title('Subunit 1'); hold off;

subplot(224); plot(AUROClin2, AUROCquad2,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0.5 1],[0.5 1],'k'); axis square; set(gca, 'Tickdir','out', 'Xlim',[0.5 1],'Ylim',[0.5 1]); 
xlabel('GLM pred'); ylabel('GQM'); title('Subunit 2'); hold off;
plot_counter = plot_counter + 1;

% GQM prediction is better than GLM prediction but the angle difference is still small