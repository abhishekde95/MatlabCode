% This script looks how a joint gaussian distribution in R and G plane
% looks in L and M plane.
close all;
clear all;
tic;
[stro,filename,no_subunit] = library();

noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
L = stro.trial(:,noisetypeidx)==1;
mu1idx = find(strcmp(stro.sum.trialFields(1,:),'mu1'));
mu2idx = find(strcmp(stro.sum.trialFields(1,:),'mu2'));
sigma1idx = find(strcmp(stro.sum.trialFields(1,:),'sigma1'));
sigma2idx = find(strcmp(stro.sum.trialFields(1,:),'sigma2'));

muvect = unique(stro.trial(L,[mu1idx mu2idx]),'rows')/1000;
sigmavect = unique(stro.trial(L,[sigma1idx sigma2idx]),'rows')/1000;
%%
% Obtaining the A matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd',[380:5:780]); % fitting a cubic spline
A = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals

%%
% Creating a bivariate normal distribution
A = A(1:2,1:2);
A = inv(A');
mu_rg = muvect;
Sigma_rg = [sigmavect(1) 0; 0 sigmavect(2)];
spacing = 0.05;
x1 = -1.5:spacing:1.5; % R 
x2 = -1.5:spacing:1.5; % G
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu_rg,Sigma_rg);
F = reshape(F,length(x2),length(x1));
F = F/max(F(:));
figure(1), subplot(121);
surf(x1,x2,F); axis([-1.5 1.5 -1.5 1.5 0 1.1]);
xlabel('R'); ylabel('G'); zlabel('Probability Density'); title('Gun space');

l = A(1,1)* X1 + A(1,2) * X2;
m = A(2,1)* X1 + A(2,2) * X2;
mu_lm = muvect*A;
Var_l2 = A(1,1).^2*sigmavect(1).^2 + A(1,2).^2*sigmavect(2).^2;
Var_m2 = A(2,1).^2*sigmavect(1).^2 + A(2,2).^2*sigmavect(2).^2;
Var_lm = A(1,1)*A(2,1)*sigmavect(1).^2 + A(1,2)*A(2,2)*sigmavect(2).^2;
Sigma_lm =[sqrt(Var_l2) sqrt(Var_lm); sqrt(Var_lm) sqrt(Var_m2)];

spacing_lm = min([spacing*(A(1,1)+A(1,2)) spacing*(A(2,1)+A(2,2))]);
x3 = min([l(:), m(:)]): spacing_lm :max([l(:), m(:)]); % L
x4 = x3; % M
[L,M] = meshgrid(x3,x4);
F_lm = mvnpdf([L(:) M(:)],mu_lm,Sigma_lm);
F_lm = reshape(F_lm,length(x3),length(x4));
F_lm = F_lm/max(F_lm(:));
figure(1), subplot(122);
surf(x3,x4,F_lm);
xlabel('L'); ylabel('M'); zlabel('Probability Density'); title('Cone space');

%%
% Calculate standard deviation along L+M and L-M direction.
ax = x3 + x4; % L+M
lum_sig = zeros(1,2*numel(ax)-1);
cont_sig = zeros(1,2*numel(ax)-1);
vec = linspace(-1,1,2*numel(ax)-1);
start = 2;
F_lm_flp = flipud(F_lm);
Mat = zeros(2*numel(ax),2*numel(ax)); 
Mat_flp = Mat;
Mat(1:numel(ax),1:numel(ax)) = F_lm;
Mat_flp(1:numel(ax),1:numel(ax)) = F_lm_flp;
while(start<=2*numel(ax))
    n = start-1;
    while(n>=1)
        lum_sig(start-1) = lum_sig(start-1) + Mat(n,start-n);
        cont_sig(start-1) = cont_sig(start-1) + Mat_flp(n,start-n);
        n = n-1;
    end
    start = start + 1;
end

lum_sig = lum_sig/sum(lum_sig);
cont_sig = cont_sig/sum(cont_sig);
lum_std = std(vec.*lum_sig);
cont_std = std(vec.*cont_sig);
lum_2_cont_std_ratio = lum_std/cont_std;
fprintf('The L+M:L-M standard deviation is %1.2d \n',lum_2_cont_std_ratio);

figure(2),plot(lum_sig/sum(lum_sig),'Linewidth',2); hold on
plot(cont_sig/sum(cont_sig),'r','Linewidth',2); legend('L+M','L-M'); hold off;
toc;