% Continuing my research on chromatic and luminance edges from natural
% images.
% Author - Abhishek De, 7/19

% Creating different daylight spectra from the CIE daylight functions
close all; clearvars;
plot_counter = 1;
global illuminants T_xyz fundamentals
wave = 400:10:720; % Taking this bizzare range to match the wavelength ranges of the natural images 
dayBasis = ieReadSpectra('cieDaylightBasis',wave); % Daylight spectra basis functions from isetbio
num_spectras = 100;
x = linspace(0.25,0.40,num_spectras);
y = 2.870*x - 3.000*(x.*x) - 0.275;
coeff1 = (-1.3515-1.7703*x+5.9114*y)./(0.0241+0.2562*x-0.7341*y);
coeff2 = (0.0300-31.4424*x+30.0717*y)./(0.0241+0.2562*x-0.7341*y);
coeffs = cat(2,ones(num_spectras,1),coeff1',coeff2'); % Limiting the coefficients between 0 and 1
illuminants = coeffs * dayBasis';

load T_xyz1964.mat; % 1964 CIE color matching functions
T_xyz = T_xyz1964;
T_xyY = T_xyz./(repmat(sum(T_xyz),[3 1])+0.0001);
rgb = reshape([T_xyY(1,1:end-20)',T_xyY(2,1:end-20)',T_xyY(3,1:end-20)'],[size(T_xyY(:,1:end-20),2) 1 3]);
figure(plot_counter),subplot(231),patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on; plot(0.33,0.34,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot(x,y,'k','Linewidth',2),xlabel('x'), ylabel('y'); title('CIE daylight');  hold off;
subplot(232),plot(wave,illuminants'), xlabel('Wavelength'), ylabel('Energy'); set(gca,'Xlim',[390 780]); title('Daylight spectra');
subplot(233),plot(coeff1,coeff2,'Linewidth',2), xlabel('M1'),ylabel('M2'); title('Coefficients');

% Loading the cone action spectra and moni tor spectral distributions
load fundamentals.mat;
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
lo = find(380:5:780==400);
hi = find(380:5:780==720);
fundamentals = fundamentals(lo:2:hi,:); % Starting the fundamentals from 390 nm
mon_spd = mon_spd(:,lo:2:hi);
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations
figure(plot_counter),subplot(234), plot(wave,fundamentals','Linewidth',2), set(gca,'Xlim',[390 780]); xlabel('Wavelength'); title('fundamentals');
subplot(235), plot(wave,mon_spd','Linewidth',2), set(gca,'Xlim',[390 780]); xlabel('Wavelength'); title('Mon spd');
subplot(236), plot(380:5:780,T_xyz','Linewidth',2), set(gca,'Xlim',[380 780]); xlabel('Wavelength'); title('CIE color matching functions');
plot_counter = plot_counter + 1;


%% Beginning to use a natural hyperspectral image from Nascimento and Foster repository. Giving this a try, don't know as of now where I want to go with this
% but I haven't dealt with natural images before, so might be interesting
load('ref_cyflower1bb_reg1.mat'); % Loading the image file, size is 1019x1337x33
% load('ref_ribeira1bbb_reg1.mat');
% load('ref_ruivaes1bb_reg1.mat');
[R,C,D] = size(reflectances);
reflectances = RGB2XWFormat(reflectances);
idx = [1 100];
im = cell(1,numel(idx));
cexc = cell(1,numel(idx));
imXYZ = cell(1,numel(idx));
imxyY = cell(1,numel(idx));
bkgndxyY = [];
bkgndLMSexc = [];
Lconecon = cell(1,numel(idx)); % Essential for storing the cone contrasts
Mconecon = cell(1,numel(idx));
Sconecon = cell(1,numel(idx));
LMconecon = cell(1,numel(idx)); % for storing L-M cone contrast of the selected pixels
Lumconecon = cell(1,numel(idx));
BYconecon = cell(1,numel(idx));
angulardiffinchromaticities = cell(1,numel(idx));
for ii = 1:numel(idx)
    [imtmp,cexctmp,imXYZtmp] = calcRGB(reflectances,idx(ii),lo,hi,R,C);
    im{ii} = imtmp;
    cexc{ii} = cexctmp;
    imXYZ{ii} = imXYZtmp;
    imxyYtmp = transpose(XYZToxyY(imXYZtmp'));
    imxyY{ii} = imxyYtmp;
    bkgndLMSexc = [bkgndLMSexc mean(cexctmp',2)];
    bkgndxyY = [bkgndxyY XYZToxyY(mean(imXYZtmp',2))];
    Lconecon{ii} = (cexctmp(:,1)-bkgndLMSexc(1,end))/bkgndLMSexc(1,end);
    Mconecon{ii} = (cexctmp(:,2)-bkgndLMSexc(2,end))/bkgndLMSexc(2,end);
    Sconecon{ii} = (cexctmp(:,3)-bkgndLMSexc(3,end))/bkgndLMSexc(3,end);
    LMconecon{ii} = reshape(Lconecon{ii} - Mconecon{ii},R,C); % Useful for analysis later on, storing the L-M cone contrast of the selected pixels
    Lumconecon{ii} = reshape(Lconecon{ii} + Mconecon{ii},R,C);
    BYconecon{ii} = reshape(Sconecon{ii} - 0.707*Lconecon{ii} - 0.707*Mconecon{ii},R,C);
    figure(plot_counter),subplot(2,2,1),image(im{ii}.^0.23), set(gca,'XTick',[],'YTick',[]); axis square;
    subplot(2,2,2),imagesc(edge(LMconecon{ii},'Canny'));set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
    subplot(2,2,3),imagesc(edge(Lumconecon{ii},'Canny'));set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
    subplot(2,2,4),imagesc(edge(BYconecon{ii},'Canny'));set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
    plot_counter = plot_counter + 1;
end


figure(plot_counter); set(gcf,'Name','Illuminants');
plot(wave,illuminants(idx,:)','Linewidth',2); xlabel('Wavelength');ylabel('Energy'); title('Illuminant spectras');
plot_counter = plot_counter + 1;