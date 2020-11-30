% In this script I am trying to fit the cells to a modified DOG model called the Crescent shaped model 
% Author - Abhishek De, 1/19
% Additing a Gaussian model as well: important for segregating single
% Gaussians vs Difference-of-Gaussians, 4/20, AD
% Introducing a new model called DoGgamma, based on a Jneurophys reviewer's suggestion

close all; clearvars;
load Output_ListWN2.mat
num_rows = 10; % Number of cells in a figure
C = 9;
kk = 1;

% calculating the M matrix
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

N = 100;
bins = -100:10:300;
numcrescentparams = 8;
numgaborparams = 8;
numDOGparams = 6;
numSGparams = 4;
numS2Gparams = 6;
numDOGgammaparams = 7;
crit = chi2inv(0.9999,300);
Crescenterror = [];
Gaborerror = [];
DOGerror = [];
SGerror = [];
S2Gerror = [];
DOGgammaerror = [];

% Pruning out the rows in the cell where signal to noise ratio is low
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];

% After storing the parameters, I can use it to classify single-opponent vs non-single opponent neurons
numcells = size(Output_List,1);
Gaborparams = cell(numcells,1);
DOGparams = cell(numcells,1);
SGparams = cell(numcells,1);
S2Gparams = cell(numcells,1);
Crescentparams = cell(numcells,1);
DOGgammaparams = cell(numcells,1);
Singleopponent = zeros(numcells,1);
modelfits = cell(numcells,5);
ori_tuning = cell(numcells,3);
sf_tuning = cell(numcells,3);

RGBsinglepixel = [];
RGB_svd = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
R_square = []; % Storing coefficient of determination
Pearson_R = []; % Pearson's correlation coefficient
lambda = [20 10 5 10/3 10/4 10/5];
theta = pi*(0:90:360)/180; theta(end) = [];
plotresults = 0;
for ii = 1:numcells
    disp(ii);
    [out2,fittedGabor,Rsquare_Gabor,Pearson_RGabor] = gaborfit_AD(Output_List{ii,4}); % fitting Gabor
    [out4,fittedSGaussian,Rsquare_SG,Pearson_RSG] = SingleGaussianfit(Output_List{ii,4}); % fitting the Single_Gaussian
    [out5,fittedS2Gaussian,Rsquare_S2G,Pearson_RS2G] = SingleGaussian2Dfit(Output_List{ii,4}); % fitting the Single_Gaussian with variable variances 
    [out3,fittedDOG,Rsquare_DOG,Pearson_RDOG] = DOGfit(Output_List{ii,4},[]); % Fitting the DoG
    
    [tmp_out1,fittedDOGwgamma1,Rsquare_DOGwgamma1,Pearson_RDOGwgamma1] = DOGfit_wgamma(Output_List{ii,4},out3); % Fitting the DoGwgamma
    [tmp_out2,fittedDOGwgamma2,Rsquare_DOGwgamma2,Pearson_RDOGwgamma2] = DOGfit_wgamma(Output_List{ii,4},[]);
     if tmp_out1.fval < tmp_out2.fval
        fittedDOGwgamma = fittedDOGwgamma1;
        out6 = tmp_out1;
        Rsquare_DOGwgamma = Rsquare_DOGwgamma1;
        Pearson_RDOGwgamma = Pearson_RDOGwgamma1;
    else
        fittedDOGwgamma = fittedDOGwgamma2;
        out6 = tmp_out2;
        Rsquare_DOGwgamma = Rsquare_DOGwgamma2;
        Pearson_RDOGwgamma = Pearson_RDOGwgamma2;
    end
    
    
    [tmp_out1,fittedCrescent1,Rsquare_Crescent1,Pearson_RCrescent1] = Crescentfit_AD(Output_List{ii,4},out3); % Fitting the crescent shaped model, takes input from the DOG model
    [tmp_out2,fittedCrescent2,Rsquare_Crescent2,Pearson_RCrescent2] = Crescentfit_AD(Output_List{ii,4},[]);
    if tmp_out1.fval < tmp_out2.fval
        fittedCrescent = fittedCrescent1;
        out1 = tmp_out1;
        Rsquare_Crescent = Rsquare_Crescent1;
        Pearson_RCrescent = Pearson_RCrescent1;
    else
        fittedCrescent = fittedCrescent2;
        out1 = tmp_out2;
        Rsquare_Crescent = Rsquare_Crescent2;
        Pearson_RCrescent = Pearson_RCrescent2;
    end
    
    modelfits{ii,1} = fittedCrescent;
    modelfits{ii,2} = fittedGabor;
    modelfits{ii,3} = fittedDOG;
    modelfits{ii,4} = fittedSGaussian;
    modelfits{ii,5} = fittedS2Gaussian;
    modelfits{ii,6} = fittedDOGwgamma;
    [ori_tuning1, sf_tuning1] = resp_to_sinusoids(out1,lambda,theta,1); % Crescent
    [ori_tuning2, sf_tuning2] = resp_to_sinusoids(out2,lambda,theta,2); % Gabor
    [ori_tuning3, sf_tuning3] = resp_to_sinusoids(out3,lambda,theta,3); % DOG
    
    ori_tuning{ii,1} = ori_tuning1; ori_tuning{ii,2} = ori_tuning2; ori_tuning{ii,3} = ori_tuning3;
    sf_tuning{ii,1} = ori_tuning1; sf_tuning{ii,2} = ori_tuning2; sf_tuning{ii,3} = sf_tuning3; 
    
    Crescenterror = [Crescenterror; out1.fval];
    Gaborerror = [Gaborerror; out2.fval];
    DOGerror = [DOGerror; out3.fval];
    SGerror = [SGerror; out4.fval];
    S2Gerror = [S2Gerror; out5.fval];
    DOGgammaerror = [DOGgammaerror; out6.fval];
    
    tmpCrescentBIC = N*log(Crescenterror(end)/N) + numcrescentparams*log(N);
    tmpGaborBIC = N*log(Gaborerror(end)/N) + numgaborparams*log(N);
    tmpDOGBIC = N*log(DOGerror(end)/N) + numDOGparams*log(N);
    tmpSGBIC = N*log(SGerror(end)/N) + numSGparams*log(N);
    tmpS2GBIC = N*log(S2Gerror(end)/N) + numS2Gparams*log(N);
    tmpDOGBIC = N*log(DOGgammaerror(end)/N) + numDOGgammaparams*log(N);
    [~,minidx] = min([tmpCrescentBIC tmpGaborBIC tmpDOGBIC tmpSGBIC tmpS2GBIC]);
    R_square = [R_square; Rsquare_Crescent Rsquare_Gabor Rsquare_DOG Rsquare_SG Rsquare_S2G Rsquare_DOGwgamma];
    Pearson_R = [Pearson_R; Pearson_RCrescent Pearson_RGabor Pearson_RDOG Pearson_RSG Pearson_RS2G Pearson_RDOGwgamma];
    
    % Storing the model parameters
    Crescentparams{ii} = out1;
    Gaborparams{ii} = out2;
    DOGparams{ii} = out3;
    SGparams{ii} = out4;
    S2Gparams{ii} = out5;
    DOGgammaparams{ii} = out6;
    
    % Storing the RGB of the pixel with max energy
    RGBallpixels = Output_List{ii,2};
    [~,ind] = max(sum(RGBallpixels.^2,2));
    RGBsinglepixel = [RGBsinglepixel; RGBallpixels(ind,:)];
    
    % Storing the RGB calculated using SVD
    RGB_svd = [RGB_svd; Output_List{ii,5}'];
        
    % Checking if single-opponent or not
    ft = fft2(Output_List{ii,4});
    if max(abs(ft(:)))== abs(ft(1))
        Singleopponent(ii) = 1;
    else
        Singleopponent(ii) = 0;
    end
    
    % Calculating unexplained variance according to Ringach et al., 2002  
%     keyboard;
%     im = Output_List{ii,4};
%     factor = max(abs(im(:)));
%     noiseimage = Output_List{ii,14};
%     unexplainedvar = [unexplainedvar; var(im(:)-factor*fittedCrescent(:))/var(im(:)) var(im(:)-factor*fittedGabor(:))/var(im(:)) var(im(:)-factor*fittedDOG(:))/var(im(:))];

if plotresults
    figure(kk); colormap('gray');
    % STA
    tmp_vec_gun = Output_List{ii,2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+1); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1
        title('STA');
    end
    
    % Peak Frame
    im = Output_List{ii,4};
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+2); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    if mod(ii,num_rows) ==1
        title('RF');
    end
    
    % Crescent fit to the RF
    im = fittedCrescent;
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+3); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    if minidx == 1
        set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
    end
    if mod(ii,num_rows) ==1
        title('Crescent');
    end
    
    % Gabor fit to the RF
    im = fittedGabor;
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+4); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    if minidx == 2
        set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
    end
    if mod(ii,num_rows) ==1
        title('Gabor');
    end
    
    % DOG fit to the RF
    im = fittedDOG;
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+5); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    
    if minidx == 3
        set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
    end
    if mod(ii,num_rows) ==1
        title('DOG');
    end
    
    % SG fit to the RF
    im = fittedSGaussian;
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+6); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    
    if minidx == 4
        set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
    end
    if mod(ii,num_rows) ==1
        title('SG');
    end
    
    % S2G fit to the RF
    im = fittedS2Gaussian;
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+7); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    
    if minidx == 5
        set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
    end
    if mod(ii,num_rows) ==1
        title('S2G');
    end
    
    % DOG w gamma fit to the RF
    im = fittedDOGwgamma;
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+8); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    
    if minidx == 6
        set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
    end
    if mod(ii,num_rows) ==1
        title('DOG gamma');
    end
    
    % Energy
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+9); plot(Output_List{ii,7}); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1
        title('Energy');
    end
    
    if mod(ii,num_rows) == 0
        kk = kk + 1;
    end
end
end
plot_counter = kk + 1;


savevariables = 1;
if savevariables == 1
    save Crescentparams Crescentparams
    save Gaborparams Gaborparams
    save DOGparams DOGparams
    save SGparams SGparams
    save S2Gparams S2Gparams
    save DOGgammaparams DOGgammaparams
    save Singleopponent Singleopponent
    save Crescenterror Crescenterror
    save Gaborerror Gaborerror
    save DOGerror DOGerror
    save SGerror SGerror
    save S2Gerror S2Gerror
    save DOGgammaerror DOGgammaerror
    save modelfits modelfits
    save R_square R_square
    save ori_tuning ori_tuning
    save sf_tuning sf_tuning
    save Pearson_R Pearson_R
end

%%
CrescentBIC = N*log(Crescenterror(~Singleopponent & simplecells,:)/N) + numcrescentparams*log(N);
GaborBIC = N*log(Gaborerror(~Singleopponent & simplecells,:)/N) + numgaborparams*log(N);
DOGBIC = N*log(DOGerror(~Singleopponent & simplecells,:)/N) + numDOGparams*log(N);
BICs = [CrescentBIC GaborBIC DOGBIC];
[minval,I] = min(BICs,[],2);
medianBICval = median(BICs - repmat(minval,[1 3]),2);
RF = cell2mat(Output_List(~Singleopponent & simplecells,8));

% storing gabor phases of the cells
gaborphases = zeros(1,numcells);
aspectratio = zeros(1,numcells);
oribandwidth = zeros(1,numcells); % in octaves
for ii = 1:numcells
    gaborphases(ii) = rem(Gaborparams{ii}.phi,pi)*180/pi;
    aspectratio(ii) = Gaborparams{ii}.gamma;
end
Singleopponent = logical(Singleopponent);
gaborphases_relevantcells = gaborphases(~Singleopponent & simplecells);
aspectratio_relevantcells = aspectratio(~Singleopponent & simplecells);

bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
% Calculating conewts of the cells of interest: SVD 

conewts_svd = Mrgbtocc * RGB_svd(~Singleopponent & simplecells,:)';
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% Calculating conewts of the cells of interest: single pixel
conewts_singlepixel = Mrgbtocc * RGBsinglepixel(~Singleopponent & simplecells,:)';
conewts_singlepixel = conewts_singlepixel./repmat(sum(abs(conewts_singlepixel),1),[3 1]);
conewts_singlepixel = conewts_singlepixel .* repmat(sign(conewts_singlepixel(2,:)),[3 1]);

% Calculating conewts of the SO cells: SVD 
conewts_svdSO = Mrgbtocc * RGB_svd(Singleopponent & simplecells,:)';
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

% Calculating conewts of the SO cells: single pixel
conewts_singlepixelSO = Mrgbtocc * RGBsinglepixel(Singleopponent & simplecells,:)';
conewts_singlepixelSO = conewts_singlepixelSO./repmat(sum(abs(conewts_singlepixelSO),1),[3 1]);
conewts_singlepixelSO = conewts_singlepixelSO .* repmat(sign(conewts_singlepixelSO(2,:)),[3 1]);

% Just visualizing the cone weights of all cells: single opponent cells & luminance and color cells
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconedominated_conewtsLM = find(abs(conewts_svd(3,:))>1-thresh & conewts_svd(3,:)<0);
Sconedominated_conewtsOC = find(abs(conewts_svd(3,:))>1-thresh & conewts_svd(3,:)>0);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
LumIds_conewtssinglepixel = find(conewts_singlepixel(1,:) + conewts_singlepixel(2,:) >thresh & sum(sign(conewts_singlepixel(1:2,:)),1)==2 & conewts_singlepixel(1,:)>0.1 & conewts_singlepixel(2,:)>0.1);
ColorOpponentIds_conewtssinglepixel = find(conewts_singlepixel(2,:) - conewts_singlepixel(1,:) >thresh & sum(sign(conewts_singlepixel(1:2,:)),1)==0 & sqrt((conewts_singlepixel(2,:)-0.5).^2 + (conewts_singlepixel(1,:)+0.5).^2)<0.3);
Sconedominated_conewtssinglepixel = find(abs(conewts_singlepixel(3,:))>1-thresh);
Other_conewtssinglepixel = 1:size(conewts_singlepixel,2); Other_conewtssinglepixel([LumIds_conewtssinglepixel ColorOpponentIds_conewtssinglepixel Sconedominated_conewtssinglepixel]) = [];

figure(plot_counter); set(gcf,'Name','Cone wts') 
subplot(321); plot(conewts_svd(1,:),conewts_svd(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svdSO(1,:),conewts_svdSO(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); title('Cone wts: SVD'); axis equal; xlabel('L'), ylabel('M');
subplot(322); plot(conewts_singlepixel(1,:),conewts_singlepixel(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; 
plot(conewts_singlepixelSO(1,:),conewts_singlepixelSO(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); title('Cone wts: single pixel'); axis equal; xlabel('L'), ylabel('M');
subplot(323); plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis equal; xlabel('L'), ylabel('M');
subplot(324); plot(conewts_singlepixel(1,LumIds_conewtssinglepixel),conewts_singlepixel(2,LumIds_conewtssinglepixel),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_singlepixel(1,ColorOpponentIds_conewtssinglepixel),conewts_singlepixel(2,ColorOpponentIds_conewtssinglepixel),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_singlepixel(1,Sconedominated_conewtssinglepixel),conewts_singlepixel(2,Sconedominated_conewtssinglepixel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_singlepixel(1,Other_conewtssinglepixel),conewts_singlepixel(2,Other_conewtssinglepixel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis equal; xlabel('L'), ylabel('M');
subplot(325); plot(conewts_svd(3,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(3,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(3,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(3,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis equal; xlabel('S'), ylabel('M');
subplot(326); plot(conewts_singlepixel(3,LumIds_conewtssinglepixel),conewts_singlepixel(2,LumIds_conewtssinglepixel),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_singlepixel(3,ColorOpponentIds_conewtssinglepixel),conewts_singlepixel(2,ColorOpponentIds_conewtssinglepixel),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_singlepixel(3,Sconedominated_conewtssinglepixel),conewts_singlepixel(2,Sconedominated_conewtssinglepixel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_singlepixel(3,Other_conewtssinglepixel),conewts_singlepixel(2,Other_conewtssinglepixel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis equal; xlabel('S'), ylabel('M');
plot_counter = plot_counter + 1;

% plotting the 3-D view of the cone weights 
figure(plot_counter); set(gcf,'Name','Cone wts: 3D view');
subplot(121); plot3(conewts_svd(1,LumIds_conewts),conewts_svd(3,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot3(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(3,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot3(conewts_svd(1,Sconedominated_conewts),conewts_svd(3,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on;
plot3(conewts_svd(1,Other_conewts),conewts_svd(3,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis equal; xlabel('L'); ylabel('S'); zlabel('M');
subplot(122); plot3(conewts_singlepixel(1,LumIds_conewtssinglepixel),conewts_singlepixel(3,LumIds_conewtssinglepixel),conewts_singlepixel(2,LumIds_conewtssinglepixel),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot3(conewts_singlepixel(1,ColorOpponentIds_conewtssinglepixel),conewts_singlepixel(3,ColorOpponentIds_conewtssinglepixel),conewts_singlepixel(2,ColorOpponentIds_conewtssinglepixel),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot3(conewts_singlepixel(1,Sconedominated_conewtssinglepixel),conewts_singlepixel(3,Sconedominated_conewtssinglepixel),conewts_singlepixel(2,Sconedominated_conewtssinglepixel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); hold on; 
plot3(conewts_singlepixel(1,Other_conewtssinglepixel),conewts_singlepixel(3,Other_conewtssinglepixel),conewts_singlepixel(2,Other_conewtssinglepixel),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis equal; xlabel('L'); ylabel('S'); zlabel('M');
plot_counter = plot_counter + 1;

% plotting the R square value for all the model fits: Gabor, Crescent & DoG
bins = 0:0.05:1;
figure(plot_counter); set(gcf,'Name','Quality of fits');
subplot(131),histogram(R_square(~Singleopponent & simplecells,1),bins,'Normalization','probability'); xlabel('R^2'); ylabel('proportion'); title('Crescent'); set(gca,'Tickdir','out','Ylim',[0 1]); axis square; 
subplot(132),histogram(R_square(~Singleopponent & simplecells,2),bins,'Normalization','probability'); xlabel('R^2'); ylabel('proportion'); title('Gabor'); set(gca,'Tickdir','out','Ylim',[0 1]); axis square; 
subplot(133),histogram(R_square(~Singleopponent & simplecells,3),bins,'Normalization','probability'); xlabel('R^2'); ylabel('proportion'); title('DoG'); set(gca,'Tickdir','out','Ylim',[0 1]); axis square; 
plot_counter = plot_counter + 1;

% Significance test for bimodality (Hartigan's dip test) 
[~, p_Lumsvd] = HartigansDipSignifTest(gaborphases_relevantcells(LumIds_conewts(I(LumIds_conewts)==2)), 500);
[~, p_COsvd] = HartigansDipSignifTest(gaborphases_relevantcells(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2)), 500);
[~, p_Ssvd] = HartigansDipSignifTest(gaborphases_relevantcells(Sconedominated_conewts(I(Sconedominated_conewts)==2)), 500);
[~, p_Othersvd] = HartigansDipSignifTest(gaborphases_relevantcells(Other_conewts(I(Other_conewts)==2)), 500);
[~, p_Lumsinglepixel] = HartigansDipSignifTest(gaborphases_relevantcells(LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)==2)), 500);
[~, p_COsinglepixel] = HartigansDipSignifTest(gaborphases_relevantcells(ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)==2)), 500);
[~, p_Ssinglepixel] = HartigansDipSignifTest(gaborphases_relevantcells(Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)==2)), 500);
[~, p_Othersinglepixel] = HartigansDipSignifTest(gaborphases_relevantcells(Other_conewtssinglepixel(I(Other_conewtssinglepixel)==2)), 500);

% Plotting phases from gabor fits for Luminance and CO cells 
figure(plot_counter); set(gcf,'Name','Phases from Gabor fit');  
subplot(421); histogram(gaborphases_relevantcells(LumIds_conewts(I(LumIds_conewts)==2)),0:10:180,'Normalization','probability','FaceColor',[0 1 0]); set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('SVD:lum',num2str(p_Lumsvd,2)));
subplot(422); histogram(gaborphases_relevantcells(LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)==2)),0:10:180,'Normalization','probability','FaceColor',[0 1 0]);  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:lum',num2str(p_Lumsinglepixel,2)));
subplot(423); histogram(gaborphases_relevantcells(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2)),0:10:180,'Normalization','probability','FaceColor',[1 0.5 0]);  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('SVD:CO',num2str(p_COsvd,2)));
subplot(424); histogram(gaborphases_relevantcells(ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)==2)),0:10:180,'Normalization','probability','FaceColor',[1 0.5 0]);  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:CO',num2str(p_COsinglepixel,2)));
subplot(425); histogram(gaborphases_relevantcells(Sconedominated_conewts(I(Sconedominated_conewts)==2)),0:10:180,'Normalization','probability','FaceColor',[0 0.5 1.0]); xlabel('phase (degrees)');  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); ylabel('count'); title(strcat('SVD:Other',num2str(p_Ssvd,2)));
subplot(426); histogram(gaborphases_relevantcells(Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)==2)),0:10:180,'Normalization','probability','FaceColor',[0 0.5 1.0]);  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:Other',num2str(p_Ssinglepixel,2)));
subplot(427); histogram(gaborphases_relevantcells(Other_conewts(I(Other_conewts)==2)),0:10:180,'Normalization','probability','FaceColor',[0 0 0]);  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('SVD:Other',num2str(p_Othersvd,2)));
subplot(428); histogram(gaborphases_relevantcells(Other_conewtssinglepixel(I(Other_conewtssinglepixel)==2)),0:10:180,'Normalization','probability','FaceColor',[0 0 0]);  set(gca,'Tickdir','out','Xlim',[0 180],'XTick',0:90:180); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:Other',num2str(p_Othersinglepixel,2)));
plot_counter = plot_counter + 1;

% Another representation of plotting phases from gabor fits for Luminance and CO cells
figure(plot_counter); set(gcf,'Name','Phases from Gabor fit: Representation adapted from Ringach et al., 2002');  
subplot(421); histogram(90-abs(90-gaborphases_relevantcells(LumIds_conewts(I(LumIds_conewts)==2))),0:10:90,'Normalization','probability','FaceColor',[0 1 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('SVD:lum',num2str(p_Lumsvd,2)));
subplot(422); histogram(90-abs(90-gaborphases_relevantcells(LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)==2))),0:10:90,'Normalization','probability','FaceColor',[0 1 0]);  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:lum',num2str(p_Lumsinglepixel,2)));
subplot(423); histogram(90-abs(90-gaborphases_relevantcells(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2))),0:10:90,'Normalization','probability','FaceColor',[1 0.5 0]);  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('SVD:CO',num2str(p_COsvd,2)));
subplot(424); histogram(90-abs(90-gaborphases_relevantcells(ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)==2))),0:10:90,'Normalization','probability','FaceColor',[1 0.5 0]);  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:CO',num2str(p_COsinglepixel,2)));
subplot(425); histogram(90-abs(90-gaborphases_relevantcells(Sconedominated_conewts(I(Sconedominated_conewts)==2))),0:10:90,'Normalization','probability','FaceColor',[0 0.5 1.0]); xlabel('phase (degrees)');  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); ylabel('count'); title(strcat('SVD:Other',num2str(p_Ssvd,2)));
subplot(426); histogram(90-abs(90-gaborphases_relevantcells(Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)==2))),0:10:90,'Normalization','probability','FaceColor',[0 0.5 1.0]);  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:Other',num2str(p_Ssinglepixel,2)));
subplot(427); histogram(90-abs(90-gaborphases_relevantcells(Other_conewts(I(Other_conewts)==2))),0:10:90,'Normalization','probability','FaceColor',[0 0 0]);  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('SVD:Other',num2str(p_Othersvd,2)));
subplot(428); histogram(90-abs(90-gaborphases_relevantcells(Other_conewtssinglepixel(I(Other_conewtssinglepixel)==2))),0:10:90,'Normalization','probability','FaceColor',[0 0 0]);  set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90); xlabel('phase (degrees)'); ylabel('count'); title(strcat('single-pixel:Other',num2str(p_Othersinglepixel,2)));
plot_counter = plot_counter + 1;

% Plotting aspect ratios for different cells
bins = 0:0.25:5;
figure(plot_counter); set(gcf,'Name','Aspect ratios');
subplot(221); histogram(aspectratio_relevantcells(LumIds_conewts(I(LumIds_conewts)==2)),bins,'Normalization','probability','FaceColor',[0 1 0]);hold on; plot(median(aspectratio_relevantcells(LumIds_conewts(I(LumIds_conewts)==2))),0,'kv'); xlabel('Aspect ratio'); ylabel('Count'); set(gca,'Tickdir','out'); title('Lum'); 
subplot(222); histogram(aspectratio_relevantcells(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2)),bins,'Normalization','probability','FaceColor',[1 0.5 0]); hold on; plot(median(aspectratio_relevantcells(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2))),0,'kv'); xlabel('Aspect ratio'); ylabel('Count'); set(gca,'Tickdir','out'); title('L-M');
subplot(223); histogram(aspectratio_relevantcells(Sconedominated_conewts(I(Sconedominated_conewts)==2)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; plot(median(aspectratio_relevantcells(Sconedominated_conewts(I(Sconedominated_conewts)==2))),0,'kv'); xlabel('Aspect ratio'); ylabel('Count'); set(gca,'Tickdir','out'); title('S');
subplot(224); histogram(aspectratio_relevantcells(Other_conewts(I(Other_conewts)==2)),bins,'Normalization','probability','FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Other_conewts(I(Other_conewts)==2))),0,'kv'); xlabel('Aspect ratio'); ylabel('Count'); set(gca,'Tickdir','out'); title('Other');
plot_counter = plot_counter + 1;

% Plotting phases from gabor fits for Luminance and CO cells as a function of RF eccentricities 
figure(plot_counter); set(gcf,'Name','Phases from Gabor fit vs eccentricity');  
subplot(421); plot(sqrt(sum(RF(LumIds_conewts(I(LumIds_conewts)==2),:).^2,2)),gaborphases_relevantcells(LumIds_conewts(I(LumIds_conewts)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('SVD:lum');
subplot(422); plot(sqrt(sum(RF(LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)==2),:).^2,2)),gaborphases_relevantcells(LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('single-pixel:lum');
subplot(423); plot(sqrt(sum(RF(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2),:).^2,2)),gaborphases_relevantcells(ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('SVD:CO');
subplot(424); plot(sqrt(sum(RF(ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)==2),:).^2,2)),gaborphases_relevantcells(ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('single-pixel:CO');
subplot(425); plot(sqrt(sum(RF(Sconedominated_conewts(I(Sconedominated_conewts)==2),:).^2,2)),gaborphases_relevantcells(Sconedominated_conewts(I(Sconedominated_conewts)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('SVD:S');
subplot(426); plot(sqrt(sum(RF(Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)==2),:).^2,2)),gaborphases_relevantcells(Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('single-pixel:S');
subplot(427); plot(sqrt(sum(RF(Other_conewts(I(Other_conewts)==2),:).^2,2)),gaborphases_relevantcells(Other_conewts(I(Other_conewts)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('SVD:other');
subplot(428); plot(sqrt(sum(RF(Other_conewtssinglepixel(I(Other_conewtssinglepixel)==2),:).^2,2)),gaborphases_relevantcells(Other_conewtssinglepixel(I(Other_conewtssinglepixel)==2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); lsline; xlabel('ecc'); ylabel('phase'); title('single-pixel:other');
plot_counter = plot_counter + 1;

% Plotting the distribution of best model fits for luminance, CO and other cells classified based on 2 criteria: SVD and single pixel
figure(plot_counter); set(gcf,'Name','Distribution of dominant model: Crescent, Gabor and DoG');
I1 = I([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts Other_conewts]);
subplot(321); bar([numel(find(I1==1)) numel(find(I1==2)) numel(find(I1==3))]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells'); 
subplot(322); bar([numel(find(I1==1))/numel(I1) numel(find(I1==2))/numel(I1) numel(find(I1==3))/numel(I1)]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells');
subplot(323); bar([sum(I(LumIds_conewts)==1) sum(I(ColorOpponentIds_conewts)==1) sum(I(Sconedominated_conewts)==1) sum(I(Other_conewts)==1); sum(I(LumIds_conewts)==2) sum(I(ColorOpponentIds_conewts)==2) sum(I(Sconedominated_conewts)==2) sum(I(Other_conewts)==2);  sum(I(LumIds_conewts)==3) sum(I(ColorOpponentIds_conewts)==3) sum(I(Sconedominated_conewts)==3) sum(I(Other_conewts)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); legend({'Lum','CO','S','Other'},'FontSize',8); ylabel('No. of cells');  title('SVD');
subplot(324); bar([sum(I(LumIds_conewts)==1)/numel(LumIds_conewts) sum(I(ColorOpponentIds_conewts)==1)/numel(ColorOpponentIds_conewts) sum(I(Sconedominated_conewts)==1)/numel(Sconedominated_conewts) sum(I(Other_conewts)==1)/numel(Other_conewts); sum(I(LumIds_conewts)==2)/numel(LumIds_conewts) sum(I(ColorOpponentIds_conewts)==2)/numel(ColorOpponentIds_conewts) sum(I(Sconedominated_conewts)==2)/numel(Sconedominated_conewts) sum(I(Other_conewts)==2)/numel(Other_conewts);  sum(I(LumIds_conewts)==3)/numel(LumIds_conewts) sum(I(Sconedominated_conewts)==3)/numel(Sconedominated_conewts) sum(I(ColorOpponentIds_conewts)==3)/numel(ColorOpponentIds_conewts) sum(I(Other_conewts)==3)/numel(Other_conewts)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('SVD');
subplot(325); bar([sum(I(LumIds_conewtssinglepixel)==1) sum(I(ColorOpponentIds_conewtssinglepixel)==1) sum(I(Sconedominated_conewtssinglepixel)==1) sum(I(Other_conewtssinglepixel)==1); sum(I(LumIds_conewtssinglepixel)==2) sum(I(ColorOpponentIds_conewtssinglepixel)==2) sum(I(Sconedominated_conewtssinglepixel)==2) sum(I(Other_conewtssinglepixel)==2);  sum(I(LumIds_conewtssinglepixel)==3) sum(I(ColorOpponentIds_conewtssinglepixel)==3) sum(I(Sconedominated_conewtssinglepixel)==3) sum(I(Other_conewtssinglepixel)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells');  title('single pixel');
subplot(326); bar([sum(I(LumIds_conewtssinglepixel)==1)/numel(LumIds_conewtssinglepixel) sum(I(ColorOpponentIds_conewtssinglepixel)==1)/numel(ColorOpponentIds_conewtssinglepixel) sum(I(Sconedominated_conewtssinglepixel)==1)/numel(Sconedominated_conewtssinglepixel) sum(I(Other_conewtssinglepixel)==1)/numel(Other_conewtssinglepixel); sum(I(LumIds_conewtssinglepixel)==2)/numel(LumIds_conewtssinglepixel) sum(I(ColorOpponentIds_conewtssinglepixel)==2)/numel(ColorOpponentIds_conewtssinglepixel) sum(I(Sconedominated_conewtssinglepixel)==2)/numel(Sconedominated_conewtssinglepixel) sum(I(Other_conewtssinglepixel)==2)/numel(Other_conewtssinglepixel);  sum(I(LumIds_conewtssinglepixel)==3)/numel(LumIds_conewtssinglepixel) sum(I(ColorOpponentIds_conewtssinglepixel)==3)/numel(ColorOpponentIds_conewtssinglepixel) sum(I(Sconedominated_conewtssinglepixel)==3)/numel(Sconedominated_conewtssinglepixel) sum(I(Other_conewtssinglepixel)==3)/numel(Other_conewtssinglepixel)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('single pixel');
plot_counter = plot_counter + 1;

% Plotting the distribution of cells better fit by Gabor, Crescent or DoG as a function  of L+M 
bins = -1:0.1:1;
figure(plot_counter); set(gcf,'Name','Model fits as a function of L+M');
subplot(221); histogram(conewts_svd(1,I1==1)+conewts_svd(2,I1==1),bins,'Normalization','probability'); hold on; histogram(conewts_svd(1,I1==2)+conewts_svd(2,I1==2),bins,'Normalization','probability');
histogram(conewts_svd(1,I1==3)+conewts_svd(2,I1==3),bins,'Normalization','probability'); xlabel('L+M'); ylabel('counts'); title('SVD'); legend({'Cres','Gabor','DoG'},'FontSize',8); hold off;
subplot(222); histogram(conewts_singlepixel(1,I1==1)+conewts_singlepixel(2,I1==1),bins,'Normalization','probability'); hold on; histogram(conewts_singlepixel(1,I1==2)+conewts_singlepixel(2,I1==2),bins,'Normalization','probability');
histogram(conewts_singlepixel(1,I1==3)+conewts_singlepixel(2,I1==3),bins,'Normalization','probability'); xlabel('L+M'); ylabel('counts'); title('single-pixel'); legend({'Cres','Gabor','DoG'},'FontSize',8); hold off;
subplot(223);cdfplot(conewts_svd(1,I1==1)+conewts_svd(2,I1==1)); hold on;cdfplot(conewts_svd(1,I1==2)+conewts_svd(2,I1==2));
cdfplot(conewts_svd(1,I1==3)+conewts_svd(2,I1==3)); legend({'Cres','Gabor','DoG'},'FontSize',8); xlabel('L+M'); title('SVD'); 
subplot(224);cdfplot(conewts_singlepixel(1,I1==1)+conewts_singlepixel(2,I1==1)); hold on;cdfplot(conewts_singlepixel(1,I1==2)+conewts_singlepixel(2,I1==2));
cdfplot(conewts_singlepixel(1,I1==3)+conewts_singlepixel(2,I1==3)); legend({'Cres','Gabor','DoG'},'FontSize',8); xlabel('L+M'); title('single-pixel');
plot_counter = plot_counter + 1;

% Plotting the distribution of cells better fit by Gabor, Crescent or DoG as a function  of S
bins = 0:0.05:1;
figure(plot_counter); set(gcf,'Name','Model fits as a function of S');
subplot(221); histogram(abs(conewts_svd(3,find(I1==1))),bins,'Normalization','probability'); hold on; histogram(abs(conewts_svd(3,find(I1==2))),bins,'Normalization','probability');
histogram(abs(conewts_svd(3,find(I1==3))),bins,'Normalization','probability'); xlabel('S'); ylabel('counts'); title('SVD'); legend({'Cres','Gabor','DoG'},'FontSize',8); hold off;
subplot(222); histogram(abs(conewts_singlepixel(3,find(I1==1))),bins,'Normalization','probability'); hold on; histogram(abs(conewts_singlepixel(3,find(I1==2))),bins,'Normalization','probability');
histogram(abs(conewts_singlepixel(3,find(I1==3))),bins,'Normalization','probability'); xlabel('S'); ylabel('counts'); title('single-pixel'); legend({'Cres','Gabor','DoG'},'FontSize',8); hold off;
subplot(223);cdfplot(abs(conewts_svd(3,find(I1==1)))); hold on;cdfplot(abs(conewts_svd(3,find(I1==2))));
cdfplot(abs(conewts_svd(3,find(I1==3)))); xlabel('S'); title('SVD'); 
subplot(224);cdfplot(abs(conewts_singlepixel(3,find(I1==1)))); hold on;cdfplot(abs(conewts_singlepixel(3,find(I1==2))));
cdfplot(abs(conewts_singlepixel(3,find(I1==3)))); xlabel('S'); title('single-pixel');
plot_counter = plot_counter + 1;

% Plotting the distribution of cells better fit by Gabor, Crescent or DoG as a function  of S+M 
bins = -1:0.1:1;
figure(plot_counter); set(gcf,'Name','Model fits as a function of S+M');
subplot(221); histogram(conewts_svd(3,find(I1==1))+conewts_svd(2,find(I1==1)),bins,'Normalization','probability'); hold on; histogram(conewts_svd(3,find(I1==2))+conewts_svd(2,find(I1==2)),bins,'Normalization','probability');
histogram(conewts_svd(3,find(I1==3))+conewts_svd(2,find(I1==3)),bins,'Normalization','probability'); xlabel('S+M'); ylabel('counts'); title('SVD'); legend({'Cres','Gabor','DoG'},'FontSize',8); hold off;
subplot(222); histogram(conewts_singlepixel(3,find(I1==1))+conewts_singlepixel(2,find(I1==1)),bins,'Normalization','probability'); hold on; histogram(conewts_singlepixel(3,find(I1==2))+conewts_singlepixel(2,find(I1==2)),bins,'Normalization','probability');
histogram(conewts_singlepixel(3,find(I1==3))+conewts_singlepixel(2,find(I1==3)),bins,'Normalization','probability'); xlabel('S+M'); ylabel('counts'); title('single-pixel'); legend({'Cres','Gabor','DoG'},'FontSize',8); hold off;
subplot(223);cdfplot(conewts_svd(3,find(I1==1))+conewts_svd(2,find(I1==1))); hold on;cdfplot(conewts_svd(3,find(I1==2))+conewts_svd(2,find(I1==2)));
cdfplot(conewts_svd(3,find(I1==3))+conewts_svd(2,find(I1==3))); xlabel('S+M'); title('SVD'); 
subplot(224);cdfplot(conewts_singlepixel(3,find(I1==1))+conewts_singlepixel(2,find(I1==1))); hold on;cdfplot(conewts_singlepixel(3,find(I1==2))+conewts_singlepixel(2,find(I1==2)));
cdfplot(conewts_singlepixel(3,find(I1==3))+conewts_singlepixel(2,find(I1==3))); xlabel('S+M'); title('single-pixel');
plot_counter = plot_counter + 1;

% Plotting the distribution of deltaBIC (Crescent & Gabor) for luminance, CO and other cells classified based on 2 criteria: SVD and single-pixel
deltaBIC = CrescentBIC - GaborBIC;
bins = -100:10:300;
Lum_GaborCrescent_svd = LumIds_conewts(I(LumIds_conewts)~=3)'; 
Color_GaborCrescent_svd = ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)~=3)';
Scone_GaborCrescent_svd = Sconedominated_conewts(I(Sconedominated_conewts)~=3)';
Other_GaborCrescent_svd = Other_conewts(I(Other_conewts)~=3)';
Lum_GaborCrescent_singlepixel = LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)~=3)'; 
Color_GaborCrescent_singlepixel = ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)~=3)';
Scone_GaborCrescent_singlepixel = Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)~=3)';
Other_GaborCrescent_singlepixel = Other_conewtssinglepixel(I(Other_conewtssinglepixel)~=3)';
[p1,h1] = ranksum(deltaBIC(Lum_GaborCrescent_svd),deltaBIC(Color_GaborCrescent_svd));
[p2,h2] = ranksum(deltaBIC(Lum_GaborCrescent_singlepixel),deltaBIC(Color_GaborCrescent_singlepixel));
figure(plot_counter); set(gcf,'Name','Crescent vs Gabor');
subplot(421); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Lum_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Lum_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Lum_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]);  xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','lum'); title('Lum: SVD'); hold off;
subplot(422); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Lum_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Lum_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','lum'); title('Lum: single-pixel'); hold off;
subplot(423); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Color_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Color_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]);  xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','CO'); title('CO: SVD'); hold off;
subplot(424); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Color_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Color_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','CO'); title('CO: single-pixel'); hold off;
subplot(425); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Scone_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Scone_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count');  set(gca,'TickDir','out'); legend('All','S'); title('S: SVD'); hold off;
subplot(426); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Scone_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Scone_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','S'); title('S: single-pixel'); hold off;
subplot(427); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Other_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Other_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count');  set(gca,'TickDir','out'); legend('All','Other'); title('Other: SVD'); hold off;
subplot(428); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Other_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Other_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','Other'); title('Other: single-pixel'); hold off;
plot_counter = plot_counter + 1;

% Plotting the distribution of deltaBICs(DoG & Gabor) for luminance, CO and other cells classified based on 2 criteria: SVD and single-pixel
deltaBIC = DOGBIC - GaborBIC;
bins = -100:10:300;
Lum_GaborCrescent_svd = LumIds_conewts(I(LumIds_conewts)~=3)'; 
Color_GaborCrescent_svd = ColorOpponentIds_conewts(I(ColorOpponentIds_conewts)~=3)';
Scone_GaborCrescent_svd = Sconedominated_conewts(I(Sconedominated_conewts)~=3)';
Other_GaborCrescent_svd = Other_conewts(I(Other_conewts)~=3)';
Lum_GaborCrescent_singlepixel = LumIds_conewtssinglepixel(I(LumIds_conewtssinglepixel)~=3)'; 
Color_GaborCrescent_singlepixel = ColorOpponentIds_conewtssinglepixel(I(ColorOpponentIds_conewtssinglepixel)~=3)';
Scone_GaborCrescent_singlepixel = Sconedominated_conewtssinglepixel(I(Sconedominated_conewtssinglepixel)~=3)';
Other_GaborCrescent_singlepixel = Other_conewtssinglepixel(I(Other_conewtssinglepixel)~=3)';
[p1,h1] = ranksum(deltaBIC(Lum_GaborCrescent_svd),deltaBIC(Color_GaborCrescent_svd));
[p2,h2] = ranksum(deltaBIC(Lum_GaborCrescent_singlepixel),deltaBIC(Color_GaborCrescent_singlepixel));
figure(plot_counter); set(gcf,'Name','DOG vs Gabor');
subplot(421); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Lum_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Lum_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Lum_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]);  xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','lum'); title('Lum: SVD'); hold off;
subplot(422); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Lum_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Lum_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','lum'); title('Lum: single-pixel'); hold off;
subplot(423); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Color_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Color_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]);  xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','CO'); title('CO: SVD'); hold off;
subplot(424); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Color_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Color_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','CO'); title('CO: single-pixel'); hold off;
subplot(425); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Scone_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Scone_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count');  set(gca,'TickDir','out'); legend('All','S'); title('S: SVD'); hold off;
subplot(426); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Scone_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Scone_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','S'); title('S: single-pixel'); hold off;
subplot(427); histogram(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd]),bins); hold on; histogram(deltaBIC(Other_GaborCrescent_svd),bins); plot(median(deltaBIC([Lum_GaborCrescent_svd; Color_GaborCrescent_svd; Scone_GaborCrescent_svd; Other_GaborCrescent_svd])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Other_GaborCrescent_svd)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count');  set(gca,'TickDir','out'); legend('All','Other'); title('Other: SVD'); hold off;
subplot(428); histogram(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel]),bins); hold on; histogram(deltaBIC(Other_GaborCrescent_singlepixel),bins); plot(median(deltaBIC([Lum_GaborCrescent_singlepixel; Color_GaborCrescent_singlepixel; Scone_GaborCrescent_singlepixel; Other_GaborCrescent_singlepixel])),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(Other_GaborCrescent_singlepixel)),0,'kv','Markerfacecolor',[1 0 0]); xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out'); legend('All','Other'); title('Other: single-pixel'); hold off;
plot_counter = plot_counter + 1;

%  Scatterplot figure as requested by Greg
figure(plot_counter); set(gcf,'Name','Scatter plot DoG BIC vs Gabor BIC: simple, CO, S-cone dominated'); 
fnames = Output_List(~Singleopponent & simplecells,1);
subplot(121); hold on;
for ii = 1:length(fnames)
    
    if any(ismember(Lum_GaborCrescent_svd,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Color_GaborCrescent_svd,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Scone_GaborCrescent_svd,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    end
end
xlabel('DoG BIC'); ylabel('Gabor BIC'); set(gca,'TickDir','out'); axis square; hold off;
subplot(122); hold on;
for ii = 1:length(fnames)
    if any(ismember(Lum_GaborCrescent_svd,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Color_GaborCrescent_svd,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    end
end
xlabel('DoG BIC'); ylabel('Gabor BIC'); set(gca,'TickDir','out'); axis square; hold off;
plot_counter = plot_counter + 1;

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts); Lumsubplot = ceil(sqrt(numel(Lumind)));
COind = ind(ColorOpponentIds_conewts); COsubplot = ceil(sqrt(numel(COind)));
Sconeind = ind(Sconedominated_conewts); Sconesubplot = ceil(sqrt(numel(Sconeind)));
Limemagentaind = ind(Sconedominated_conewtsLM); Limemagentasubplot = ceil(sqrt(numel(Limemagentaind)));
OrangeCyanind = ind(Sconedominated_conewtsOC); OrangeCyansubplot = ceil(sqrt(numel(OrangeCyanind)));
Otherind = ind(Other_conewts); Othersubplot = ceil(sqrt(numel(Otherind)));
SOind = find(Singleopponent) ; SOsubplot = ceil(sqrt(numel(SOind)));
complexind = find(~simplecells) ; complexsubplot = ceil(sqrt(numel(complexind)));
% plotting STAs of lumninance simple cells
figure(plot_counter); set(gcf,'Name','Luminance cells: SVD');
for ii = 1:numel(Lumind)
    tmp_vec_gun = Output_List{Lumind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(Lumsubplot,Lumsubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of L-M simple cells
figure(plot_counter); set(gcf,'Name','CO: SVD');
for ii = 1:numel(COind)
    tmp_vec_gun = Output_List{COind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(COsubplot,COsubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of S cone-dominated simple cells
figure(plot_counter); set(gcf,'Name','S: SVD');
for ii = 1:numel(Sconeind)
    tmp_vec_gun = Output_List{Sconeind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(Sconesubplot,Sconesubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of S: Lime-Magenta cells
figure(plot_counter); set(gcf,'Name','S: Lime-Magenta');
for ii = 1:numel(Limemagentaind)
    tmp_vec_gun = Output_List{Limemagentaind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(Limemagentasubplot,Limemagentasubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of S: Orange-Cyan cells
figure(plot_counter); set(gcf,'Name','S: Orange-Cyan');
for ii = 1:numel(OrangeCyanind)
    tmp_vec_gun = Output_List{OrangeCyanind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(OrangeCyansubplot,OrangeCyansubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;


% plotting STAs of Other category simple cells
figure(plot_counter); set(gcf,'Name','Other: SVD');
for ii = 1:numel(Otherind)
    tmp_vec_gun = Output_List{Otherind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(Othersubplot,Othersubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of Single-opponent cells
figure(plot_counter); set(gcf,'Name','SO: SVD');
for ii = 1:numel(SOind)
    tmp_vec_gun = Output_List{SOind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(SOsubplot,SOsubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of complex cells
figure(plot_counter); set(gcf,'Name','complex: SVD');
for ii = 1:numel(complexind)
    tmp_vec_gun = Output_List{complexind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(complexsubplot,complexsubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% Seeing if there any relation of color tuning with spatial location
figure(plot_counter); set(gcf,'Name','RF locations of lum, CO and Other cells'); 
subplot(121); plot(RF(LumIds_conewts,1), RF(LumIds_conewts,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RF(ColorOpponentIds_conewts,1), RF(ColorOpponentIds_conewts,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
plot(RF(Sconedominated_conewts,1), RF(Sconedominated_conewts,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
plot(RF(Other_conewts,1), RF(Other_conewts,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square; set(gca,'Xlim',[-100 100],'Ylim',[-100 100]); 
xlabel('X'); ylabel('Y'); title('SVD'); axis square; grid on;
subplot(122); plot(RF(LumIds_conewtssinglepixel,1), RF(LumIds_conewtssinglepixel,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(RF(ColorOpponentIds_conewtssinglepixel,1), RF(ColorOpponentIds_conewtssinglepixel,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
plot(RF(Sconedominated_conewtssinglepixel,1), RF(Sconedominated_conewtssinglepixel,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
plot(RF(Other_conewtssinglepixel,1), RF(Other_conewtssinglepixel,2),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square; set(gca,'Xlim',[-100 100],'Ylim',[-100 100]); 
xlabel('X'); ylabel('Y'); title('single-pixel'); axis square; grid on;
plot_counter = plot_counter + 1;

% Now I want to look at the latency of the luminance, CO and Other cells classified according to 2 criteria: SVD and single-pixel 
bins = 0:10:150;
latencies = cell2mat(Output_List(:,6));
figure(plot_counter); set(gcf,'Name','Latency');
subplot(221); histogram(latencies(ind(LumIds_conewts)),bins,'Normalization','probability'); hold on; histogram(latencies(ind(ColorOpponentIds_conewts)),bins,'Normalization','probability'); histogram(latencies(ind(Sconedominated_conewts)),bins,'Normalization','probability'); histogram(latencies(ind(Other_conewts)),bins,'Normalization','probability'); 
plot(mean(latencies(ind(LumIds_conewts))),0,'kv','Markerfacecolor',[0 0 1]); plot(mean(latencies(ind(ColorOpponentIds_conewts))),0,'kv','Markerfacecolor',[1 0 0]);  plot(mean(latencies(ind(Sconedominated_conewts))),0,'kv','Markerfacecolor',[0 0.5 1]); plot(mean(latencies(ind(Other_conewts))),0,'kv','Markerfacecolor',[1 1 0]);legend({'lum','CO','S','Other'},'FontSize',8); xlabel('latency (ms)'); title('SVD')
subplot(222); histogram(latencies(ind(LumIds_conewtssinglepixel)),bins,'Normalization','probability'); hold on; histogram(latencies(ind(ColorOpponentIds_conewtssinglepixel)),bins,'Normalization','probability'); histogram(latencies(ind(Other_conewtssinglepixel)),bins,'Normalization','probability'); 
plot(mean(latencies(ind(LumIds_conewtssinglepixel))),0,'kv','Markerfacecolor',[0 0 1]); plot(mean(latencies(ind(ColorOpponentIds_conewtssinglepixel))),0,'kv','Markerfacecolor',[1 0 0]); plot(mean(latencies(ind(Sconedominated_conewtssinglepixel))),0,'kv','Markerfacecolor',[0 0.5 1.0]); plot(mean(latencies(ind(Other_conewtssinglepixel))),0,'kv','Markerfacecolor',[1 1 0]); xlabel('latency (ms)'); title('single-pixel')
subplot(223); cdfplot(latencies(ind(LumIds_conewts))); hold on; cdfplot(latencies(ind(ColorOpponentIds_conewts))); cdfplot(latencies(ind(Sconedominated_conewts))); cdfplot(latencies(ind(Other_conewts))); xlabel('latency'); legend({'Lum','CO','S','Other'},'FontSize',8); title('SVD');
subplot(224); cdfplot(latencies(ind(LumIds_conewtssinglepixel))); hold on; cdfplot(latencies(ind(ColorOpponentIds_conewtssinglepixel))); cdfplot(latencies(ind(Sconedominated_conewtssinglepixel))); cdfplot(latencies(ind(Other_conewtssinglepixel))); xlabel('latency'); title('single-pixel');
plot_counter = plot_counter + 1;

% I want to analyze the latencies of the cells classified according to their spatial structure
figure(plot_counter); set(gcf,'Name','Latency vs spatial structure');
subplot(121); histogram(latencies(ind(I==1)),bins,'Normalization','probability'); hold on; histogram(latencies(ind(I==2)),bins,'Normalization','probability'); histogram(latencies(ind(I==3)),bins,'Normalization','probability'); 
plot(mean(latencies(ind(I==1))),0,'kv','Markerfacecolor',[0 0 1]); plot(mean(latencies(ind(I==2))),0,'kv','Markerfacecolor',[1 0 0]); plot(mean(latencies(ind(I==3))),0,'kv','Markerfacecolor',[1 1 0]);legend({'Crescent','Gabor','DoG'},'FontSize',8); xlabel('latency (ms)'); 
subplot(122); cdfplot(latencies(ind(I==1))); hold on; cdfplot(latencies(ind(I==2))); cdfplot(latencies(ind(I==3 )));  xlabel('latency'); legend({'Crescent','Gabor','DoG'},'FontSize',8);
plot_counter = plot_counter + 1;

% I want to analyze latency as function of L+M and S signals 
[r1,p1] = corr(sum(conewts_svd(1:2,:))',latencies(ind),'type','Spearman');
[r2,p2] = corr(sum(conewts_singlepixel(1:2,:))',latencies(ind),'type','Spearman');
[r3,p3] = corr(abs(conewts_svd(3,:))',latencies(ind),'type','Spearman');
[r4,p4] = corr(abs(conewts_singlepixel(3,:))',latencies(ind),'type','Spearman');
figure(plot_counter); set(gcf,'Name','Latency vs cone signals');
subplot(221); plot(sum(conewts_svd(1:2,:)),latencies(ind),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline; xlabel('L+M'); ylabel('latency'); title(strcat('p=',num2str(p1,4)));
subplot(222); plot(sum(conewts_singlepixel(1:2,:)),latencies(ind),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline; xlabel('L+M'); ylabel('latency'); title(strcat('p=',num2str(p2,4)));
subplot(223); plot(abs(conewts_svd(3,:)),latencies(ind),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline; xlabel('S'); ylabel('latency'); title(strcat('p=',num2str(p3,4)));
subplot(224); plot(abs(conewts_singlepixel(3,:)),latencies(ind),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; lsline; xlabel('S'); ylabel('latency'); title(strcat('p=',num2str(p4,4)));
plot_counter = plot_counter + 1;


%% Next I want to classify based on Cross validation performances: percentage accuracy for predicting spike/no-spike and SSE
load Peraccuracy.mat
load SSE.mat
load Deviation.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
meanR = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

[~,I_peracc] = max(meanperaccuracy,[],2); 
figure(plot_counter); set(gcf,'Name','Distribution of dominant model: Per accuracy CV');
I1 = I_peracc([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts Other_conewts]);
subplot(321); bar([numel(find(I1==1)) numel(find(I1==2)) numel(find(I1==3))]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells'); 
subplot(322); bar([numel(find(I1==1))/numel(I1) numel(find(I1==2))/numel(I1) numel(find(I1==3))/numel(I1)]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells');
subplot(323); bar([sum(I_peracc(LumIds_conewts)==1) sum(I_peracc(ColorOpponentIds_conewts)==1) sum(I_peracc(Sconedominated_conewts)==1) sum(I_peracc(Other_conewts)==1); sum(I_peracc(LumIds_conewts)==2) sum(I_peracc(ColorOpponentIds_conewts)==2) sum(I_peracc(Sconedominated_conewts)==2) sum(I_peracc(Other_conewts)==2);  sum(I_peracc(LumIds_conewts)==3) sum(I_peracc(ColorOpponentIds_conewts)==3) sum(I_peracc(Sconedominated_conewts)==3) sum(I_peracc(Other_conewts)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); legend({'Lum','CO','S','Other'},'FontSize',8); ylabel('No. of cells');  title('SVD');
subplot(324); bar([sum(I_peracc(LumIds_conewts)==1)/numel(LumIds_conewts) sum(I_peracc(ColorOpponentIds_conewts)==1)/numel(ColorOpponentIds_conewts) sum(I_peracc(Sconedominated_conewts)==1)/numel(Sconedominated_conewts) sum(I_peracc(Other_conewts)==1)/numel(Other_conewts); sum(I_peracc(LumIds_conewts)==2)/numel(LumIds_conewts) sum(I_peracc(ColorOpponentIds_conewts)==2)/numel(ColorOpponentIds_conewts) sum(I_peracc(Sconedominated_conewts)==2)/numel(Sconedominated_conewts) sum(I_peracc(Other_conewts)==2)/numel(Other_conewts);  sum(I_peracc(LumIds_conewts)==3)/numel(LumIds_conewts) sum(I_peracc(Sconedominated_conewts)==3)/numel(Sconedominated_conewts) sum(I_peracc(ColorOpponentIds_conewts)==3)/numel(ColorOpponentIds_conewts) sum(I_peracc(Other_conewts)==3)/numel(Other_conewts)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('SVD');
subplot(325); bar([sum(I_peracc(LumIds_conewtssinglepixel)==1) sum(I_peracc(ColorOpponentIds_conewtssinglepixel)==1) sum(I_peracc(Sconedominated_conewtssinglepixel)==1) sum(I_peracc(Other_conewtssinglepixel)==1); sum(I_peracc(LumIds_conewtssinglepixel)==2) sum(I_peracc(ColorOpponentIds_conewtssinglepixel)==2) sum(I_peracc(Sconedominated_conewtssinglepixel)==2) sum(I_peracc(Other_conewtssinglepixel)==2);  sum(I_peracc(LumIds_conewtssinglepixel)==3) sum(I_peracc(ColorOpponentIds_conewtssinglepixel)==3) sum(I_peracc(Sconedominated_conewtssinglepixel)==3) sum(I_peracc(Other_conewtssinglepixel)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells');  title('single pixel');
subplot(326); bar([sum(I_peracc(LumIds_conewtssinglepixel)==1)/numel(LumIds_conewtssinglepixel) sum(I_peracc(ColorOpponentIds_conewtssinglepixel)==1)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_peracc(Sconedominated_conewtssinglepixel)==1)/numel(Sconedominated_conewtssinglepixel) sum(I_peracc(Other_conewtssinglepixel)==1)/numel(Other_conewtssinglepixel); sum(I_peracc(LumIds_conewtssinglepixel)==2)/numel(LumIds_conewtssinglepixel) sum(I_peracc(ColorOpponentIds_conewtssinglepixel)==2)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_peracc(Sconedominated_conewtssinglepixel)==2)/numel(Sconedominated_conewtssinglepixel) sum(I_peracc(Other_conewtssinglepixel)==2)/numel(Other_conewtssinglepixel);  sum(I_peracc(LumIds_conewtssinglepixel)==3)/numel(LumIds_conewtssinglepixel) sum(I_peracc(ColorOpponentIds_conewtssinglepixel)==3)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_peracc(Sconedominated_conewtssinglepixel)==3)/numel(Sconedominated_conewtssinglepixel) sum(I_peracc(Other_conewtssinglepixel)==3)/numel(Other_conewtssinglepixel)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('single pixel');
plot_counter = plot_counter + 1;

[~,I_SSE] = min(meanSSE,[],2); 
figure(plot_counter); set(gcf,'Name','Distribution of dominant model: mean SSE CV');
I1 = I_SSE([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts Other_conewts]);
subplot(321); bar([numel(find(I1==1)) numel(find(I1==2)) numel(find(I1==3))]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells'); 
subplot(322); bar([numel(find(I1==1))/numel(I1) numel(find(I1==2))/numel(I1) numel(find(I1==3))/numel(I1)]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells');
subplot(323); bar([sum(I_SSE(LumIds_conewts)==1) sum(I_SSE(ColorOpponentIds_conewts)==1) sum(I_SSE(Sconedominated_conewts)==1) sum(I_SSE(Other_conewts)==1); sum(I_SSE(LumIds_conewts)==2) sum(I_SSE(ColorOpponentIds_conewts)==2) sum(I_SSE(Sconedominated_conewts)==2) sum(I_SSE(Other_conewts)==2);  sum(I_SSE(LumIds_conewts)==3) sum(I_SSE(ColorOpponentIds_conewts)==3) sum(I_SSE(Sconedominated_conewts)==3) sum(I_SSE(Other_conewts)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); legend({'Lum','CO','S','Other'},'FontSize',8); ylabel('No. of cells');  title('SVD');
subplot(324); bar([sum(I_SSE(LumIds_conewts)==1)/numel(LumIds_conewts) sum(I_SSE(ColorOpponentIds_conewts)==1)/numel(ColorOpponentIds_conewts) sum(I_SSE(Sconedominated_conewts)==1)/numel(Sconedominated_conewts) sum(I_SSE(Other_conewts)==1)/numel(Other_conewts); sum(I_SSE(LumIds_conewts)==2)/numel(LumIds_conewts) sum(I_SSE(ColorOpponentIds_conewts)==2)/numel(ColorOpponentIds_conewts) sum(I_SSE(Sconedominated_conewts)==2)/numel(Sconedominated_conewts) sum(I_SSE(Other_conewts)==2)/numel(Other_conewts);  sum(I_SSE(LumIds_conewts)==3)/numel(LumIds_conewts) sum(I_SSE(Sconedominated_conewts)==3)/numel(Sconedominated_conewts) sum(I_SSE(ColorOpponentIds_conewts)==3)/numel(ColorOpponentIds_conewts) sum(I_SSE(Other_conewts)==3)/numel(Other_conewts)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('SVD');
subplot(325); bar([sum(I_SSE(LumIds_conewtssinglepixel)==1) sum(I_SSE(ColorOpponentIds_conewtssinglepixel)==1) sum(I_SSE(Sconedominated_conewtssinglepixel)==1) sum(I_SSE(Other_conewtssinglepixel)==1); sum(I_SSE(LumIds_conewtssinglepixel)==2) sum(I_SSE(ColorOpponentIds_conewtssinglepixel)==2) sum(I_SSE(Sconedominated_conewtssinglepixel)==2) sum(I_SSE(Other_conewtssinglepixel)==2);  sum(I_SSE(LumIds_conewtssinglepixel)==3) sum(I(ColorOpponentIds_conewtssinglepixel)==3) sum(I_SSE(Sconedominated_conewtssinglepixel)==3) sum(I_SSE(Other_conewtssinglepixel)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells');  title('single pixel');
subplot(326); bar([sum(I_SSE(LumIds_conewtssinglepixel)==1)/numel(LumIds_conewtssinglepixel) sum(I_SSE(ColorOpponentIds_conewtssinglepixel)==1)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_SSE(Sconedominated_conewtssinglepixel)==1)/numel(Sconedominated_conewtssinglepixel) sum(I_SSE(Other_conewtssinglepixel)==1)/numel(Other_conewtssinglepixel); sum(I_SSE(LumIds_conewtssinglepixel)==2)/numel(LumIds_conewtssinglepixel) sum(I_SSE(ColorOpponentIds_conewtssinglepixel)==2)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_SSE(Sconedominated_conewtssinglepixel)==2)/numel(Sconedominated_conewtssinglepixel) sum(I_SSE(Other_conewtssinglepixel)==2)/numel(Other_conewtssinglepixel);  sum(I_SSE(LumIds_conewtssinglepixel)==3)/numel(LumIds_conewtssinglepixel) sum(I_SSE(ColorOpponentIds_conewtssinglepixel)==3)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_SSE(Sconedominated_conewtssinglepixel)==3)/numel(Sconedominated_conewtssinglepixel) sum(I_SSE(Other_conewtssinglepixel)==3)/numel(Other_conewtssinglepixel)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('single pixel');
plot_counter = plot_counter + 1;

[~,I_R] = max(meanR,[],2); 
figure(plot_counter); set(gcf,'Name','Distribution of dominant model: mean R CV');
I1 = I_R([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts Other_conewts]);
subplot(321); bar([numel(find(I1==1)) numel(find(I1==2)) numel(find(I1==3))]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells'); 
subplot(322); bar([numel(find(I1==1))/numel(I1) numel(find(I1==2))/numel(I1) numel(find(I1==3))/numel(I1)]); set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells');
subplot(323); bar([sum(I_R(LumIds_conewts)==1) sum(I_R(ColorOpponentIds_conewts)==1) sum(I_R(Sconedominated_conewts)==1) sum(I_R(Other_conewts)==1); sum(I_R(LumIds_conewts)==2) sum(I_R(ColorOpponentIds_conewts)==2) sum(I_R(Sconedominated_conewts)==2) sum(I_R(Other_conewts)==2);  sum(I_R(LumIds_conewts)==3) sum(I_R(ColorOpponentIds_conewts)==3) sum(I_R(Sconedominated_conewts)==3) sum(I_R(Other_conewts)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); legend({'Lum','CO','S','Other'},'FontSize',8); ylabel('No. of cells');  title('SVD');
subplot(324); bar([sum(I_R(LumIds_conewts)==1)/numel(LumIds_conewts) sum(I_R(ColorOpponentIds_conewts)==1)/numel(ColorOpponentIds_conewts) sum(I_R(Sconedominated_conewts)==1)/numel(Sconedominated_conewts) sum(I_R(Other_conewts)==1)/numel(Other_conewts); sum(I_R(LumIds_conewts)==2)/numel(LumIds_conewts) sum(I_R(ColorOpponentIds_conewts)==2)/numel(ColorOpponentIds_conewts) sum(I_R(Sconedominated_conewts)==2)/numel(Sconedominated_conewts) sum(I_R(Other_conewts)==2)/numel(Other_conewts);  sum(I_R(LumIds_conewts)==3)/numel(LumIds_conewts) sum(I_R(Sconedominated_conewts)==3)/numel(Sconedominated_conewts) sum(I_R(ColorOpponentIds_conewts)==3)/numel(ColorOpponentIds_conewts) sum(I_R(Other_conewts)==3)/numel(Other_conewts)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('SVD');
subplot(325); bar([sum(I_R(LumIds_conewtssinglepixel)==1) sum(I_R(ColorOpponentIds_conewtssinglepixel)==1) sum(I_R(Sconedominated_conewtssinglepixel)==1) sum(I_R(Other_conewtssinglepixel)==1); sum(I_R(LumIds_conewtssinglepixel)==2) sum(I_R(ColorOpponentIds_conewtssinglepixel)==2) sum(I_R(Sconedominated_conewtssinglepixel)==2) sum(I_R(Other_conewtssinglepixel)==2);  sum(I_R(LumIds_conewtssinglepixel)==3) sum(I(ColorOpponentIds_conewtssinglepixel)==3) sum(I_R(Sconedominated_conewtssinglepixel)==3) sum(I_R(Other_conewtssinglepixel)==3)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('No. of cells');  title('single pixel');
subplot(326); bar([sum(I_R(LumIds_conewtssinglepixel)==1)/numel(LumIds_conewtssinglepixel) sum(I_R(ColorOpponentIds_conewtssinglepixel)==1)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_R(Sconedominated_conewtssinglepixel)==1)/numel(Sconedominated_conewtssinglepixel) sum(I_R(Other_conewtssinglepixel)==1)/numel(Other_conewtssinglepixel); sum(I_R(LumIds_conewtssinglepixel)==2)/numel(LumIds_conewtssinglepixel) sum(I_R(ColorOpponentIds_conewtssinglepixel)==2)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_R(Sconedominated_conewtssinglepixel)==2)/numel(Sconedominated_conewtssinglepixel) sum(I_R(Other_conewtssinglepixel)==2)/numel(Other_conewtssinglepixel);  sum(I_R(LumIds_conewtssinglepixel)==3)/numel(LumIds_conewtssinglepixel) sum(I_R(ColorOpponentIds_conewtssinglepixel)==3)/numel(ColorOpponentIds_conewtssinglepixel) sum(I_R(Sconedominated_conewtssinglepixel)==3)/numel(Sconedominated_conewtssinglepixel) sum(I_R(Other_conewtssinglepixel)==3)/numel(Other_conewtssinglepixel)]);  
set(gca,'xticklabel',{'Crescent','Gabor','DOG'}); ylabel('Proportion of cells'); title('single pixel');
plot_counter = plot_counter + 1;

figure(plot_counter); set(gcf,'Name','Cross validation & BICs: Lum');
subplot(341),plot(meanperaccuracy(LumIds_conewts,2),meanperaccuracy(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(342),plot(meanSSE(LumIds_conewts,2),meanSSE(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 3],[0 3]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(343),plot(BICs(LumIds_conewts,2),BICs(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(344),plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
subplot(345),plot(meanperaccuracy(LumIds_conewts,1),meanperaccuracy(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(346),plot(meanSSE(LumIds_conewts,1),meanSSE(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 3],[0 3]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(347),plot(BICs(LumIds_conewts,1),BICs(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(348),plot(meanR(LumIds_conewts,1),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
subplot(349),plot(meanperaccuracy(LumIds_conewts,2),meanperaccuracy(LumIds_conewts,1),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(3,4,10),plot(meanSSE(LumIds_conewts,2),meanSSE(LumIds_conewts,1),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 3],[0 3]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:30,'YTick',0:1.5:30); title('CV SSE');
subplot(3,4,11),plot(BICs(LumIds_conewts,2),BICs(LumIds_conewts,1),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(3,4,12),plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');

set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

figure(plot_counter); set(gcf,'Name','Cross validation & BICs: L-M');
subplot(341),plot(meanperaccuracy(ColorOpponentIds_conewts,2),meanperaccuracy(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(342), plot(meanSSE(ColorOpponentIds_conewts,2),meanSSE(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 3],[0 3]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(343),plot(BICs(ColorOpponentIds_conewts,2),BICs(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(344),plot(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
subplot(345), plot(meanperaccuracy(ColorOpponentIds_conewts,1),meanperaccuracy(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(346), plot(meanSSE(ColorOpponentIds_conewts,1),meanSSE(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 3],[0 3]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(347), plot(BICs(ColorOpponentIds_conewts,1),BICs(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on
line([-1200 -200],[-1200 -200]); xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(348),plot(meanR(ColorOpponentIds_conewts,1),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]); xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
subplot(349), plot(meanperaccuracy(ColorOpponentIds_conewts,2),meanperaccuracy(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(3,4,10), plot(meanSSE(ColorOpponentIds_conewts,2),meanSSE(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 3],[0 3]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(3,4,11), plot(BICs(ColorOpponentIds_conewts,2),BICs(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;  
line([-1200 -200],[-1200 -200]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(3,4,12),plot(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

figure(plot_counter); set(gcf,'Name','Cross validation & BICs: S');
subplot(341),plot(meanperaccuracy(Sconedominated_conewts,2),meanperaccuracy(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(342), plot(meanSSE(Sconedominated_conewts,2),meanSSE(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 3],[0 3]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(343), plot(BICs(Sconedominated_conewts,2),BICs(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(344),plot(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
subplot(345), plot(meanperaccuracy(Sconedominated_conewts,1),meanperaccuracy(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(346), plot(meanSSE(Sconedominated_conewts,1),meanSSE(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 3],[0 3]);  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(347), plot(BICs(Sconedominated_conewts,1),BICs(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(348),plot(meanR(Sconedominated_conewts,1),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]); xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
subplot(349), plot(meanperaccuracy(Sconedominated_conewts,2),meanperaccuracy(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(3,4,10), plot(meanSSE(Sconedominated_conewts,2),meanSSE(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 3],[0 3]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 3],'Ylim',[0 3],'XTick',0:1.5:3.0,'YTick',0:1.5:3.0); title('CV SSE');
subplot(3,4,11), plot(BICs(Sconedominated_conewts,2),BICs(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1200 -200],[-1200 -200]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:500:-200,'YTick',-1200:500:-200); title('BIC');
subplot(3,4,12),plot(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1]); xlabel('Gabor'); ylabel('Crescent');axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.25:1.0,'YTick',0:0.25:1.0); title('CV R');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;


% Going into the histogram land and comparing Gabor and DOG 
All_ind = [LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts];
meanperaccuracy_GaborDOG = meanperaccuracy(:,2)-meanperaccuracy(:,3);
meanSSE_GaborDOG = meanSSE(:,3)-meanSSE(:,2);
BIC_GaborDOG = DOGBIC-GaborBIC;
bins1 = -0.05:0.01:0.1;
bins2 = -0.5:0.05:0.5;
bins3 = -100:20:300;
figure(plot_counter); set(gcf,'Name','DOG vs Gabor');
subplot(331), histogram(meanperaccuracy_GaborDOG(All_ind),bins1); hold on; plot(median(meanperaccuracy_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanperaccuracy_GaborDOG(LumIds_conewts),bins1); plot(median(meanperaccuracy_GaborDOG(LumIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('Lum: CV per accuracy'); set(gca,'Tickdir','out'); 
subplot(332), histogram(meanSSE_GaborDOG(All_ind),bins2); hold on; plot(median(meanSSE_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanSSE_GaborDOG(LumIds_conewts),bins2); plot(median(meanSSE_GaborDOG(LumIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('Lum: CV SSE'); set(gca,'Tickdir','out'); 
subplot(333), histogram(BIC_GaborDOG(All_ind),bins3); hold on; plot(median(BIC_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(BIC_GaborDOG(LumIds_conewts),bins3); plot(median(BIC_GaborDOG(LumIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('Lum: BIC'); set(gca,'Tickdir','out'); 
subplot(334), histogram(meanperaccuracy_GaborDOG(All_ind),bins1); hold on; plot(median(meanperaccuracy_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanperaccuracy_GaborDOG(ColorOpponentIds_conewts),bins1); plot(median(meanperaccuracy_GaborDOG(ColorOpponentIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('L-M: CV per accuracy'); set(gca,'Tickdir','out'); 
subplot(335), histogram(meanSSE_GaborDOG(All_ind),bins2); hold on; plot(median(meanSSE_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanSSE_GaborDOG(ColorOpponentIds_conewts),bins2); plot(median(meanSSE_GaborDOG(ColorOpponentIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('L-M: CV SSE'); set(gca,'Tickdir','out'); 
subplot(336), histogram(BIC_GaborDOG(All_ind),bins3); hold on; plot(median(BIC_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(BIC_GaborDOG(ColorOpponentIds_conewts),bins3); plot(median(BIC_GaborDOG(ColorOpponentIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('L-M: BIC'); set(gca,'Tickdir','out'); 
subplot(337), histogram(meanperaccuracy_GaborDOG(All_ind),bins1); hold on; plot(median(meanperaccuracy_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanperaccuracy_GaborDOG(Sconedominated_conewts),bins1); plot(median(meanperaccuracy_GaborDOG(Sconedominated_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('S: CV per accuracy'); set(gca,'Tickdir','out'); 
subplot(338), histogram(meanSSE_GaborDOG(All_ind),bins2); hold on; plot(median(meanSSE_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanSSE_GaborDOG(Sconedominated_conewts),bins2); plot(median(meanSSE_GaborDOG(Sconedominated_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('S: CV SSE'); set(gca,'Tickdir','out'); 
subplot(339), histogram(BIC_GaborDOG(All_ind),bins3); hold on; plot(median(BIC_GaborDOG(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(BIC_GaborDOG(Sconedominated_conewts),bins3); plot(median(BIC_GaborDOG(Sconedominated_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - DOG'); title('S: BIC'); set(gca,'Tickdir','out'); 
plot_counter = plot_counter + 1;

% Gabor vs Crescent  
meanperaccuracy_GaborCrescent = meanperaccuracy(:,2)-meanperaccuracy(:,1);
meanSSE_GaborCrescent = meanSSE(:,1)-meanSSE(:,2);
BIC_GaborCrescent = CrescentBIC-GaborBIC;
bins1 = -0.05:0.01:0.1;
bins2 = -0.5:0.05:0.5;
bins3 = -100:20:300;
figure(plot_counter); set(gcf,'Name','Crescent vs Gabor');
subplot(331), histogram(meanperaccuracy_GaborCrescent(All_ind),bins1); hold on; plot(median(meanperaccuracy_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanperaccuracy_GaborCrescent(LumIds_conewts),bins1); plot(median(meanperaccuracy_GaborCrescent(LumIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('Lum: CV per accuracy'); set(gca,'Tickdir','out'); 
subplot(332), histogram(meanSSE_GaborCrescent(All_ind),bins2); hold on; plot(median(meanSSE_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanSSE_GaborCrescent(LumIds_conewts),bins2); plot(median(meanSSE_GaborCrescent(LumIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('Lum: CV SSE'); set(gca,'Tickdir','out'); 
subplot(333), histogram(BIC_GaborCrescent(All_ind),bins3); hold on; plot(median(BIC_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(BIC_GaborCrescent(LumIds_conewts),bins3); plot(median(BIC_GaborCrescent(LumIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('Lum: BIC'); set(gca,'Tickdir','out'); 
subplot(334), histogram(meanperaccuracy_GaborCrescent(All_ind),bins1); hold on; plot(median(meanperaccuracy_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanperaccuracy_GaborCrescent(ColorOpponentIds_conewts),bins1); plot(median(meanperaccuracy_GaborCrescent(ColorOpponentIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('L-M: CV per accuracy'); set(gca,'Tickdir','out'); 
subplot(335), histogram(meanSSE_GaborCrescent(All_ind),bins2); hold on; plot(median(meanSSE_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanSSE_GaborCrescent(ColorOpponentIds_conewts),bins2); plot(median(meanSSE_GaborCrescent(ColorOpponentIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('L-M: CV SSE'); set(gca,'Tickdir','out'); 
subplot(336), histogram(BIC_GaborCrescent(All_ind),bins3); hold on; plot(median(BIC_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(BIC_GaborCrescent(ColorOpponentIds_conewts),bins3); plot(median(BIC_GaborCrescent(ColorOpponentIds_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('L-M: BIC'); set(gca,'Tickdir','out'); 
subplot(337), histogram(meanperaccuracy_GaborCrescent(All_ind),bins1); hold on; plot(median(meanperaccuracy_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanperaccuracy_GaborCrescent(Sconedominated_conewts),bins1); plot(median(meanperaccuracy_GaborCrescent(Sconedominated_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('S: CV per accuracy'); set(gca,'Tickdir','out'); 
subplot(338), histogram(meanSSE_GaborCrescent(All_ind),bins2); hold on; plot(median(meanSSE_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(meanSSE_GaborCrescent(Sconedominated_conewts),bins2); plot(median(meanSSE_GaborCrescent(Sconedominated_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('S: CV SSE'); set(gca,'Tickdir','out'); 
subplot(339), histogram(BIC_GaborCrescent(All_ind),bins3); hold on; plot(median(BIC_GaborCrescent(All_ind)),0,'kv','MarkerFaceColor',[0 0 1]); histogram(BIC_GaborCrescent(Sconedominated_conewts),bins3); plot(median(BIC_GaborCrescent(Sconedominated_conewts)),0,'kv','MarkerFaceColor',[1 0 0]); xlabel('Gabor - Crescent'); title('S: BIC'); set(gca,'Tickdir','out'); 
plot_counter = plot_counter + 1;

% Next, I am going to plot the phases of Gabor fits for Lum, L-M and S cone dominates cells
bins = 0:10:90;
Lumid = zeros(size(gaborphases_relevantcells))'; Lumid(LumIds_conewts) = 1;
COid = zeros(size(gaborphases_relevantcells))'; COid(ColorOpponentIds_conewts) = 1;
Sid = zeros(size(gaborphases_relevantcells))'; Sid(Sconedominated_conewts) = 1;
LMid = zeros(size(gaborphases_relevantcells))'; LMid(Sconedominated_conewtsLM) = 1;
OCid = zeros(size(gaborphases_relevantcells))'; OCid(Sconedominated_conewtsOC) = 1;
figure(plot_counter); set(gcf,'Name','Gabor Phases');
subplot(531),histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanperaccuracy_GaborDOG>0)),bins,'FaceColor',[0 1 0]); set(gca,'Tickdir','out'); title('Lum: CV per accuracy'); 
subplot(532),histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanSSE_GaborDOG>0)),bins,'FaceColor',[0 1 0]); set(gca,'Tickdir','out'); title('Lum: CV SSE'); 
subplot(533),histogram(90-abs(90-gaborphases_relevantcells(Lumid & BIC_GaborDOG>0)),bins,'FaceColor',[0 1 0]); set(gca,'Tickdir','out'); title('Lum: BIC'); 
subplot(534),histogram(90-abs(90-gaborphases_relevantcells(COid & meanperaccuracy_GaborDOG>0)),bins,'FaceColor',[1 0.5 0]); set(gca,'Tickdir','out'); title('L-M: CV per accuracy'); 
subplot(535),histogram(90-abs(90-gaborphases_relevantcells(COid & meanSSE_GaborDOG>0)),bins,'FaceColor',[1 0.5 0]); set(gca,'Tickdir','out'); title('L-M: CV SSE'); 
subplot(536),histogram(90-abs(90-gaborphases_relevantcells(COid & BIC_GaborDOG>0)),bins,'FaceColor',[1 0.5 0]); set(gca,'Tickdir','out'); title('L-M: BIC');
subplot(537),histogram(90-abs(90-gaborphases_relevantcells(Sid & meanperaccuracy_GaborDOG>0)),bins,'FaceColor',[0 0.5 1]); set(gca,'Tickdir','out'); title('S: CV per accuracy'); 
subplot(538),histogram(90-abs(90-gaborphases_relevantcells(Sid & meanSSE_GaborDOG>0)),bins,'FaceColor',[0 0.5 1]); set(gca,'Tickdir','out'); title('S: CV SSE');
subplot(539),histogram(90-abs(90-gaborphases_relevantcells(Sid & BIC_GaborDOG>0)),bins,'FaceColor',[0 0.5 1]); set(gca,'Tickdir','out'); title('S: BIC');
subplot(5,3,10),histogram(90-abs(90-gaborphases_relevantcells(LMid & meanperaccuracy_GaborDOG>0)),bins,'FaceColor','m'); set(gca,'Tickdir','out'); title('LM: CV per accuracy'); 
subplot(5,3,11),histogram(90-abs(90-gaborphases_relevantcells(LMid & meanSSE_GaborDOG>0)),bins,'FaceColor','m'); set(gca,'Tickdir','out'); title('LM: CV SSE');
subplot(5,3,12),histogram(90-abs(90-gaborphases_relevantcells(LMid & BIC_GaborDOG>0)),bins,'FaceColor','m'); set(gca,'Tickdir','out'); title('LM: BIC');
subplot(5,3,13),histogram(90-abs(90-gaborphases_relevantcells(OCid & meanperaccuracy_GaborDOG>0)),bins,'FaceColor',[0 0.75 0.75]); set(gca,'Tickdir','out'); title('OC: CV per accuracy'); 
subplot(5,3,14),histogram(90-abs(90-gaborphases_relevantcells(OCid & meanSSE_GaborDOG>0)),bins,'FaceColor',[0 0.75 0.75]); set(gca,'Tickdir','out'); title('OC: CV SSE');
subplot(5,3,15),histogram(90-abs(90-gaborphases_relevantcells(OCid & BIC_GaborDOG>0)),bins,'FaceColor',[0 0.75 0.75]); set(gca,'Tickdir','out'); title('OC: BIC');
plot_counter = plot_counter + 1;

% Next, I am going to plot the aspect-ratios of Gabor fits for Lum, L-M and S cone dominates cells
bins = 0:0.5:5;
figure(plot_counter); set(gcf,'Name','Aspect ratios');
subplot(531),histogram(aspectratio_relevantcells(Lumid & meanperaccuracy_GaborDOG>0),bins,'FaceColor',[0 1 0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanperaccuracy_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('Lum: CV per accuracy'); 
subplot(532),histogram(aspectratio_relevantcells(Lumid & meanSSE_GaborDOG>0),bins,'FaceColor',[0 1 0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanSSE_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('Lum: CV SSE'); 
subplot(533),histogram(aspectratio_relevantcells(Lumid & BIC_GaborDOG>0),bins,'FaceColor',[0 1 0]);  hold on; plot(median(aspectratio_relevantcells(Lumid & BIC_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('Lum: BIC'); 
subplot(534),histogram(aspectratio_relevantcells(COid & meanperaccuracy_GaborDOG>0),bins,'FaceColor',[1 0.5 0]); hold on; plot(median(aspectratio_relevantcells(COid & meanperaccuracy_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('L-M: CV per accuracy'); 
subplot(535),histogram(aspectratio_relevantcells(COid & meanSSE_GaborDOG>0),bins,'FaceColor',[1 0.5 0]); hold on; plot(median(aspectratio_relevantcells(COid & meanSSE_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('L-M: CV SSE'); 
subplot(536),histogram(aspectratio_relevantcells(COid & BIC_GaborDOG>0),bins,'FaceColor',[1 0.5 0]); hold on; plot(median(aspectratio_relevantcells(COid & BIC_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('L-M: BIC');
subplot(537),histogram(aspectratio_relevantcells(Sid & meanperaccuracy_GaborDOG>0),bins,'FaceColor',[0 0.5 1]); hold on; plot(median(aspectratio_relevantcells(Sid & meanperaccuracy_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('S: CV per accuracy'); 
subplot(538),histogram(aspectratio_relevantcells(Sid & meanSSE_GaborDOG>0),bins,'FaceColor',[0 0.5 1]);  hold on; plot(median(aspectratio_relevantcells(Sid & meanSSE_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('S: CV SSE');
subplot(539),histogram(aspectratio_relevantcells(Sid & BIC_GaborDOG>0),bins,'FaceColor',[0 0.5 1]); hold on; plot(median(aspectratio_relevantcells(Sid & BIC_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('S: BIC');
subplot(5,3,10),histogram(aspectratio_relevantcells(LMid & meanperaccuracy_GaborDOG>0),bins,'FaceColor','m'); hold on; plot(median(aspectratio_relevantcells(LMid & meanperaccuracy_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('LM: CV per accuracy'); 
subplot(5,3,11),histogram(aspectratio_relevantcells(LMid & meanSSE_GaborDOG>0),bins,'FaceColor','m');  hold on; plot(median(aspectratio_relevantcells(LMid & meanSSE_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('LM: CV SSE');
subplot(5,3,12),histogram(aspectratio_relevantcells(LMid & BIC_GaborDOG>0),bins,'FaceColor','m'); hold on; plot(median(aspectratio_relevantcells(LMid & BIC_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('LM: BIC');
subplot(5,3,13),histogram(aspectratio_relevantcells(OCid & meanperaccuracy_GaborDOG>0),bins,'FaceColor',[0 0.75 0.75]); hold on; plot(median(aspectratio_relevantcells(OCid & meanperaccuracy_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('OC: CV per accuracy'); 
subplot(5,3,14),histogram(aspectratio_relevantcells(OCid & meanSSE_GaborDOG>0),bins,'FaceColor',[0 0.75 0.75]);  hold on; plot(median(aspectratio_relevantcells(OCid & meanSSE_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('OC: CV SSE');
subplot(5,3,15),histogram(aspectratio_relevantcells(OCid & BIC_GaborDOG>0),bins,'FaceColor',[0 0.75 0.75]); hold on; plot(median(aspectratio_relevantcells(OCid & BIC_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out'); title('OC: BIC');
plot_counter = plot_counter + 1;

%% Here I am plotting the best Gabor, Crescent or DOGfit for luminance, CO and Other cells based on SVD criteria 
LumCrescentind = find(I(Lumind)==1); COCrescentind = find(I(COind)==1); OtherCrescentind = find(I(Otherind)==1); 
LumGaborind = find(I(Lumind)==2); COGaborind = find(I(COind)==2); OtherGaborind = find(I(Otherind)==2); 
LumDOGind = find(I(Lumind)==3); CODOGind = find(I(COind)==3); OtherDOGind = find(I(Otherind)==3);
[~,idx1] = max(medianBICval(Lumind(LumCrescentind))); [~,idx2] = max(medianBICval(COind(COCrescentind))); [~,idx3] = max(medianBICval(Otherind(OtherCrescentind)));
[~,idx4] = max(medianBICval(Lumind(LumGaborind))); [~,idx5] = max(medianBICval(COind(COGaborind))); [~,idx6] = max(medianBICval(Otherind(OtherGaborind)));
[~,idx7] = max(medianBICval(Lumind(LumDOGind))); [~,idx8] = max(medianBICval(COind(CODOGind))); [~,idx9] = max(medianBICval(Otherind(OtherDOGind)));
idxs = [Lumind(LumCrescentind(idx1)); COind(COCrescentind(idx2)); Otherind(OtherCrescentind(idx3));...
    Lumind(LumGaborind(idx4)); COind(COGaborind(idx5)); Otherind(OtherGaborind(idx6));...
    Lumind(LumDOGind(idx7)); COind(CODOGind(idx8)); Otherind(OtherDOGind(idx9))];
figure(plot_counter); set(gcf,'Name','Best examples from each category: SVD');
figure(plot_counter+1); set(gcf,'Name','SVD based RF');
for ii = 1:numel(idxs)
    tmp_vec_gun = Output_List{idxs(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
    figure(plot_counter); subplot(3,3,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    if ii==1 
        title('Lum'); ylabel('Crescent');
    elseif ii == 2
        title('CO');
    elseif ii == 3
        title('Other');
    elseif ii == 4
        ylabel('Gabor');
    elseif ii == 7
        ylabel('DOG');
    end
    figure(plot_counter+1); subplot(3,3,ii); imagesc(Output_List{idxs(ii),4}); set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray');
    
end
plot_counter = plot_counter + 2;


% Doing the same thing as before but for single-pixel based classification
% not plotting these as of now.
plotcellssinglepixels = 0;
if plotcellssinglepixels
    Lumind = ind(LumIds_conewtssinglepixel); Lumsubplot = ceil(sqrt(numel(Lumind)));
    COind = ind(ColorOpponentIds_conewtssinglepixel); COsubplot = ceil(sqrt(numel(COind)));
    Otherind = ind(Other_conewtssinglepixel); Othersubplot = ceil(sqrt(numel(Otherind)));
    figure(plot_counter); set(gcf,'Name','Luminance cells: single-pixel');
    for ii = 1:numel(Lumind)
        tmp_vec_gun = Output_List{Lumind(ii),2};
        normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
        im = normfactor*tmp_vec_gun + 0.5;
        im = reshape(im,[10 10 3]);
        
        subplot(Lumsubplot,Lumsubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    plot_counter = plot_counter + 1;
    
    figure(plot_counter); set(gcf,'Name','Color cells: single-pixel');
    for ii = 1:numel(COind)
        tmp_vec_gun = Output_List{COind(ii),2};
        normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
        im = normfactor*tmp_vec_gun + 0.5;
        im = reshape(im,[10 10 3]);
        
        subplot(COsubplot,COsubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    plot_counter = plot_counter + 1;
    
    figure(plot_counter); set(gcf,'Name','Other cells: single-pixel');
    for ii = 1:numel(Otherind)
        tmp_vec_gun = Output_List{Otherind(ii),2};
        normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
        im = normfactor*tmp_vec_gun + 0.5;
        im = reshape(im,[10 10 3]);
        
        subplot(Othersubplot,Othersubplot,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    plot_counter = plot_counter + 1;
    
    LumCrescentind = find(I(Lumind)==1); COCrescentind = find(I(COind)==1); OtherCrescentind = find(I(Otherind)==1);
    LumGaborind = find(I(Lumind)==2); COGaborind = find(I(COind)==2); OtherGaborind = find(I(Otherind)==2);
    LumDOGind = find(I(Lumind)==3); CODOGind = find(I(COind)==3); OtherDOGind = find(I(Otherind)==3);
    [~,idx1] = max(medianBICval(Lumind(LumCrescentind))); [~,idx2] = max(medianBICval(COind(COCrescentind))); [~,idx3] = max(medianBICval(Otherind(OtherCrescentind)));
    [~,idx4] = max(medianBICval(Lumind(LumGaborind))); [~,idx5] = max(medianBICval(COind(COGaborind))); [~,idx6] = max(medianBICval(Otherind(OtherGaborind)));
    [~,idx7] = max(medianBICval(Lumind(LumDOGind))); [~,idx8] = max(medianBICval(COind(CODOGind))); [~,idx9] = max(medianBICval(Otherind(OtherDOGind)));
    idxs = [Lumind(LumCrescentind(idx1)); COind(COCrescentind(idx2)); Otherind(OtherCrescentind(idx3));...
        Lumind(LumGaborind(idx4)); COind(COGaborind(idx5)); Otherind(OtherGaborind(idx6));...
        Lumind(LumDOGind(idx7)); COind(CODOGind(idx8)); Otherind(OtherDOGind(idx9))];
    figure(plot_counter); set(gcf,'Name','Best examples from each category');
    figure(plot_counter+1); set(gcf,'Name','SVD based RF');
    for ii = 1:numel(idxs)
        tmp_vec_gun = Output_List{idxs(ii),2};
        normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
        im = normfactor*tmp_vec_gun + 0.5;
        im = reshape(im,[10 10 3]);
        subplot(3,3,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
        if ii==1
            title('Lum'); ylabel('Crescent');
        elseif ii == 2
            title('CO');
        elseif ii == 3
            title('Other');
        elseif ii == 4
            ylabel('Gabor');
        elseif ii == 7
            ylabel('DOG');
        end
        figure(plot_counter+1); subplot(3,3,ii); imagesc(Output_List{idxs(ii),4}); set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray');
    end
    plot_counter = plot_counter + 2;
end

% Trying to figure out why some color cells are so well fit by gabor
% Selecting 10 color cells best fit by Gabor
debugfitquality = 0;
if debugfitquality
    ind = find(simplecells & ~Singleopponent);
    idx = ind(Color_GaborCrescent_svd(deltaBIC(Color_GaborCrescent_svd)>70));
    idx = ind(Scone_GaborCrescent_svd(deltaBIC(Scone_GaborCrescent_svd)>70));
%     idx = ind(Lum_GaborCrescent_svd(deltaBIC(Lum_GaborCrescent_svd)>100));
%     idx = ind(Other_GaborCrescent_svd(deltaBIC(Other_GaborCrescent_svd)>100));
    figure(plot_counter);
    for ii = 1:numel(idx)
        % STA
        tmp_vec_gun = Output_List{idx(ii),2};
        normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
        im = normfactor*tmp_vec_gun + 0.5;
        im = reshape(im,[10 10 3]);
        subplot(numel(idx),5,5*ii-4); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
        
        % Peak Frame
        subplot(numel(idx),5,5*ii-3); imagesc(Output_List{idx(ii),4}); set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
        
        % Crescent fit
        fittedCrescent = returnfittedimage(Crescentparams{idx(ii)},'Crescent');
        subplot(numel(idx),5,5*ii-2); imagesc(fittedCrescent); set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
        
        % Gabor fit
        fittedGabor = returnfittedimage(Gaborparams{idx(ii)},'Gabor');
        subplot(numel(idx),5,5*ii-1); imagesc(fittedGabor); set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
        
        % DoG fit
        fittedDOG = returnfittedimage(DOGparams{idx(ii)},'DOG');
        subplot(numel(idx),5,5*ii); imagesc(fittedDOG); set(gca,'XTick',[],'YTick',[]); colormap('gray'); axis square;
        
    end
    plot_counter = plot_counter + 1;
end

