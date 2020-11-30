% This script is about analysis of magnitude of spatial opponency 
% Author - Abhishek De, 4/20

% Section  1: Analyses of spatial opponency of different cell types: Suggestion of reviewer 1 JoV

% Section 2: Analysis of one file for which Greg has the Gratings data 

%% Section  1: Analyses of spatial opponency of different cell types: Suggestion of reviewer 1 JoV

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
% Singleopponent = logical(Singleopponent);
Singleopponent = logical(Ratio_of_power<1.2);
SNR = Zmax(Z_cellsofinterest);


% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

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
% meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
% meanSSE = meanSSE(~Singleopponent & simplecells,:);
% meanR = meanR(~Singleopponent & simplecells,:);

% Need to assess the spatial opponecny of model fits 
MODELFITS = modelfits;
SPATIALRF = Output_List(:,4);
SOI_index = []; % Spatial opponency index 
SOI_index_svd = []; % Spatial opponency index using SVD spatial RF 
for ii=1:size(MODELFITS,1)
 
    % Based on Gabor model
    tmp = MODELFITS{ii,2}; % Gabor fits
    P = abs(sum(sum(tmp(tmp>0)))); N = abs(sum(sum(tmp(tmp<0))));
    SOI_index = [SOI_index; 1-abs((P-N)/(P+N))];
    
    % Based on SVD derived maps 
    tmp_svd = SPATIALRF{ii};
    P_svd = abs(sum(sum(tmp_svd(tmp_svd>0)))); N_svd = abs(sum(sum(tmp_svd(tmp_svd<0))));
    SOI_index_svd = [SOI_index_svd; 1-abs((P_svd-N_svd)/(P_svd+N_svd))];
    
end


% Plotting the results on SOI index that is derived from Gabor fits
GABORbetterthanDOG = meanR(:,2)>meanR(:,3);
bins = 0:0.1:1.0;
figure(plot_counter); set(gcf,'Name','Spatial opponency index: based on Gabor fits');
subplot(321); plot(SOI_index(Lumind),meanR(Lumind,2)-meanR(Lumind,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SOI_index(COind),meanR(COind,2)-meanR(COind,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SOI_index(Sconeind),meanR(Sconeind,2)-meanR(Sconeind,3),'o','MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[-.5 0.5]); xlabel('SOI'); ylabel('R(Gabor)-R(DoG)');
subplot(322); histogram(SOI_index(Lumind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index(Lumind(GABORbetterthanDOG(Lumind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index(Lumind)),30,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index(Lumind(GABORbetterthanDOG(Lumind)))),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 35]); xlabel('SOI'); ylabel('# cells'); title('LUM');
subplot(323); histogram(SOI_index(COind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index(COind(GABORbetterthanDOG(COind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index(COind)),12,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index(COind(GABORbetterthanDOG(COind)))),12,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 15]); xlabel('SOI'); ylabel('# cells'); title('DO');
subplot(324); histogram(SOI_index(Sconeind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
histogram(SOI_index(Sconeind(GABORbetterthanDOG(Sconeind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index(Sconeind)),9,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index(Sconeind(GABORbetterthanDOG(Sconeind)))),9,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 10]); xlabel('SOI'); ylabel('# cells'); title('S');
subplot(325); histogram(SOI_index(Singleopponent),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
plot(median(SOI_index(Singleopponent)),20,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
axis square; set(gca,'Tickdir','out','Ylim',[0 40]); xlabel('SOI'); ylabel('# cells'); title('Single-opponent (Ratio of power<1.2)');
subplot(326); histogram(SOI_index(Singleopponent),bins,'Displaystyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on; 
histogram(SOI_index([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2);
plot(median(SOI_index(Singleopponent)),30,'v','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
plot(median(SOI_index([Lumind; COind; Sconeind])),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('SOI'); ylabel('# cells'); legend('SO','LUM & DO'); title('All cells');
plot_counter = plot_counter + 1;


% Plotting the results on SOI index that is derived from SVD based spatial RF
GABORbetterthanDOG = meanR(:,2)>meanR(:,3);
bins = 0:0.1:1.0;
figure(plot_counter); set(gcf,'Name','Spatial opponency index: based on SVD spatial RF');
subplot(321); plot(SOI_index_svd(Lumind),meanR(Lumind,2)-meanR(Lumind,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SOI_index_svd(COind),meanR(COind,2)-meanR(COind,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SOI_index_svd(Sconeind),meanR(Sconeind,2)-meanR(Sconeind,3),'o','MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.3 1],'Ylim',[-.5 0.5]); xlabel('SOI'); ylabel('R(Gabor)-R(DoG)');
subplot(322); histogram(SOI_index_svd(Lumind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(Lumind(GABORbetterthanDOG(Lumind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(Lumind)),30,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(Lumind(GABORbetterthanDOG(Lumind)))),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 50]); xlabel('SOI'); ylabel('# cells'); title('LUM');
subplot(323); histogram(SOI_index_svd(COind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(COind(GABORbetterthanDOG(COind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(COind)),20,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(COind(GABORbetterthanDOG(COind)))),20,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30]); xlabel('SOI'); ylabel('# cells'); title('DO');
subplot(324); histogram(SOI_index_svd(Sconeind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
histogram(SOI_index_svd(Sconeind(GABORbetterthanDOG(Sconeind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(Sconeind)),9,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(Sconeind(GABORbetterthanDOG(Sconeind)))),9,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30]); xlabel('SOI'); ylabel('# cells'); title('S');
subplot(325); histogram(SOI_index_svd(Singleopponent),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
plot(median(SOI_index_svd(Singleopponent)),25,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30]); xlabel('SOI'); ylabel('# cells'); title('Single-opponent (Ratio of power<1.2)');
subplot(326); histogram(SOI_index_svd(Singleopponent),bins,'Displaystyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on; 
histogram(SOI_index_svd([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2);
plot(median(SOI_index_svd(Singleopponent)),30,'v','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
plot(median(SOI_index_svd([Lumind; COind; Sconeind])),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 80]); xlabel('SOI'); ylabel('# cells'); legend('SO','LUM & DO'); title('All cells');
plot_counter = plot_counter + 1;


% Checking if there are any correlatiosn
idx = [Lumind; COind; Sconeind]';
[r1,p1] = corr(SOI_index(idx),meanR(idx,2)-meanR(idx,3),'type','Spearman');
[r2,p2] = corr(SOI_index(idx),SNR(idx),'type','Spearman');

% Checking whether the SOI is different among the cell types 
group = [ones(size(Lumind)); 2*ones(size(COind)); 3*ones(size(Sconeind))];
data1 = SOI_index(idx);
p = kruskalwallis(data1,group,'off');


% Plotting the STAs of cells with SOI <0.2

N = [Lumind(SOI_index(Lumind)<0.2); COind(SOI_index(COind)<0.2); Sconeind(SOI_index(Sconeind)<0.2)];
numplots = ceil(sqrt(numel(N)));
count = 1;
lambda = 12;
for jj = 1:numel(N)
    % STA 
    ii = N(jj);
    STA = Output_List{ii,2};
    normfactor = 0.5/(max(abs(STA(:)))+0.01);
    STA = normfactor*STA + 0.5;
    STA = reshape(STA,[10 10 3]);
    
    % Spatial RF
    im = Output_List{ii,4};
    im = sigmoid(im,lambda,0)-0.5;
    
    % Gabor fit 
    fitim = modelfits{ii,2}; 
    fitim = sigmoid(fitim,lambda,0)-0.5;
    
    % Plotting ths STAs
    figure(plot_counter); subplot(numplots,numplots,count); image(STA); set(gca,'XTick',[],'YTick',[]); axis square; 
    % PLotting the spatial RF derived using SVD
    figure(plot_counter+1); subplot(numplots,numplots,count); imagesc(255*(im./(2*max(abs(im(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]);
    
    % PLotting Gabor
    figure(plot_counter+2); subplot(numplots,numplots,count); image(255*(fitim./(2*max(abs(fitim(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    count = count + 1;
    
end
plot_counter = plot_counter + 3;

%% Section 2: Analysis of one file for which Greg has the Gratings data 

if ~exist('plot_counter')
    plot_counter = 1;
end

filename = 'K042809001.nex'; 
II = 132;
lambda = 12;

% STA
STA = Output_List{II,2};
normfactor = 0.5/(max(abs(STA(:)))+0.01);
STA = normfactor*STA + 0.5;
STA = reshape(STA,[10 10 3]);
STA = sigmoid(STA,11,0.5);
    
% SVD derived spatial RF 
imRF = Output_List{II,4};
imRF = sigmoid(imRF,lambda,0)-0.5;

% Gabor fit 
fit = modelfits{II,2};
fit = sigmoid(fit,lambda,0)-0.5;

% fft of the Gabor fit 
fitfft = fftshift(fft2(fit));

figure(plot_counter); subplot(221); image(STA); set(gca,'XTick',[],'YTick',[]); axis square; title('STA')
subplot(222); imagesc(255*(imRF./(2*max(abs(imRF(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]); title('Spatial RF');
subplot(223); imagesc(255*(fit./(2*max(abs(fit(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]); title('Gabor fit');
subplot(224); imagesc(abs(fitfft)); colormap(gray); axis square; set(gca,'XTick',[],'YTick',[]); title(strcat('Ratio of power=',num2str(max(abs(fitfft(:)))/abs(fitfft(6,6)),3)))
