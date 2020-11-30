% Here I want to analyze the linear separability of space-color maps
% Author - Abhishek De, 4/20

% Section 1: PLotting the variance accounted for different cell types and calculating cone weights from WN checkerboard stim

% Section 2: Show the STA, SVD, R, G and B significance maps for 10 cells that have the highest and lowest variance accounted for 

% Section 3: A k-means clustering on the RGB values within the RF (of significant pixels) and take the two vectors

% Section 4: A simulation of variance accounted by the first eigenvector as a function of SNR for a perfectly separable and inseparable neuron

% Section 5: A simulation of k-means clustering analysis on 10x10x3 Gaussian random noise 

% Section 6: A simulation of assesing the correlation of color tuning for 2 random 3-D vectors 

% Section 7: Filenames for 1) a high SNR DO cell 2) a high SNR simple cell 3) a low SNR DO cell 4) a low SNR simple cell (asked by Greg)

% Section 8: Using pca and projection vectors to test separability 

%% Section 1: PLotting the variance accounted for different cell types and calculating cone weights from WN checkerboard stim

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
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
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% With all the 10x10 pixels included: Analyzing the LUM, DO and SO cells  
eigenvals = cell2mat(Output_List(:,17));
expvariance = eigenvals.^2./repmat(sum(eigenvals.^2,2),[1 3]);

% Within RF 
eigenvals_max = cell2mat(Output_List(:,19));
expvariance_max = eigenvals_max.^2./repmat(sum(eigenvals_max.^2,2),[1 3]); 

figure(plot_counter); set(gcf,'Name','% explained Variance');
subplot(221); errorbar((1:3)- 0.3,100*mean(expvariance(Lumind,:),1),100*std(expvariance(Lumind,:),1),'s','MarkerSize',12,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
errorbar((1:3)- 0.15,100*mean(expvariance(COind,:),1),100*std(expvariance(COind,:),1),'s','MarkerSize',12,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
errorbar((1:3),100*mean(expvariance(Sconeind,:),1),100*std(expvariance(Sconeind,:),1),'s','MarkerSize',12,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); 
errorbar((1:3)+0.15,100*mean(expvariance(Singleopponent,:),1),100*std(expvariance(Singleopponent,:),1),'s','MarkerSize',12,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); 
errorbar((1:3)+0.3,100*mean(expvariance(~simplecells,:),1),100*std(expvariance(~simplecells,:),1),'s','MarkerSize',12,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[0 4],'Ylim',[0 100],'XTick',[1 2 3],'YTick',[0:25:100]); axis square; xlabel('Singular vectors'); ylabel('% Variance'); 
legend('LUM','DO(LM)','DO(S)','SO','NLI>0'); title('Over all pixels');  hold off;

subplot(222); errorbar((1:3)- 0.3,100*nanmean(expvariance_max(Lumind,:),1),100*nanstd(expvariance_max(Lumind,:),1),'s','MarkerSize',12,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
errorbar((1:3)- 0.15,100*nanmean(expvariance_max(COind,:),1),100*nanstd(expvariance_max(COind,:),1),'s','MarkerSize',12,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
errorbar((1:3),100*nanmean(expvariance_max(Sconeind,:),1),100*nanstd(expvariance_max(Sconeind,:),1),'s','MarkerSize',12,'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]); 
errorbar((1:3)+0.15,100*nanmean(expvariance_max(Singleopponent,:),1),100*nanstd(expvariance_max(Singleopponent,:),1),'s','MarkerSize',12,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); 
errorbar((1:3)+0.3,100*nanmean(expvariance_max(~simplecells,:),1),100*nanstd(expvariance_max(~simplecells,:),1),'s','MarkerSize',12,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[0 4],'Ylim',[0 100],'XTick',[1 2 3],'YTick',[0:25:100]); axis square; xlabel('Singular vectors'); ylabel('% Variance'); 
legend('LUM','DO(LM)','DO(S)','SO','NLI>0'); title('Within RF');  hold off;

group1 = [ones(size(Lumind)); 2*ones(size(COind)); 3*ones(size(Sconeind))];
group2 = [group1; 4*ones(sum(Singleopponent,1),1); 5*ones(sum(~simplecells,1),1)];
data1 = expvariance([Lumind; COind; Sconeind]);
data2 = expvariance([Lumind; COind; Sconeind; find(Singleopponent); find(~simplecells)]);
p1 = kruskalwallis(data1,group1,'off');
p2 = kruskalwallis(data2,group2,'off');

% Need to think of a population analysis
% Something like difference in angles 
% STA - 2nd column, Spatial RF - 4th column

ANGLEDIFF_bw_posneg = []; % for storing the angular differences between the pos and neg subunits 
ANGLEDIFF_bw_conevectors = []; % for storing the angular differences between the SVD- and gun-map derived cone wts 
for ii = 1:size(Output_List,1)
   
   % Pulling out the STA
   imSTA = Output_List{ii,2};
   
   % Pulling out the SVD derived spatial RF
   imRF = Output_List{ii,4};
   
   % Pulling out the Mrgbtocc value
   Mrgbtocc = Output_List{ii,28};
   
   % RGB values of the increments
   RGBpos = sum(repmat(imRF(imRF>0),[1 3]).*imSTA(imRF>0,:),1);
   conewts_pos = Mrgbtocc*RGBpos';
   conewts_pos = conewts_pos/sum(abs(conewts_pos));
   
   % RGB values of the decrements 
   RGBneg = sum(repmat(-1*imRF(imRF<0),[1 3]).*imSTA(imRF<0,:),1);
   conewts_neg = Mrgbtocc*RGBneg';
   conewts_neg = conewts_neg/sum(abs(conewts_neg));
   
   DATADRIVEN_CONEWT_VECTOR = [conewts_pos; conewts_neg];
   DATADRIVEN_CONEWT_VECTOR = DATADRIVEN_CONEWT_VECTOR/norm(DATADRIVEN_CONEWT_VECTOR);
   
   
   % Actual cone weights
   wt = cell2mat(Output_List(ii,23)');
   wt = wt/sum(abs(wt));
   SVDCONEWT_VECTOR = [wt; -1*wt];
   SVDCONEWT_VECTOR = SVDCONEWT_VECTOR/norm(SVDCONEWT_VECTOR);
   
   % for storing the angular differences between the pos and neg subunits
   ANGLEDIFF_bw_posneg = [ANGLEDIFF_bw_posneg; acos(dot(conewts_pos,conewts_neg))*180/pi];
   
   % for storing the angular differences between the SVD- and gun-map derived cone wts
   ANGLEDIFF_bw_conevectors = [ ANGLEDIFF_bw_conevectors; 90-abs((acos(dot(DATADRIVEN_CONEWT_VECTOR,SVDCONEWT_VECTOR))*180/pi) - 90)];
   
end

figure(plot_counter); subplot(223); histogram(ANGLEDIFF_bw_posneg(Lumind),0:10:180,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(ANGLEDIFF_bw_posneg(COind),0:10:180,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
histogram(ANGLEDIFF_bw_posneg(Sconeind),0:10:180,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[0 180]); xlabel('Angle difference'); ylabel('# cells'); legend('LUM','DO(LM)','S'); 
title('Positive vs. Negative subunit'); hold off;

subplot(224); histogram(ANGLEDIFF_bw_conevectors(Lumind),0:5:90,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(ANGLEDIFF_bw_conevectors(COind),0:5:90,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
histogram(ANGLEDIFF_bw_conevectors(Sconeind),0:5:90,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2);
axis square; set(gca,'Tickdir','out','Xlim',[0 90]); xlabel('Angle difference'); ylabel('# cells'); legend('LUM','DO(LM)','S'); 
title('Cone wt vs. prediction from RGB maps'); hold off;

plot_counter= plot_counter + 1;

%% Section 2: Show the STA, SVD, R, G and B significance maps for 10 cells that have the highest and lowest variance accounted for 

IDX = [Lumind; COind; Sconeind];
[~,I] = sort(expvariance(IDX,1));
numcells = 10;
HI = IDX(I(end:-1:end-numcells+1));
LO = IDX(I(1:numcells));
rows = numel(HI);
cols = 6;
lambda = 12;
count = 1;
IDXS = [LO; HI];

for ii = 1:numel(IDXS)
    jj = IDXS(ii);
    %STA
    im = Output_List{jj,2};
    normfactor = 0.5/(max(abs(im(:)))+0.01);
    im = normfactor*im + 0.5;
    im = reshape(im,[10 10 3]);
    figure(plot_counter);
    subplot(rows,cols,(count-1)*cols+1); image(im); set(gca,'XTick',[],'YTick',[]); axis square; title(num2str(100*expvariance(jj),3));
    
    % Calculation of significance maps
    im = reshape(Output_List{jj,2},[10 10 3]);
    mumat = reshape(Output_List{jj,26},[10 10 3]);
    sigmamat = reshape(Output_List{jj,27},[10 10 3]);
    nspikes = Output_List{jj,15};
    zmat = (im-mumat)./(sigmamat/sqrt(nspikes));
    maxzscore = max(abs([zmat(:)]));
    alpha = 0.05;
    crit = norminv(1-alpha,0,1);
    
    % R-significance map
    pmat = logical(abs(zmat(:,:,1))>crit);
    im = zmat(:,:,1)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.5; im(sigidxsneg+200) = 1;
    subplot(rows,cols,(count-1)*cols+2); image(im); axis image; set(gca,'XTick',[],'YTick',[]);

    % G-significance map
    pmat = logical(abs(zmat(:,:,2))>crit);
    im = zmat(:,:,2)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.5; im(sigidxsneg+200) = 1;
    subplot(rows,cols,(count-1)*cols+3); image(im); axis image; set(gca,'XTick',[],'YTick',[]); 
    
    % B-significance map
    pmat = logical(abs(zmat(:,:,3))>crit);
    im = zmat(:,:,3)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.5; im(sigidxsneg+200) = 1;
    % red to .5 where sig.  Looks red on dark and cyan on bright.
    subplot(rows,cols,(count-1)*cols+4); image(im); axis image; set(gca,'XTick',[],'YTick',[]); 
    
    % Cone weights
    wt = cell2mat(Output_List(jj,23)');
    wt = sign(wt(2))*wt;
    wt = wt/sum(abs(wt));
    subplot(rows,cols,(count-1)*cols+5); bar(wt','FaceColor','k'); set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;
    
    % Spatial weighting function
    im = Output_List{jj,4};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(rows,cols,(count-1)*cols+6); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    count  = count + 1;
    
    if count ==rows+1
        plot_counter = plot_counter + 1;
        count = 1;
    end
     
end
plot_counter = plot_counter + 1;

%% Section 3: A k-means clustering on the RGB values within the RF (of significant pixels) and take the two vectors 

if ~exist('plot_counter')
    plot_counter = 1;
end

CORR_SUB12_withinRF = [];
CORR_SUB12_allpixels = [];

for ii = 1:size(Output_List,1)
    sigRF = Output_List{ii,18};
    STA = Output_List{ii,2};
    
    % analysis of all the pixels
    clusterID_all = kmeans(STA,2);
    subunit1RGB_all = mean(STA(clusterID_all==1,:),1);
    subunit2RGB_all = mean(STA(clusterID_all==2,:),1);
    
    tmpcorr_all = corr(subunit1RGB_all',subunit2RGB_all');
    
    % analysis within the RF
    if ~any(isnan(sigRF(:))) 
        RGB = STA(sigRF(:),:);
        clusterID = kmeans(RGB,2);
        
        subunit1RGB = mean(RGB(clusterID==1,:),1);
        subunit2RGB = mean(RGB(clusterID==2,:),1);
        
        tmpcorr = corr(subunit1RGB',subunit2RGB');
    else
        tmpcorr = nan;
    end
    
    CORR_SUB12_allpixels = [CORR_SUB12_allpixels; tmpcorr_all];
    CORR_SUB12_withinRF = [CORR_SUB12_withinRF; tmpcorr];
end


figure(plot_counter); set(gcf,'Name','Within RF');
subplot(221); histogram(CORR_SUB12_withinRF([Lumind;COind;Sconeind]),-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(CORR_SUB12_withinRF(Singleopponent),-1:0.1:1,'Displaystyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2);
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); legend('LUM & DO','SO'); hold off;
subplot(222); histogram(CORR_SUB12_withinRF([Lumind]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title('LUM'); hold off;
subplot(223); histogram(CORR_SUB12_withinRF([COind]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title('DO(LM)'); hold off;
subplot(224); histogram(CORR_SUB12_withinRF([Sconeind]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title('DO(S)'); hold off;
plot_counter = plot_counter + 1;

figure(plot_counter); set(gcf,'Name','All pixels');
subplot(221); histogram(CORR_SUB12_allpixels([Lumind;COind;Sconeind]),-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(CORR_SUB12_allpixels(Singleopponent),-1:0.1:1,'Displaystyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2);
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); legend('LUM & DO','SO'); hold off;
subplot(222); histogram(CORR_SUB12_allpixels([Lumind]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title('LUM'); hold off;
subplot(223); histogram(CORR_SUB12_allpixels([COind]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title('DO(LM)'); hold off;
subplot(224); histogram(CORR_SUB12_allpixels([Sconeind]),-1:0.1:1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title('DO(S)'); hold off;
plot_counter = plot_counter + 1;

%% Section 4: A simulation of variance accounted by the first eigenvector as a function of SNR for a perfectly separable and inseparable neuron

if ~exist('plot_counter')
    plot_counter = 1;
end

% A perfectly seprabale neuron
Separableneuron = [1 -1 -1;...
    1 -1 -1];

% A perfectly inseparable neuron
Inseparableneuron = [1 0 -1;...
    0.5 1 0.5];

N = 41;
std_noise = logspace(-4,4,N);

SNR_separable = [];
SNR_inseparable = [];
Variance_explained_separable = [];
Variance_explained_inseparable = [];
iter = 1000;
for ii = 1:N
    
    % Analysis of a separable neuron
    SNR_separable = [SNR_separable; sum(sum(Separableneuron.^2/std_noise(ii)))]; % SNR is defined is the square of z-score 
    
    % Analysis of a separable neuron
    SNR_inseparable = [SNR_inseparable; sum(sum(Inseparableneuron.^2/std_noise(ii)))]; 
    
    VS = []; % Variance separable
    VIS = []; % Variance inseparable
    for jj = 1:iter
        tmp1 = Separableneuron + randn([2 3])*std_noise(ii);
        [~,s1,~] = svd(tmp1);
        VS = [VS; 100*(s1(1,1)^2/(s1(1,1)^2+s1(2,2)^2))];
        
        tmp2 = Inseparableneuron + randn([2 3])*std_noise(ii);
        [~,s2,~] = svd(tmp2);
        VIS = [VIS; 100*(s2(1,1)^2/(s2(1,1)^2+s2(2,2)^2))];
    end
    
    Variance_explained_separable = [Variance_explained_separable; mean(VS)];
    Variance_explained_inseparable = [Variance_explained_inseparable; mean(VIS)];
    
end

figure(plot_counter); 
subplot(121); plot(SNR_separable,Variance_explained_separable,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'Tickdir','out','XScale','log','Ylim',[50 100],'Xlim',[0.0001 10000]); axis square; xlabel('SNR'); ylabel('% Variance explained'); title('Perfectly Separable neuron'); hold off;
subplot(122); plot(SNR_inseparable,Variance_explained_inseparable,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
set(gca,'Tickdir','out','XScale','log','Ylim',[50 100],'Xlim',[0.0001 10000]); axis square; xlabel('SNR'); ylabel('% Variance explained'); title('Perfectly inseparable neuron'); hold off;
plot_counter = plot_counter + 1;

%% Section 5: A simulation of k-means clustering analysis on 10x10x3 Gaussian random noise 

if ~exist('plot_counter')
    plot_counter = 1;
end

iter = 1000; 
dimensions = [3 100 1000];

for jj = 1:numel(dimensions)
    dim = dimensions(jj);
    CORR_SUB12_allpixels = [];
    CORR_SUB12_allpixels_shift = [];
    for ii = 1:iter
        STA = randn(100,dim);
        clusterID_all = kmeans(STA,2);
        
        STAshift = STA+repmat(randn(1,dim),[100 1]);
        clusterID_all_shift = kmeans(STAshift,2);
        
        subunit1RGB_all = mean(STA(clusterID_all==1,:),1);
        subunit2RGB_all = mean(STA(clusterID_all==2,:),1);
        
        subunit1RGB_all_shift = mean(STAshift(clusterID_all_shift==1,:),1);
        subunit2RGB_all_shift = mean(STAshift(clusterID_all_shift==2,:),1);
        
        tmpcorr_all = corr(subunit1RGB_all',subunit2RGB_all');
        tmpcorr_all_shift = corr(subunit1RGB_all_shift',subunit2RGB_all_shift');
        
        CORR_SUB12_allpixels = [CORR_SUB12_allpixels; tmpcorr_all];
        CORR_SUB12_allpixels_shift = [CORR_SUB12_allpixels_shift; tmpcorr_all_shift];
    end
    
    figure(plot_counter);
    subplot(3,2,2*(jj-1)+1); histogram(CORR_SUB12_allpixels,-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
    axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title(strcat('Zero mean Gaussian noise, dim=',num2str(dim))); hold off;
    subplot(3,2,2*(jj-1)+2); histogram(CORR_SUB12_allpixels_shift,-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
    axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title(strcat('Shifted mean Gaussian noise, dim=',num2str(dim))); hold off;
end
plot_counter = plot_counter + 1;

%% Section 6: A simulation of assesing the correlation of color tuning for 2 random 3-D vectors 

if ~exist('plot_counter')
    plot_counter = 1;
end

iter = 10000; 
dimensions = [3 100 1000];
dim = dimensions(1);
CORR_SUB12_allpixels = [];

for ii = 1:iter
    
    tmpcorr_all = corr(randn(3,1),randn(3,1));
    
    CORR_SUB12_allpixels = [CORR_SUB12_allpixels; tmpcorr_all];
    
end

figure(plot_counter);
histogram(CORR_SUB12_allpixels,-1:0.1:1,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
axis square; set(gca,'Tickdir','out'); xlabel('Correlation'); ylabel('# cells'); title(strcat('Zero mean Gaussian noise, dim=',num2str(dim))); hold off;
plot_counter = plot_counter + 1;

%% Section 7: Filenames for 1) a high SNR DO cell 2) a high SNR simple cell 3) a low SNR DO cell 4) a low SNR simple cell (asked by Greg)

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load Ratio_of_power.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
SNR = Zmax(Z_cellsofinterest);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% Calculating the indices of the highest and lowest SNR cells
[~,Sim1] = max(SNR(Lumind)); [~,Sim2] = min(SNR(Lumind)); % Simple cells
[~,DLM1] = max(SNR(COind)); [~,DLM2] = min(SNR(COind)); % DO LM cells
[~,DS1] = max(SNR(Sconeind)); [~,DS2] = min(SNR(Sconeind)); % DO LM cells
idx = [Lumind(Sim1); Lumind(Sim2); COind(DLM1); COind(DLM2); Sconeind(DS1); Sconeind(DS2)];


% Plotting the STA of the highest and lowest SNR of simple, DO(L-M) and DO(S) cells 
figure(plot_counter);
count = 1;
rows = numel(idx);
cols = 4;
for ii = 1:numel(idx)
    jj = idx(ii);
    %STA
    im = Output_List{jj,2};
    normfactor = 0.5/(max(abs(im(:)))+0.01);
    im = normfactor*im + 0.5;
    im = reshape(im,[10 10 3]);
    figure(plot_counter);
    subplot(rows,cols,(count-1)*cols+1); image(im); set(gca,'XTick',[],'YTick',[]); axis square; title(num2str(SNR(jj),3));
    
    % Calculation of significance maps
    im = reshape(Output_List{jj,2},[10 10 3]);
    mumat = reshape(Output_List{jj,26},[10 10 3]);
    sigmamat = reshape(Output_List{jj,27},[10 10 3]);
    nspikes = Output_List{jj,15};
    zmat = (im-mumat)./(sigmamat/sqrt(nspikes));
    maxzscore = max(abs([zmat(:)]));
    alpha = 0.05;
    crit = norminv(1-alpha,0,1);
    
    % R-significance map
    pmat = logical(abs(zmat(:,:,1))>crit);
    im = zmat(:,:,1)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.5; im(sigidxsneg+200) = 1;
    subplot(rows,cols,(count-1)*cols+2); image(im); axis image; set(gca,'XTick',[],'YTick',[]);
    if count == 1
        title('R');
    end

    % G-significance map
    pmat = logical(abs(zmat(:,:,2))>crit);
    im = zmat(:,:,2)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.5; im(sigidxsneg+200) = 1;
    subplot(rows,cols,(count-1)*cols+3); image(im); axis image; set(gca,'XTick',[],'YTick',[]); 
    if count == 1
        title('G');
    end
    
    % B-significance map
    pmat = logical(abs(zmat(:,:,3))>crit);
    im = zmat(:,:,3)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.5; im(sigidxsneg+200) = 1;
    % red to .5 where sig.  Looks red on dark and cyan on bright.
    subplot(rows,cols,(count-1)*cols+4); image(im); axis image; set(gca,'XTick',[],'YTick',[]);
    if count == 1
        title('B');
    end
    
    count = count + 1;
end
plot_counter = plot_counter + 1;

%% Section 8: Using pca and projection vectors to test separability 

close all; clearvars;
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load Ratio_of_power.mat

crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
SNR = Zmax(Z_cellsofinterest);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

Separability_index_1 = []; % Adding a separability index, includes the calculation of angle between the pca and the projection vectors 
Separability_index_2 = []; % Adding a separability index, DOES NOT include ...
Separability_index_3 = []; % Separability index, NULL distribution for Separability_index_1 
Separability_index_4 = []; % Separability index, NULL distribution for Separability_index_2
Varianceaccounted = []; % Difference in variance between the two models, quantified as a percentage

nstixperside = 10; % # KIND OF A HACK
for ii = 1:size(Output_List,1)
    disp(ii);
    reshapedweightedSTA = Output_List{ii,2};
   
    coef = pca(reshapedweightedSTA); 
       
    v1 = coef(:,1);
    v2 = -coef(:,1);
    
    myvariance(reshapedweightedSTA,[v1 v2]); % The mean squared projection onto the PC1
    f = @(x)-1*myvariance(reshapedweightedSTA, x); % Setting up anonymous function for fmincon
    [vs,fval] = fmincon(f,[v1 v2],[],[],[],[],[],[],@unitnormcon);
    
    vec1_pca = coef(:,1)/norm(coef(:,1));
    vec1_proj = vs(:,1)/norm(vs(:,1));
    vec2_proj = vs(:,2)/norm(vs(:,2));
    
    tmp_SI = abs(dot(vec1_proj,vec2_proj))*max([abs(dot(vec1_pca,vec1_proj)) abs(dot(vec1_pca,vec2_proj))]);
    Separability_index_1 = [Separability_index_1; tmp_SI];
    
    tmp_SII = abs(dot(vec1_proj,vec2_proj));
    Separability_index_2 = [Separability_index_2; tmp_SII];
    
    % Calculating the variance accounted for by the vectors
    var_pca = myvariance(reshapedweightedSTA,[v1 v2]); % The mean squared projection onto the PC1
    var_twovector = myvariance(reshapedweightedSTA,[vec1_proj vec2_proj]); % The mean squared projection onto the PC1
    
    Varianceaccounted = [Varianceaccounted; 100*(var_pca-var_twovector)/var_pca];
    
    
    % Creating a null distribution for the separability indices and
    % repeating the analysis again: BAD coding style 
    STAs = reshapedweightedSTA(:);
    idxs = [randperm(nstixperside.^2)' randperm(nstixperside.^2)' randperm(nstixperside.^2)'];
    idxs = idxs+repmat([0 nstixperside.^2 2*nstixperside.^2],nstixperside.^2,1);
    STAs = reshape(STAs(idxs),[nstixperside.^2 3]);
    
    coef2 = pca(STAs); 
       
    v3 = coef2(:,1);
    v4 = -coef2(:,1);
   
    f = @(x)-1*myvariance(STAs, x); % Setting up anonymous function for fmincon
    [vs2,fval] = fmincon(f,[v3 v4],[],[],[],[],[],[],@unitnormcon);
    
    vec1_pca = coef2(:,1)/norm(coef2(:,1));
    vec1_proj = vs2(:,1)/norm(vs2(:,1));
    vec2_proj = vs2(:,2)/norm(vs2(:,2));
    
    tmp_SI = abs(dot(vec1_proj,vec2_proj))*max([abs(dot(vec1_pca,vec1_proj)) abs(dot(vec1_pca,vec2_proj))]);
    Separability_index_3 = [Separability_index_3; tmp_SI];
    
    tmp_SII = abs(dot(vec1_proj,vec2_proj));
    Separability_index_4 = [Separability_index_4; tmp_SII];
   
end

bins = 0:0.05:1.0;
figure(plot_counter); 
subplot(221); histogram(Separability_index_1([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(Separability_index_3([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); 
xlabel('Separability index 1'); ylabel('# cells'); set(gca,'Tickdir','out'); axis square; legend('Simple & DO','NULL dist');
subplot(222); histogram(Separability_index_2([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on; 
histogram(Separability_index_4([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2);
xlabel('Separability index 2'); ylabel('# cells'); set(gca,'Tickdir','out'); axis square; legend('Simple & DO','NULL dist');
subplot(223); histogram(Separability_index_1(Lumind),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(Separability_index_1(COind),bins,'Displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
histogram(Separability_index_1(Sconeind),bins,'Displaystyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2);
xlabel('Separability index 1'); ylabel('# cells'); set(gca,'Tickdir','out'); axis square; legend('Simple','DO(LM)','DO(S)');
subplot(224); histogram(Separability_index_2(Lumind),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(Separability_index_2(COind),bins,'Displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
histogram(Separability_index_2(Sconeind),bins,'Displaystyle','stairs','EdgeColor',[0 0.5 1.0],'Linewidth',2);
xlabel('Separability index 2'); ylabel('# cells'); set(gca,'Tickdir','out'); axis square; legend('Simple','DO(LM)','DO(S)');
plot_counter = plot_counter + 1;

% Doing some stats (just looking at the Separability index 2 and its NULL distribution
p1 = ranksum(Separability_index_2([Lumind; COind; Sconeind]), Separability_index_4([Lumind; COind; Sconeind]));
p2 = ranksum(Separability_index_2([Lumind]), Separability_index_4([Lumind]));
p3 = ranksum(Separability_index_2([COind]), Separability_index_4([COind]));
p4 = ranksum(Separability_index_2([Sconeind]), Separability_index_4([Sconeind]));

% Note: this function does not subtract the mean from the data, which makes
% sense in this context, so the "variance" it returns is slightly different 
% from the variance of projections onto the PC1, even when "basis" = [PC1
% -PC1].
function v = myvariance(data, basis)
projs = data*basis;
projs(projs<0) = 0;
projs = max(projs,[],2);
v = mean(projs.^2);
end

function [c,ceq] = unitnormcon(x)
c = [];
ceq = sum(x.^2)-1;
end





