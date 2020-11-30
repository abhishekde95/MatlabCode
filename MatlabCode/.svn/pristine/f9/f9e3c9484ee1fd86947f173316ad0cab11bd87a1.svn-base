% figures for reviewers' eyes only
% Author - Abhishek De, 4/20

close all; clearvars;
%% Figure 1: Show the STA, SVD, R, G and B significance maps for 10 cells that have the highest and lowest variance accounted for 
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
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

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
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.7; im(sigidxsneg+200) = 1;
    subplot(rows,cols,(count-1)*cols+2); image(im); axis image; set(gca,'XTick',[],'YTick',[]);

    % G-significance map
    pmat = logical(abs(zmat(:,:,2))>crit);
    im = zmat(:,:,2)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.7; im(sigidxsneg+200) = 1;
    subplot(rows,cols,(count-1)*cols+3); image(im); axis image; set(gca,'XTick',[],'YTick',[]); 
    
    % B-significance map
    pmat = logical(abs(zmat(:,:,3))>crit);
    im = zmat(:,:,3)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.7; im(sigidxsneg+200) = 1;
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

%% Figure 2: A k-means clustering on the RGB values within the RF (of significant pixels) and take the two vectors 

if ~exist('plot_counter')
    plot_counter = 1;
end
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
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

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

%% Figure 3: Plotting the feature vectors of simple and complex cells based on the NLI criterion
if ~exist('plot_counter')
    plot_counter = 1;
end

load Output_ListWN2.mat
load Singleopponent.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
FV_simple = [];
FV_complex = [];

for ii = 1:size(Output_List,1) 
   fr = Output_List{ii,24}; 
   bins = Output_List{ii,25};
   fr = (fr-min(fr));
   fr = fr/max(abs(fr));
   
   if fr(1)>fr(end)
       fr = fliplr(fr);
   end
   
   if simplecells(ii)
       FV_simple = [FV_simple; fr'];
   else 
       FV_complex = [FV_complex; fr'];
   end
 
end
subplot(121); errorbar(mean(FV_simple,1),std(FV_simple,1)/sqrt(size(FV_simple,1)),'-ko'); axis square; xlabel('proj'); ylabel('Normalized firing rate'); 
set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0); title('NLI<0'); text(3,0.75,strcat('N=',num2str(size(FV_simple,1))));
subplot(122); errorbar(mean(FV_complex,1),std(FV_complex,1)/sqrt(size(FV_complex,1)),'-ko'); axis square; xlabel('proj'); ylabel('Normalized firing rate'); 
set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0); title('NLI>0'); text(3,0.75,strcat('N=',num2str(size(FV_complex,1))));
plot_counter = plot_counter + 1;

%% Figure 4: Identifying the distribution of simple and DO cells for each monkey
if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);

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

MonkeyID = [];
for ii = 1:size(Output_List,1)
    MonkeyID = [MonkeyID; Output_List{ii,1}(1)];
end

Monkey_K = [sum(ismember(LumIds_conewts,find(MonkeyID=='K'))) sum(ismember(ColorOpponentIds_conewts,find(MonkeyID=='K'))) sum(ismember(Sconedominated_conewts,find(MonkeyID=='K'))) sum(ismember(Other_conewts,find(MonkeyID=='K')))];
Monkey_S = [sum(ismember(LumIds_conewts,find(MonkeyID=='S'))) sum(ismember(ColorOpponentIds_conewts,find(MonkeyID=='S'))) sum(ismember(Sconedominated_conewts,find(MonkeyID=='S'))) sum(ismember(Other_conewts,find(MonkeyID=='S')))];
Monkey_M = [sum(ismember(LumIds_conewts,find(MonkeyID=='M'))) sum(ismember(ColorOpponentIds_conewts,find(MonkeyID=='M'))) sum(ismember(Sconedominated_conewts,find(MonkeyID=='M'))) sum(ismember(Other_conewts,find(MonkeyID=='M')))];
Monkey_P = [sum(ismember(LumIds_conewts,find(MonkeyID=='P'))) sum(ismember(ColorOpponentIds_conewts,find(MonkeyID=='P'))) sum(ismember(Sconedominated_conewts,find(MonkeyID=='P'))) sum(ismember(Other_conewts,find(MonkeyID=='P')))];

%% Plotting the spatial RF and their ffts for all single opponent cells:
% Data driven approach 
% Ratio_of_power<2 
if ~exist('plot_counter')
    plot_counter = 1;
end
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

lambda = 12;
CRITERIA = Ratio_of_power<2 & simplecells;
N = find(CRITERIA)';
[~,I] = sort(Ratio_of_power(N));
N = N(I);
count = 1;
numplots = ceil(sqrt(numel(N)));
AllIdxs = ind([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]);
for ii = N 
    
    % STA 
    STA = Output_List{ii,2};
    normfactor = 0.5/(max(abs(STA(:)))+0.01);
    STA = normfactor*STA + 0.5;
    STA = reshape(STA,[10 10 3]);
    STA = sigmoid(STA,lambda/2,0.5);
    
    % Spatial RF
    im = Output_List{ii,4};
    fim = fftshift(fft2(im));
    
    % Modifying im for plotting purposes
    im = sigmoid(im,lambda/2,0)-0.5;
    
   
    if CRITERIA(ii)
        figure(plot_counter); subplot(numplots,numplots,count); image(STA); set(gca,'XTick',[],'YTick',[]); axis square; title(num2str(Ratio_of_power(ii),3));
        if ismember(ii,AllIdxs)
            set(gca,'XColor','r','YColor','r');
        end
        figure(plot_counter+1); subplot(numplots,numplots,count); imagesc(abs(fim)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title(num2str(Ratio_of_power(ii),3));
        if ismember(ii,AllIdxs)
            set(gca,'XColor','r','YColor','r');
        end
        count = count + 1;
    end
end
plot_counter = plot_counter + 2;