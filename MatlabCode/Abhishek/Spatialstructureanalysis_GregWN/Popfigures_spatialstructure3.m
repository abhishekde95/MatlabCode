% A replica of Popfigures_spatialstructure2 with minor corrections. Writing this script after JoV rejection. 
% Author - Abhishek De, 8/20
% Sending the paper to J Neurophysiology 
% Also check more codes on Popfigures_spatialstructure3.m

close all; clearvars;
plot_counter = 1;

%% Figure 1: Derivation of cone weights and spatial weighting functions
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

figure(plot_counter); set(gcf,'Name','Stim->STA');
for ii = 1:3
    im = rand(10,10,3);
    avg = mean(mean(im));
    im = im - repmat(avg,[10 10 1]);
    subplot(2,4,ii),image(im./(6*max(abs(im(:))))+.5); set(gca,'XTick',[],'YTick',[]); axis square; 
end
idx = COind(23);
tmp_vec_gun = Output_List{idx,2};
normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
im = normfactor*tmp_vec_gun + 0.5;
im = reshape(im,[10 10 3]);
subplot(2,4,4); image(im); set(gca,'XTick',[],'YTick',[]); axis square;

% Next, I am going to convert STA to cone and spatial weighting functions
lambda = 12;
im = Output_List{idx,4};
im = sigmoid(im,lambda,0)-0.5;
subplot(2,4,5); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
wt = Output_List{idx,23}; wt = sign(wt(2))*wt;
wt = wt/sum(abs(wt));
subplot(2,4,6); bar(wt,'FaceColor','k'); set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;


% Calculating the residue
STA = Output_List{idx,2}; 
STA = STA/norm(STA(:));
color_wts = Output_List{idx,5};
spatial_wts = Output_List{idx,4};
predictedSTA = spatial_wts(:)*color_wts'; 
predictedSTA = predictedSTA/norm(predictedSTA(:));
residue = STA - predictedSTA;

% plotting the result
normfactor = norm(residue(:))*0.5/(max(abs(residue(:)))+0.01);
residue = 0.8*normfactor*residue + 0.5;
residue = reshape(residue,[10 10 3]);
normfactor = norm(predictedSTA(:))*0.5/(max(abs(predictedSTA(:)))+0.01);
predictedSTA = 0.7*normfactor*predictedSTA + 0.5;
predictedSTA = reshape(predictedSTA,[10 10 3]);
subplot(247);imagesc(residue); set(gca,'XTick',[],'YTick',[]); axis square;
subplot(248);imagesc(predictedSTA); set(gca,'XTick',[],'YTick',[]); axis square;
plot_counter = plot_counter + 1;

% Now quantifying the per variance captured by the eigenevectors
eigenvals = cell2mat(Output_List([Lumind;COind;Sconeind],17));
expvariance = eigenvals.^2./repmat(sum(eigenvals.^2,2),[1 3]); 
mean_expvariance = 100*mean(expvariance,1);
std_expvariance = 100*std(expvariance,1);
eigenvals_max = cell2mat(Output_List([Lumind;COind;Sconeind],19));
expvariance_max = eigenvals_max.^2./repmat(sum(eigenvals_max.^2,2),[1 3]); 
mean_expvariance_max = 100*nanmean(expvariance_max,1);
std_expvariance_max = 100*nanstd(expvariance_max,1);
figure(plot_counter); set(gcf,'Name','% explained Variance');
errorbar((1:3)+ 0.2,mean_expvariance_max,std_expvariance_max,'s','MarkerSize',12,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
errorbar((1:3),mean_expvariance,std_expvariance,'o','MarkerSize',10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot((1:3) - 0.2,100*Output_List{idx,17}.^2/sum(Output_List{idx,17}.^2),'o','MarkerSize',10,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[0 4],'Ylim',[0 100],'XTick',[1 2 3],'YTick',[0:25:100]); axis square; xlabel('eigenvectors'); ylabel('% Variance'); hold off;
plot_counter = plot_counter + 1;

% Some additional analyses on the number of spikes and the eccentricities of the neurons that we considered 
ecc = cell2mat(Output_List([Lumind;COind;Sconeind],8));
ecc = sqrt(sum(ecc.^2,2))/10;
spikes = cell2mat(Output_List([Lumind;COind;Sconeind],15));

% Also plotting the other singular vectors - Would be good to see what they look like
[u,s,v] = svd(STA');
vec2 = v(:,2)*u(:,2)'; vec2 = vec2(:);
vec3 = v(:,3)*u(:,3)'; vec3 = vec3(:);
vec2 = (0.5*0.5*vec2./(max(abs(vec2))+eps)) + 0.5;
vec3 = (0.5*0.5*vec3./(max(abs(vec3))+eps)) + 0.5;
figure(plot_counter); 
subplot(221); image(predictedSTA); set(gca,'XTick',[],'YTick',[]); axis square;
subplot(222); image(reshape(vec2,[10 10 3])); set(gca,'XTick',[],'YTick',[]); axis square;
subplot(223); image(reshape(vec3,[10 10 3])); set(gca,'XTick',[],'YTick',[]); axis square;
plot_counter = plot_counter + 1;

%% Figure 2: Plot of cone weights of all the screened cells
% Calculating conewts of the cells of interest: SVD 
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


figure(plot_counter); set(gcf,'Name','Cone wts') 
plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

% Further classification of S-cone sensitive cells
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
SvsLM = Sconesensitive(:,sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==-1);
SMvsL = Sconesensitive(:,sign(Sconesensitive(1,:))==-1 & sign(Sconesensitive(3,:))==1);
SLvsM = Sconesensitive(:,sign(Sconesensitive(1,:))==-1 & sign(Sconesensitive(3,:))==-1);

%% Figure 3: Model fits for each cell type
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
ind = find(~Singleopponent & simplecells);

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
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% idx = [Lumind(8); COind(30); Sconeind(20)];
idx = [Lumind(18); Lumind(2); COind(23); COind(18); Sconeind(22); Sconeind(11)];
figure(plot_counter);
L = numel(idx);
swt_tmp = [-1 -1 1 1 -1 1];
gf_tmp = [-1 1 1 -1 -1 1];
dogf_tmp = [-1 1 1 1 -1 -1];
lambda = 12;
for ii = 1:numel(idx)
    jj = idx(ii);
    %STA
    im = Output_List{jj,2};
    normfactor = 0.5/(max(abs(im(:)))+0.01);
    im = normfactor*im + 0.5;
    im = reshape(im,[10 10 3]);
    subplot(8,L,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    
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
    subplot(8,L,ii+L); image(im); axis image; set(gca,'XTick',[],'YTick',[]);

   
    % G-significance map
    pmat = logical(abs(zmat(:,:,2))>crit);
    im = zmat(:,:,2)./(2*maxzscore)+.5;
    pmat = pmat.*sign(im-0.5);
    im = repmat(im,[1 1 3]);
    im = sigmoid(im,10,0.5);
    sigidxspos = find(pmat==1); sigidxsneg = find(pmat==-1);
    im(sigidxspos) = 1;  im(sigidxspos+100) = 0.2; im(sigidxspos+200) = 0;
    im(sigidxsneg) = 0;  im(sigidxsneg+100) = 0.7; im(sigidxsneg+200) = 1;
    subplot(8,L,ii+2*L); image(im); axis image; set(gca,'XTick',[],'YTick',[]); 
    
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
    subplot(8,L,ii+3*L); image(im); axis image; set(gca,'XTick',[],'YTick',[]); 
    
    % Cone weights
    wt = cell2mat(Output_List(jj,23)');
    wt = sign(wt(2))*wt;
    wt = wt/sum(abs(wt));
    subplot(8,L,ii+4*L); bar(wt','FaceColor','k'); set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;
    
    % Spatial weighting function
    im = swt_tmp(ii)*Output_List{jj,4};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(8,L,ii+5*L); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    
    % Gabor fit
    im = gf_tmp(ii)*modelfits{jj,2};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(8,L,ii+6*L);image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    
    % DoG fit
    im = -1*dogf_tmp(ii)*modelfits{jj,3};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(8,L,ii+7*L); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
end
plot_counter = plot_counter + 1;

% Adding some new results: AD, 8/2020
% Loading R_square and Pearson R for the entire dataset
load R_square.mat
load Pearson_R.mat

% Calculating the fraction of unexplained variance for Gabor and DOG models
FUV_Gabor = 1-R_square([Lumind; COind; Sconeind],2);
FUV_DOG = 1-R_square([Lumind; COind; Sconeind],3);

% Calculating the Pearson correlation coefficient for Gabor and DOG models
Pearson_Gabor = abs(cos(Pearson_R([Lumind; COind; Sconeind],2)));
Pearson_DOG = abs(cos(Pearson_R([Lumind; COind; Sconeind],3)));

%% Figure 4: Model comparison, Gabor vs DOG  
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
SNR = Zmax(Z_cellsofinterest)./Z(Z_cellsofinterest,1);
SNR = SNR(~Singleopponent & simplecells);

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

idx = [18; 2; 23; 18; 22; 11];
bins2 = -10:1:10;
% Based on CV-R
figure(plot_counter); set(gcf,'Name','Model Comparison: DoG vs Gabor');
subplot(311);  hold on;
for ii = 1:length(LumIds_conewts)
    if ii ~=18 | ii ~=2
        h(ii) = plot(meanR(LumIds_conewts(ii),2),meanR(LumIds_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 18
        h(ii) = plot(meanR(LumIds_conewts(ii),2),meanR(LumIds_conewts(ii),3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 2
        h(ii) = plot(meanR(LumIds_conewts(ii),2),meanR(LumIds_conewts(ii),3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
end
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');

subplot(312), hold on;
for ii = 1:length(ColorOpponentIds_conewts)
    if ii ~=23 | ii ~=18
        h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),2),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 23
        h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),2),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 18
        h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),2),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
end
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

subplot(313), hold on;
for ii = 1:length(Sconedominated_conewts)
    if ii ~=22 | ii ~=11
        h(ii) = plot(meanR(Sconedominated_conewts(ii),2),meanR(Sconedominated_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 22
        h(ii) = plot(meanR(Sconedominated_conewts(ii),2),meanR(Sconedominated_conewts(ii),3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 11
        h(ii) = plot(meanR(Sconedominated_conewts(ii),2),meanR(Sconedominated_conewts(ii),3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
end
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('S');
plot_counter = plot_counter + 1;

% Some more analyses for comparing the results between Lum, L-M and S cells: Gabor vs DOG
% sign rank sum
p1 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3));
p2 = signrank(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,3));
p3 = signrank(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,3));

% Mann-Whitney U test
p4 = ranksum(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,3),meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3));
p5 = ranksum(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,3),meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,3));
p6 = ranksum(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3),meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,3));

% Kruskal Wallis test for comparing the difference in R
data1 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data2 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p7 = kruskalwallis(data1-data2,group,'off');

% Pooling the DO cells to check whether the p<0.05, A concern brought to
% attention by Reviewer 1
DOwts = [ColorOpponentIds_conewts Sconedominated_conewts];
p8 = signrank(meanR(DOwts,2),meanR(DOwts,3)); % Wilcoxon signed test
[~,p9] = ttest(meanR(DOwts,2)-meanR(DOwts,3)); % T-test
[~,p10] = ttest(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3)); % T-test

%% Figure 5: Plotting R square as a function of SNR (energy of peak frame/energy of noise frame)
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

SNR = Zmax(Z_cellsofinterest)/300; % Dividing by the number of elements
SNR = SNR(~Singleopponent & simplecells);

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

bins = 0:0.1:1.2;
figure(plot_counter); set(gcf,'Name','R=f(SNR)');
subplot(121); hold on;
plot(SNR(ColorOpponentIds_conewts),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(SNR(LumIds_conewts),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1],'XTick',[1 2 4 8 16],'YTick',0:0.25:1.0, 'XScale','log'); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(122); hold on;
plot(SNR(ColorOpponentIds_conewts),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(SNR(LumIds_conewts),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1],'XTick',[1 2 4 8 16],'YTick',0:0.25:1.0,'XScale','log'); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Kruskal Wallis test for comparing SNR between groups 
data1 = [log10(SNR(LumIds_conewts)); log10(SNR(ColorOpponentIds_conewts)); log10(SNR(Sconedominated_conewts))];
data2 = [meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]; % Gabor 
data3 = [meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]; % DoG
group = [ones(size(log10(SNR(LumIds_conewts)))); 2*ones(size(log10(SNR(ColorOpponentIds_conewts)))); 3*ones(size(log10(SNR(Sconedominated_conewts))))];
p1 = kruskalwallis(data1,group,'off');
p2 = kruskalwallis(data2,group,'off');
p3 = kruskalwallis(data3,group,'off');
idxs = data1>median(data1);
p4 = kruskalwallis(data3(idxs),group(idxs),'off');
p5 = kruskalwallis(data2(idxs),group(idxs),'off');

% Z-transforming the R values
data1 = [SNR(LumIds_conewts); SNR(ColorOpponentIds_conewts); SNR(Sconedominated_conewts)];
data1log = log10(data1);
data2 = atanh([meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]); % Gabor 
data3 = atanh([meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]); % DoG
group = [ones(size(SNR(LumIds_conewts))); 2*ones(size(SNR(ColorOpponentIds_conewts))) ; 3*ones(size(SNR(Sconedominated_conewts)))];

% Need to quantify the the differences in DoG as a function of SNR:
celltype_dv = zeros(numel(group),2);
for ii = 1:size(celltype_dv,1)
    celltype_dv(ii,group(ii)) = 1;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
colors = {'black','red','cyan'};

numcells = size(celltype_dv,2);
 
X = [repmat(data1log,1,numcells).*celltype_dv celltype_dv]; % Individual model
[b, bint, r_Gabor,~, stats] = regress(data2,X);
[b_full, bint_full, r_full_Gabor] = regress(data2,[data1log ones(size(data1log))]);

% Trying linear regression
figure(plot_counter);
subplot(121); plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,2)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SNR(ColorOpponentIds_conewts),atanh(meanR(ColorOpponentIds_conewts,2)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),atanh(meanR(Sconedominated_conewts,2)),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6); 
lims = [1 2 4 8 16];

for i = 1:numcells
    tmp = [b(i) b(i+numcells)]*[log10(lims); ones(size(lims))];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[log10(lims); ones(size(lims))],'g', 'Linewidth', 2); hold on;
axis square; xlabel('SNR'); ylabel('Z(R)'); title('Gabor'); set(gca,'XScale','log'); 
hold off;

subplot(122); plot(SNR(LumIds_conewts),atanh(meanR(LumIds_conewts,3)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SNR(ColorOpponentIds_conewts),atanh(meanR(ColorOpponentIds_conewts,3)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SNR(Sconedominated_conewts),atanh(meanR(Sconedominated_conewts,3)),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[1 16],'Ylim',[0 1.6],'XTick',[1 2 4 8 16],'YTick',0:0.8:1.6,'XScale','log'); 

X = [repmat(data1log,1,numcells).*celltype_dv celltype_dv];
[b, bint, r_DoG,~, stats1_DoG] = regress(data3,X);
[b_full, bint_full, r_full_DoG,~,stats0_DoG] = regress(data3,[data1log ones(size(data1log))]);
for i = 1:numcells
    tmp = [b(i) b(i+numcells)]*[log10(lims); ones(size(lims))];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[log10(lims); ones(size(lims))],'g', 'Linewidth', 2); hold on;
axis square; xlabel('SNR'); ylabel('Z(R)'); title('DoG'); hold off;
plot_counter = plot_counter + 1;

% Doing some stats: F test
dp = numcells*2-2;
RSSE_Gabor = sum(r_Gabor.^2);
RSSE_full_Gabor = sum(r_full_Gabor.^2);
Fstat_Gabor = ((RSSE_full_Gabor - RSSE_Gabor)/dp)/(RSSE_Gabor/(numel(data1)-numcells*2));
p_Gabor = 1-fcdf(Fstat_Gabor,dp,(numel(data1)-numcells*2));

RSSE_DoG = sum(r_DoG.^2);
RSSE_full_DoG = sum(r_full_DoG.^2);
Fstat_DoG = ((RSSE_full_DoG - RSSE_DoG)/dp)/(RSSE_DoG/(numel(data1)-numcells*2));
p_DoG = 1-fcdf(Fstat_DoG,dp,(numel(data1)-numcells*2));


%% Figure 6: Gabor phases and aspect ratios of Luminance, L-M and S cone dominated cells 
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

numcells = size(Output_List,1);
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
meanR_GaborDOG = meanR(:,2)-meanR(:,3); % Gabor - DoG

% storing gabor phases of the cells
load Gaborparams
gaborphases = zeros(1,numcells);
aspectratio = zeros(1,numcells);
for ii = 1:numcells
    gaborphases(ii) = rem(Gaborparams{ii}.phi,pi)*180/pi;
    aspectratio(ii) = Gaborparams{ii}.gamma;
end
Singleopponent = logical(Singleopponent);
gaborphases_relevantcells = gaborphases(~Singleopponent & simplecells);
aspectratio_relevantcells = aspectratio(~Singleopponent & simplecells);

bins1 = 0:10:90;
bins2 = logspace(log10(0.3),log10(10),10);
Lumid = zeros(size(gaborphases_relevantcells))'; Lumid(LumIds_conewts) = 1;
COid = zeros(size(gaborphases_relevantcells))'; COid(ColorOpponentIds_conewts) = 1;
Sid = zeros(size(gaborphases_relevantcells))'; Sid(Sconedominated_conewts) = 1;
DOid = zeros(size(gaborphases_relevantcells))'; DOid([ColorOpponentIds_conewts Sconedominated_conewts]) = 1;
Lumid = logical(Lumid); COid = logical(COid); Sid = logical(Sid); DOid = logical(DOid); 

% All cells and Best fitting Gabor cells
figure(plot_counter); set(gcf,'Name','Gabor Phases & aspect ratios');
subplot(321),histogram(90-abs(90-gaborphases_relevantcells(Lumid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 30],'YTick',0:10:30); axis square; title('Lum'); 
subplot(322),histogram(aspectratio_relevantcells(Lumid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid)),30,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Ylim',[0 30],'YTick',0:15:30,'Xlim',[0.3 10],'XTick',[0.3 1 3 10], 'Xscale','log'); axis square; title('Lum'); 
subplot(323),histogram(90-abs(90-gaborphases_relevantcells(COid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 30],'YTick',0:10:30); axis square; title('L-M'); 
subplot(324),histogram(aspectratio_relevantcells(COid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);  hold on; plot(median(aspectratio_relevantcells(COid)),30,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(COid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(COid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0.3 10],'XTick',[0.3 1 3 10],'Ylim',[0 30],'YTick',0:15:30,'Xscale','log'); axis square;title('L-M'); 
subplot(325),histogram(90-abs(90-gaborphases_relevantcells(Sid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 10],'YTick',0:5:10); axis square; title('S');
subplot(326),histogram(aspectratio_relevantcells(Sid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Sid)),15,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(Sid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Sid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0.3 10],'XTick',[0.3 1 3 10],'Ylim',[0 15],'YTick',0:5:15,'Xscale','log'); axis square; title('S');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% For cells that are better fit by Gabor 
% Mann-Whitney U test for phases 
p1 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)));
p2 = ranksum(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)),90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG>0)));
p3 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG>0)));

% Mann-Whitney U test for aspect ratios
p4 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),aspectratio_relevantcells(COid & meanR_GaborDOG>0));
p5 = ranksum(aspectratio_relevantcells(Sid & meanR_GaborDOG>0),aspectratio_relevantcells(COid & meanR_GaborDOG>0));
p6 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),aspectratio_relevantcells(Sid & meanR_GaborDOG>0));

data = [aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)';aspectratio_relevantcells(COid & meanR_GaborDOG>0)';aspectratio_relevantcells(Sid & meanR_GaborDOG>0)'];
data3 = [90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0))';90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0))';90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG>0))'];
group = [ones(size(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)')); 2*ones(size(aspectratio_relevantcells(COid & meanR_GaborDOG>0)')); 3*ones(size(aspectratio_relevantcells(Sid & meanR_GaborDOG>0)'))];
p = kruskalwallis(data,group,'off');
p25 = kruskalwallis(data3,group,'off');

% Hartigan's dip test for Gabor phases.
[dip1, p7,~,~]=HartigansDipSignifTest(sort(90-abs(90-gaborphases_relevantcells(Lumid))),1000);
[dip2, p8,~,~]=HartigansDipSignifTest(sort(90-abs(90-gaborphases_relevantcells(COid))),1000);
[dip3, p9,~,~]=HartigansDipSignifTest(sort(90-abs(90-gaborphases_relevantcells(Sid))),1000);

% For ALL cells 
% Mann-Whitney U test for phases 
p10 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid)),90-abs(90-gaborphases_relevantcells(COid)));
p11 = ranksum(90-abs(90-gaborphases_relevantcells(COid)),90-abs(90-gaborphases_relevantcells(Sid)));
p12 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid)),90-abs(90-gaborphases_relevantcells(Sid)));

% Mann-Whitney U test for aspect ratios
p13 = ranksum(aspectratio_relevantcells(Lumid),aspectratio_relevantcells(COid));
p14 = ranksum(aspectratio_relevantcells(Sid),aspectratio_relevantcells(COid));
p15 = ranksum(aspectratio_relevantcells(Lumid),aspectratio_relevantcells(Sid));
p21 = ranksum(aspectratio_relevantcells(Lumid),aspectratio_relevantcells(DOid));

data1 = [aspectratio_relevantcells(Lumid)';aspectratio_relevantcells(COid)';aspectratio_relevantcells(Sid)'];
data2 = [90-abs(90-gaborphases_relevantcells(Lumid))';90-abs(90-gaborphases_relevantcells(COid))';90-abs(90-gaborphases_relevantcells(Sid))'];
group = [ones(size(aspectratio_relevantcells(Lumid)')); 2*ones(size(aspectratio_relevantcells(COid)')); 3*ones(size(aspectratio_relevantcells(Sid)'))];
p16 = kruskalwallis(data1,group,'off');
p20 = kruskalwallis(data2,group,'off');

% Hartigan's dip test for Gabor phases.
[dip1, p17,~,~]=HartigansDipSignifTest(sort(90-abs(90-gaborphases_relevantcells(Lumid))),1000);
[dip2, p18,~,~]=HartigansDipSignifTest(sort(90-abs(90-gaborphases_relevantcells(COid))),1000);
[dip3, p19,~,~]=HartigansDipSignifTest(sort(90-abs(90-gaborphases_relevantcells(Sid))),1000);

%% Figure 7: Model comaprison, Crescent vs Gabor

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

bins2 = -10:1:10;
% Based on CV-SSE
figure(plot_counter); set(gcf,'Name','Model Comparison: Crescent vs Gabor');
subplot(311),plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');
subplot(312), plot(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');
subplot(313), plot(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('S');
plot_counter = plot_counter + 1;

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

idx = [Lumind(18); Lumind(2); COind(23); COind(18); Sconeind(22); Sconeind(11)];
figure(plot_counter);
L = numel(idx);
tmp = [-1 -1 -1 -1 1 1];
lambda = 12;
for ii = 1:numel(idx)
    jj = idx(ii);
    %STA
    im = Output_List{jj,2};
    normfactor = 0.5/(max(abs(im(:)))+0.01);
    im = normfactor*im + 0.5;
    im = reshape(im,[10 10 3]);
    subplot(2,L,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    % Crescent fit
    im = tmp(ii)*modelfits{jj,1};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(2,L,ii+L);image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
end
plot_counter = plot_counter + 1;

% Some more analyses for comparing the results between Lum, L-M and S cells: Gabor vs Crescent
% sign rank sum
p1 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1));
p2 = signrank(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,1));
p3 = signrank(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,1));

% Mann-Whitney U test
p4 = ranksum(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,1),meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,1));
p5 = ranksum(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,1),meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,1));
p6 = ranksum(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,1),meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,1));

%% Figure 8: Analysis of spatial opponency - non-paramteric method using spatial weighting function
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
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

load Deviation.mat
meanR = zeros(size(Deviation));
for ii = 1:size(Deviation,1)
    for jj = 1:size(Deviation,2)
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end

% Need to assess the spatial opponecny of model fits 
SPATIALRF = Output_List(:,4);
SOI_index_svd = []; % Spatial opponency index using SVD spatial RF 
for ii=1:size(SPATIALRF,1)
    % Based on SVD derived maps 
    tmp_svd = SPATIALRF{ii};
    P_svd = abs(sum(sum(tmp_svd(tmp_svd>0)))); N_svd = abs(sum(sum(tmp_svd(tmp_svd<0))));
    SOI_index_svd = [SOI_index_svd; 1-abs((P_svd-N_svd)/(P_svd+N_svd))];
    
end

examplecell_idx = [Lumind(18); Lumind(2); COind(23); COind(18); Sconeind(22); Sconeind(11)];

% Plotting the results on SOI index that is derived from SVD based spatial RF
GABORbetterthanDOG = meanR(:,2)>meanR(:,3);
bins = 0:0.1:1.0;
figure(plot_counter); set(gcf,'Name','Spatial opponency index: based on SVD spatial RF');
subplot(321); plot(SOI_index_svd(Lumind),meanR(Lumind,2)-meanR(Lumind,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SOI_index_svd(COind),meanR(COind,2)-meanR(COind,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(SOI_index_svd(Sconeind),meanR(Sconeind,2)-meanR(Sconeind,3),'o','MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.3 1],'Ylim',[-.5 0.5],'XTick',0.3:0.35:1.0); xlabel('SOI'); ylabel('R(Gabor)-R(DoG)');
subplot(322); histogram(SOI_index_svd(Lumind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(Lumind(GABORbetterthanDOG(Lumind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(Lumind)),30,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(Lumind(GABORbetterthanDOG(Lumind)))),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 50],'Xlim',[0.3 1],'XTick',0.3:0.35:1.0); xlabel('SOI'); ylabel('# cells'); title('LUM');
subplot(323); histogram(SOI_index_svd(COind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(COind(GABORbetterthanDOG(COind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(COind)),20,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(COind(GABORbetterthanDOG(COind)))),20,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30],'Xlim',[0.3 1],'XTick',0.3:0.35:1.0); xlabel('SOI'); ylabel('# cells'); title('DO');
subplot(324); histogram(SOI_index_svd(Sconeind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
histogram(SOI_index_svd(Sconeind(GABORbetterthanDOG(Sconeind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(Sconeind)),9,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(Sconeind(GABORbetterthanDOG(Sconeind)))),9,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30],'Xlim',[0.3 1],'XTick',0.3:0.35:1.0); xlabel('SOI'); ylabel('# cells'); title('S');
subplot(325); histogram(SOI_index_svd(Singleopponent),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
plot(median(SOI_index_svd(Singleopponent)),25,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30],'XTick',0:0.25:1.0); xlabel('SOI'); ylabel('# cells'); title('Single-opponent');
subplot(326); histogram(SOI_index_svd(Singleopponent),bins,'Displaystyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on; 
histogram(SOI_index_svd([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2);
plot(median(SOI_index_svd(Singleopponent)),30,'v','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
plot(median(SOI_index_svd([Lumind; COind; Sconeind])),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 80],'XTick',0:0.25:1.0); xlabel('SOI'); ylabel('# cells'); legend('SO','LUM & DO'); title('All cells');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Checking if there are any correlatiosn
idx = [Lumind; COind; Sconeind]';
[r1,p1] = corr(SOI_index_svd(idx),meanR(idx,2)-meanR(idx,3),'type','Spearman');
[r2,p2] = corr(SOI_index_svd(idx),SNR(idx),'type','Spearman');

% Checking whether the SOI is different among the cell types 
group = [ones(size(Lumind)); 2*ones(size(COind)); 3*ones(size(Sconeind))];
data1 = SOI_index_svd(idx);
data2 = meanR([Lumind; COind; Sconeind],2)- meanR([Lumind; COind; Sconeind],3);
p = kruskalwallis(data1,group,'off');
medidxs = data1>mean(data1);
p3 = kruskalwallis(data2(medidxs),group(medidxs),'off');

p4 = ranksum(SOI_index_svd(ind),SOI_index_svd(Singleopponent));

%% Figure 9: Eye movement analysis for population of cells

close all; 
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
Singleopponent = logical(Singleopponent);
ind = find(~Singleopponent & simplecells);

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
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% Pulling the eye movements
eyeAmp = Output_List(:,29);
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

% Calculating the spatial opponency index  
SPATIALRF = Output_List(:,4);
SOI_index_svd = []; % Spatial opponency index using SVD spatial RF 
for ii=1:size(SPATIALRF,1)
    % Based on SVD derived maps 
    tmp_svd = SPATIALRF{ii};
    P_svd = abs(sum(sum(tmp_svd(tmp_svd>0)))); N_svd = abs(sum(sum(tmp_svd(tmp_svd<0))));
    SOI_index_svd = [SOI_index_svd; 1-abs((P_svd-N_svd)/(P_svd+N_svd))];
    
end

% Plotting the population results
figure(plot_counter);
Amp = cell2mat(Output_List(:,29));
subplot(311); 
plot(Amp(Lumind), meanR(LumIds_conewts,2)-meanR(LumIds_conewts,3), 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Amp(COind), meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3), 'o', 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Amp(Sconeind), meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,3), 'o', 'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0 0.40], 'Ylim', [-0.5 0.5]); 
xlabel('median eye displacement (in degrees)'); ylabel('R(Gabor)-R(DoG)');

subplot(312); 
plot(Amp(Lumind), meanR(LumIds_conewts,2)-meanR(LumIds_conewts,1), 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Amp(COind), meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,1), 'o', 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Amp(Sconeind), meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,1), 'o', 'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0 0.40], 'Ylim', [-0.5 0.5]); 
xlabel('median eye displacement (in degrees)'); ylabel('R(Gabor)-R(modified DoG)');

% Plotting the relationship between eye movements and spatial opponency index 
subplot(313); 
plot(Amp(Lumind), SOI_index_svd(Lumind), 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Amp(COind), SOI_index_svd(COind), 'o', 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(Amp(Sconeind), SOI_index_svd(Sconeind), 'o', 'MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0.0 0.40], 'Ylim', [0.3 1.0]); 
xlabel('median eye displacement (in degrees)'); ylabel('Spatial Opponency Index');
plot_counter = plot_counter + 1;


% Doing some stats for eye movement amplitude vs. R(Gabor)-R(DoG)
combined_ind = [LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts];
[r1,p1] = corr(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,3), Amp(Lumind),'type','Spearman');
[r2,p2] = corr(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,3), Amp(COind),'type','Spearman','rows','complete');
[r3,p3] = corr(meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,3), Amp(Sconeind),'type','Spearman','rows','complete');
[r_combined,p_combined] = corr(meanR(combined_ind,2)-meanR(combined_ind,3), Amp([Lumind; COind; Sconeind]),'type','Spearman','rows','complete');

% Eye movement amplitude vs. R(Gabor)-R(modified DoG)
[r4,p4] = corr(meanR(LumIds_conewts,2)-meanR(LumIds_conewts,1), Amp(Lumind),'type','Spearman');
[r5,p5] = corr(meanR(ColorOpponentIds_conewts,2)-meanR(ColorOpponentIds_conewts,1), Amp(COind),'type','Spearman','rows','complete');
[r6,p6] = corr(meanR(Sconedominated_conewts,2)-meanR(Sconedominated_conewts,1), Amp(Sconeind),'type','Spearman','rows','complete');
[r_combined2, p_combined2] = corr(meanR(combined_ind,2)-meanR(combined_ind,1), Amp([Lumind; COind; Sconeind]),'type','Spearman','rows','complete');

% Eye movement amplitude vs. SOI index
[r7,p7] = corr(SOI_index_svd(Lumind), Amp(Lumind),'type','Spearman');
[r8,p8] = corr(SOI_index_svd(COind), Amp(COind),'type','Spearman','rows','complete');
[r9,p9] = corr(SOI_index_svd(Sconeind), Amp(Sconeind),'type','Spearman','rows','complete');
[r_combined3, p_combined3] = corr(SOI_index_svd([Lumind; COind; Sconeind]), Amp([Lumind; COind; Sconeind]),'type','Spearman','rows','complete');

%% Old addtional figure 1: Analysis of spatial opponency - paramteric method using Gabor fits 
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
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

load Deviation.mat
meanR = zeros(size(Deviation));
for ii = 1:size(Deviation,1)
    for jj = 1:size(Deviation,2)
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end

% Need to assess the spatial opponecny of model fits 
MODELFITS = modelfits;
SOI_index = []; % Spatial opponency index 

for ii=1:size(MODELFITS,1)
 
    % Based on Gabor model
    tmp = MODELFITS{ii,2}; % Gabor fits
    P = abs(sum(sum(tmp(tmp>0)))); N = abs(sum(sum(tmp(tmp<0))));
    SOI_index = [SOI_index; 1-abs((P-N)/(P+N))];
        
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
axis square; set(gca,'Tickdir','out','Ylim',[0 20]); xlabel('SOI'); ylabel('# cells'); title('DO');
subplot(324); histogram(SOI_index(Sconeind),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
histogram(SOI_index(Sconeind(GABORbetterthanDOG(Sconeind))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index(Sconeind)),9,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index(Sconeind(GABORbetterthanDOG(Sconeind)))),9,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 15]); xlabel('SOI'); ylabel('# cells'); title('S');
subplot(325); histogram(SOI_index(Singleopponent),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; 
plot(median(SOI_index(Singleopponent)),20,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
axis square; set(gca,'Tickdir','out','Ylim',[0 40]); xlabel('SOI'); ylabel('# cells'); title('Single-opponent (Ratio of power<1.2)');
subplot(326); histogram(SOI_index(Singleopponent),bins,'Displaystyle','stairs','EdgeColor',[1 0.5 0],'Linewidth',2); hold on; 
histogram(SOI_index([Lumind; COind; Sconeind]),bins,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2);
plot(median(SOI_index(Singleopponent)),30,'v','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
plot(median(SOI_index([Lumind; COind; Sconeind])),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('SOI'); ylabel('# cells'); legend('SO','LUM & DO'); title('All cells');
plot_counter = plot_counter + 1;

% Checking if there are any correlatiosn
idx = [Lumind; COind; Sconeind]';
[r1,p1] = corr(SOI_index(idx),meanR(idx,2)-meanR(idx,3),'type','Spearman');
[r2,p2] = corr(SOI_index(idx),SNR(idx),'type','Spearman');

% Checking whether the SOI is different among the cell types 
group = [ones(size(Lumind)); 2*ones(size(COind)); 3*ones(size(Sconeind))];
data1 = SOI_index(idx);
p = kruskalwallis(data1,group,'off');

p1 = ranksum(SOI_index(ind),SOI_index(Singleopponent));

%% Additional figure 1: Relationship between the RGB color space, cone-contrast space and cone-oppponent color (DKL) space
if ~exist('plot_counter')
    plot_counter = 1;
end
stro=nex2stro(findfile('K033108004'));
gammaTable = stro.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);
% Reconstructing the M matrix
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
% % Getting the background rgb/lms
ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*[.5; .5; .5];
nstim = 10000;
rgb = zeros(nstim, 3);
for i = 1:nstim
    rgb(i,:) = normrnd(bkgndrgb, [.15 .15 .15]')';
end
lms = M*rgb';
cc = (lms'-repmat(bkgndlms',nstim,1))./repmat(bkgndlms',nstim,1);
LvM = cc(:,1)-cc(:,2);
LM = cc(:,1)+cc(:,2);
S = cc(:,3);
figure(plot_counter);
subplot(331); plot(rgb(:,1)-bkgndrgb(1),rgb(:,2)-bkgndrgb(2),'k.');
xlabel('R'); ylabel('G'); axis square; set(gca,'xlim',[-1 1],'ylim',[-1 1],'Tickdir','out');
subplot(332); plot(rgb(:,2)-bkgndrgb(2),rgb(:,3)-bkgndrgb(3),'k.');
xlabel('G'); ylabel('B'); axis square; set(gca,'xlim',[-1 1],'ylim',[-1 1],'Tickdir','out');
subplot(333); plot(rgb(:,1)-bkgndrgb(1),rgb(:,3)-bkgndrgb(3),'k.');
xlabel('R'); ylabel('B'); axis square; set(gca,'xlim',[-1 1],'ylim',[-1 1],'Tickdir','out');
subplot(334); plot(LvM,S,'k.');
xlabel('L-M'); ylabel('S'); axis square; set(gca,'xlim',[-2 2],'ylim',[-2 2],'Tickdir','out');
subplot(335); plot(LvM,LM,'k.');
xlabel('L-M'); ylabel('L+M'); axis square; set(gca,'xlim',[-2 2],'ylim',[-2 2],'Tickdir','out');
subplot(336); plot(S,LM,'k.');
xlabel('S'); ylabel('L+M'); axis square; set(gca,'xlim',[-2 2],'ylim',[-2 2],'Tickdir','out');
subplot(337); plot(cc(:,1),cc(:,2),'k.');
xlabel('L'); ylabel('M'); axis square; set(gca,'xlim',[-2 2],'ylim',[-2 2],'Tickdir','out');
subplot(338); plot(cc(:,1),cc(:,3),'k.');
xlabel('L'); ylabel('S'); axis square; set(gca,'xlim',[-2 2],'ylim',[-2 2],'Tickdir','out');
subplot(339); plot(cc(:,2),cc(:,3),'k.');
xlabel('M'); ylabel('S'); set(gca,'xlim',[-2 2],'ylim',[-2 2],'Tickdir','out'); axis square;
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

%% Additional figure 2: Elliptical DOG (gamma) model fits for example cells (could be a reviewer figure)
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
ind = find(~Singleopponent & simplecells);

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
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% idx = [Lumind(8); COind(30); Sconeind(20)];
idx = [Lumind(19); Lumind(2); COind(23); COind(18); Sconeind(22); Sconeind(11)];
figure(plot_counter);
L = numel(idx);
swt_tmp = [-1 -1 1 1 -1 1];
gf_tmp = [-1 1 1 -1 -1 1];
dogf_tmp = [-1 1 1 1 -1 -1];
lambda = 12;
for ii = 1:numel(idx)
    jj = idx(ii);
    %STA
    im = Output_List{jj,2};
    normfactor = 0.5/(max(abs(im(:)))+0.01);
    im = normfactor*im + 0.5;
    im = reshape(im,[10 10 3]);
    subplot(5,L,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    
      
    % Spatial weighting function
    im = swt_tmp(ii)*Output_List{jj,4};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+L); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    
    % Gabor fit
    im = gf_tmp(ii)*modelfits{jj,2};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+2*L);image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    
    % DoG fit
    im = -1*dogf_tmp(ii)*modelfits{jj,3};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+3*L); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));


     % DoG gamma fit
    im = -1*dogf_tmp(ii)*modelfits{jj,6};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+4*L); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));

end
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

%% Additional figure 3: Scatterplot showing the DoG and the crescent shaped model (reviewer's eyes only)

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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
SNR = Zmax(Z_cellsofinterest)./Z(Z_cellsofinterest,1);
SNR = SNR(~Singleopponent & simplecells);

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

idx = [18; 2; 23; 18; 22; 11];
bins2 = -10:1:10;
% Based on CV-R
figure(plot_counter); set(gcf,'Name','Model Comparison: Crescent vs. DoG');
subplot(131);  hold on;
for ii = 1:length(LumIds_conewts)   
    h(ii) = plot(meanR(LumIds_conewts(ii),1),meanR(LumIds_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);  
end
line([0 1],[0 1],'Color','k');  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');

subplot(132), hold on;
for ii = 1:length(ColorOpponentIds_conewts)
    h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),1),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
line([0 1],[0 1],'Color','k');  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

subplot(133), hold on;
for ii = 1:length(Sconedominated_conewts)
    h(ii) = plot(meanR(Sconedominated_conewts(ii),1),meanR(Sconedominated_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
    set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
end
line([0 1],[0 1],'Color','k');  xlabel('Crescent'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('S');
plot_counter = plot_counter + 1;

% Wilcoxon signed rank test between non-concentric DoG and DoG
p1 = signrank(meanR(LumIds_conewts,1),meanR(LumIds_conewts,3));
p2 = signrank(meanR(ColorOpponentIds_conewts,1),meanR(ColorOpponentIds_conewts,3));
p3 = signrank(meanR(Sconedominated_conewts,1),meanR(Sconedominated_conewts,3));

%% Additional figure 4: Calculating a new RF area suggested by Greg (pi*sigma.^2/gamma: area of 1 STD ellipse)

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
load Ratio_of_power.mat
load Gaborparams.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);
ind = find(~Singleopponent & simplecells);

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
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% Extracting the RF area from the Gabor fit and the eccentricity of cells
ecc = [];
RFarea_gaborfit = [];
for ii = 1:numel(Gaborparams)
    %SNR from the Gabor fit
    RFarea_gaborfit = [RFarea_gaborfit; 0.04*pi*Gaborparams{ii}.sigma.^2/(Gaborparams{ii}.gamma)];
    
    % RF eccentricity is stored in column 8
    tmpecc = cell2mat(Output_List(ii,8));
    ecc = [ecc; sqrt(sum(tmpecc.^2,2))/10];
end

figure(plot_counter);
subplot(131); plot(ecc(Lumind),RFarea_gaborfit(Lumind),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; set(gca,'Tickdir','out', 'Xlim', [0 10],'YScale','log','Ylim',[0.005 10]);  ylabel('RF area'); 
title('LUM'); xlabel('Eccentricity in degrees');
subplot(132); plot(ecc(COind),RFarea_gaborfit(COind),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out', 'Xlim', [0 10],'YScale','log','Ylim',[0.005 10]);  ylabel('RF area'); 
title('L-M'); xlabel('Eccentricity in degrees');
subplot(133); plot(ecc(Sconeind),RFarea_gaborfit(Sconeind),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out', 'Xlim', [0 10],'YScale','log','Ylim',[0.005 10]);  ylabel('RF area'); 
title('S-cone sensitive'); xlabel('Eccentricity in degrees');
plot_counter = plot_counter + 1;

% Some stats
[r1,p1] = corr(ecc(Lumind),log(RFarea_gaborfit(Lumind)),'type','Spearman');
[r2,p2] = corr(ecc(COind),log(RFarea_gaborfit(COind)),'type','Spearman');
[r3,p3] = corr(ecc(Sconeind),log(RFarea_gaborfit(Sconeind)),'type','Spearman');


%% Additional figure 5: Control Plot of cone weights of all the screened cells
% Calculating conewts of the cells of interest: SVD 
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
LumIds_conewts = find(conewts_svd(2,:) + conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)-0.5).^2)<0.3);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & conewts_svd(1,:)<-0.1 & conewts_svd(2,:)>0.1);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

figure(plot_counter); set(gcf,'Name','Cone wts') 
plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

% Further classification of S-cone sensitive cells
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Scone_nonopponent = Sconesensitive(:,sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1);
SvsLM = Sconesensitive(:,sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==-1);
SMvsL = Sconesensitive(:,sign(Sconesensitive(1,:))==-1 & sign(Sconesensitive(3,:))==1);
SLvsM = Sconesensitive(:,sign(Sconesensitive(1,:))==-1 & sign(Sconesensitive(3,:))==-1);

%% Additional figure 6: Model comparison, Gabor vs DOG, Gabor phases, aspect ratios, non-concentric DoG vs. Gabor  
if ~exist('plot_counter')
    plot_counter = 1;
end

% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
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
LumIds_conewts = find(conewts_svd(2,:) + conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)-0.5).^2)<0.3);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & conewts_svd(1,:)<-0.1 & conewts_svd(2,:)>0.1);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

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
meanR_GaborDOG = meanR(:,2)-meanR(:,3); % Gabor - DoG

% Need to assess the spatial opponecny of model fits 
SPATIALRF = Output_List(:,4);
SOI_index_svd = []; % Spatial opponency index using SVD spatial RF 
for ii=1:size(SPATIALRF,1)
    % Based on SVD derived maps 
    tmp_svd = SPATIALRF{ii};
    P_svd = abs(sum(sum(tmp_svd(tmp_svd>0)))); N_svd = abs(sum(sum(tmp_svd(tmp_svd<0))));
    SOI_index_svd = [SOI_index_svd; 1-abs((P_svd-N_svd)/(P_svd+N_svd))];
    
end
SOI_index_svd = SOI_index_svd(~Singleopponent & simplecells);

% Plotting the results on SOI index that is derived from SVD based spatial RF
GABORbetterthanDOG = meanR(:,2)>meanR(:,3);

% storing gabor phases of the cells
load Gaborparams
gaborphases = zeros(1,numcells);
aspectratio = zeros(1,numcells);
for ii = 1:numcells
    gaborphases(ii) = rem(Gaborparams{ii}.phi,pi)*180/pi;
    aspectratio(ii) = Gaborparams{ii}.gamma;
end
Singleopponent = logical(Singleopponent);
gaborphases_relevantcells = gaborphases(~Singleopponent & simplecells);
aspectratio_relevantcells = aspectratio(~Singleopponent & simplecells);

bins1 = 0:10:90;
bins2 = logspace(log10(0.3),log10(10),10);
Lumid = zeros(size(gaborphases_relevantcells))'; Lumid(LumIds_conewts) = 1;
COid = zeros(size(gaborphases_relevantcells))'; COid(ColorOpponentIds_conewts) = 1;
Sid = zeros(size(gaborphases_relevantcells))'; Sid(Sconedominated_conewts) = 1;
Lumid = logical(Lumid); COid = logical(COid); Sid = logical(Sid);

% Plotting
figure(plot_counter); set(gcf,'Name','Model Comparison: DoG vs Gabor');
subplot(521);  plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum'); hold off;

subplot(522), 
plot(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M'); hold off;

subplot(523),histogram(90-abs(90-gaborphases_relevantcells(Lumid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 20],'YTick',0:10:20); axis square; title('Lum'); 

subplot(524),histogram(90-abs(90-gaborphases_relevantcells(COid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 40],'YTick',0:20:40); axis square; title('L-M'); 

subplot(525),histogram(aspectratio_relevantcells(Lumid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid)),20,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Ylim',[0 20],'YTick',0:10:20,'Xlim',[0.3 10],'XTick',[0.3 1 3 10], 'Xscale','log'); axis square; title('Lum'); 

subplot(526),histogram(aspectratio_relevantcells(COid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);  hold on; plot(median(aspectratio_relevantcells(COid)),40,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(COid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(COid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0.3 10],'XTick',[0.3 1 3 10],'Ylim',[0 40],'YTick',0:20:40,'Xscale','log'); axis square;title('L-M'); 

subplot(527),plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');

subplot(528), plot(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

bins = 0:0.1:1.0;
subplot(529); histogram(SOI_index_svd(LumIds_conewts),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(LumIds_conewts(GABORbetterthanDOG(LumIds_conewts))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(LumIds_conewts)),35,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(LumIds_conewts(GABORbetterthanDOG(LumIds_conewts)))),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 40]); xlabel('SOI'); ylabel('# cells'); title('LUM');

subplot(5,2,10); histogram(SOI_index_svd(ColorOpponentIds_conewts),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(ColorOpponentIds_conewts(GABORbetterthanDOG(ColorOpponentIds_conewts))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(ColorOpponentIds_conewts)),25,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(ColorOpponentIds_conewts(GABORbetterthanDOG(ColorOpponentIds_conewts)))),20,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 30]); xlabel('SOI'); ylabel('# cells'); title('DO');
plot_counter = plot_counter + 1;

% Some more analyses for comparing the results between Lum, L-M and S cells: Gabor vs DOG
% sign rank sum
p1 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3));
p2 = signrank(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,3));

% Mann-Whitney U test for phases 
p3 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)));
p4 = ranksum(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)),90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG>0)));

% Mann-Whitney U test for aspect ratios
p5 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),aspectratio_relevantcells(COid & meanR_GaborDOG>0));
p6 = ranksum(aspectratio_relevantcells(Sid & meanR_GaborDOG>0),aspectratio_relevantcells(COid & meanR_GaborDOG>0));

% Comparing non-concentric DoG and Gabor 
p7 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1));
p8 = signrank(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,1));


%% Additional figure 7: Cone weight plot with a relaxed criteria 
% Trying a new criteria, simple cells -> cone non-opponent, DO cells -> cone opponent 
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
Singleopponent = logical(Singleopponent);
% Singleopponent = logical(Ratio_of_power<1.2);

numcells = size(Output_List,1);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
SNR = Zmax(Z_cellsofinterest)./Z(Z_cellsofinterest,1);
SNR = SNR(~Singleopponent & simplecells);

% calculating the cone weights
conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

LumIds_conewts = find(conewts_svd(1,:)>0 & conewts_svd(3,:)>0);
DOIds_conewts = find(conewts_svd(1,:)<0 | conewts_svd(3,:)<0);

figure(plot_counter); set(gcf,'Name','Cone wts') 
plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,DOIds_conewts),conewts_svd(2,DOIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

%% Additional figure 8: Model comparison, Gabor vs DOG, Gabor phases, aspect ratios, non-concentric DoG vs. Gabor with the relaxed criteria 
if ~exist('plot_counter')
    plot_counter = 1;
end

% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
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

LumIds_conewts = find(conewts_svd(1,:)>0 & conewts_svd(3,:)>0);
DOIds_conewts = find(conewts_svd(1,:)<0 | conewts_svd(3,:)<0);

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

% Need to assess the spatial opponecny of model fits 
SPATIALRF = Output_List(:,4);
SOI_index_svd = []; % Spatial opponency index using SVD spatial RF 
for ii=1:size(SPATIALRF,1)
    % Based on SVD derived maps 
    tmp_svd = SPATIALRF{ii};
    P_svd = abs(sum(sum(tmp_svd(tmp_svd>0)))); N_svd = abs(sum(sum(tmp_svd(tmp_svd<0))));
    SOI_index_svd = [SOI_index_svd; 1-abs((P_svd-N_svd)/(P_svd+N_svd))];
    
end
SOI_index_svd = SOI_index_svd(~Singleopponent & simplecells);

% Plotting the results on SOI index that is derived from SVD based spatial RF
GABORbetterthanDOG = meanR(:,2)>meanR(:,3);
meanR_GaborDOG = meanR(:,2)-meanR(:,3); % Gabor - DoG

% storing gabor phases of the cells
load Gaborparams
gaborphases = zeros(1,numcells);
aspectratio = zeros(1,numcells);
for ii = 1:numcells
    gaborphases(ii) = rem(Gaborparams{ii}.phi,pi)*180/pi;
    aspectratio(ii) = Gaborparams{ii}.gamma;
end
Singleopponent = logical(Singleopponent);
gaborphases_relevantcells = gaborphases(~Singleopponent & simplecells);
aspectratio_relevantcells = aspectratio(~Singleopponent & simplecells);

bins1 = 0:10:90;
bins2 = logspace(log10(0.3),log10(10),10);
Lumid = zeros(size(gaborphases_relevantcells))'; Lumid(LumIds_conewts) = 1;
DOid = zeros(size(gaborphases_relevantcells))'; DOid(DOIds_conewts) = 1;

Lumid = logical(Lumid); DOid = logical(DOid); 

% Plotting
figure(plot_counter); set(gcf,'Name','Model Comparison: DoG vs Gabor');
subplot(521);  plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum'); hold off;

subplot(522), 
plot(meanR(DOIds_conewts,2),meanR(DOIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M'); hold off;

subplot(523),histogram(90-abs(90-gaborphases_relevantcells(Lumid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 20],'YTick',0:10:20); axis square; title('Lum'); 

subplot(524),histogram(90-abs(90-gaborphases_relevantcells(DOid)),bins1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(DOid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 60],'YTick',0:30:60); axis square; title('L-M'); 

subplot(525),histogram(aspectratio_relevantcells(Lumid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid)),20,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Ylim',[0 20],'YTick',0:10:20,'Xlim',[0.3 10],'XTick',[0.3 1 3 10], 'Xscale','log'); axis square; title('Lum'); 

subplot(526),histogram(aspectratio_relevantcells(DOid),bins2,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);  hold on; plot(median(aspectratio_relevantcells(DOid)),50,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(DOid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(DOid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0.3 10],'XTick',[0.3 1 3 10],'Ylim',[0 60],'YTick',0:30:60,'Xscale','log'); axis square;title('L-M'); 

subplot(527),plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');

subplot(528), plot(meanR(DOIds_conewts,2),meanR(DOIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

bins = 0:0.1:1.0;
subplot(529); histogram(SOI_index_svd(LumIds_conewts),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(LumIds_conewts(GABORbetterthanDOG(LumIds_conewts))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(LumIds_conewts)),35,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(LumIds_conewts(GABORbetterthanDOG(LumIds_conewts)))),30,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 50],'YTick',0:25:50); xlabel('SOI'); ylabel('# cells'); title('LUM');

subplot(5,2,10); histogram(SOI_index_svd(DOIds_conewts),bins,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]); hold on;
histogram(SOI_index_svd(DOIds_conewts(GABORbetterthanDOG(DOIds_conewts))),bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
plot(median(SOI_index_svd(DOIds_conewts)),25,'v','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
plot(median(SOI_index_svd(DOIds_conewts(GABORbetterthanDOG(DOIds_conewts)))),20,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 60],'YTick',0:30:60); xlabel('SOI'); ylabel('# cells'); title('DO');

set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

% Some more analyses for comparing the results between Lum, L-M and S cells: Gabor vs DOG
% sign rank sum
p1 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3));
p2 = signrank(meanR(DOIds_conewts,2),meanR(DOIds_conewts,3));

% Mann-Whitney U test for phases 
p3 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),90-abs(90-gaborphases_relevantcells(DOid & meanR_GaborDOG>0)));

% Mann-Whitney U test for aspect ratios
p4 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),aspectratio_relevantcells(DOid & meanR_GaborDOG>0));

% Comparing non-concentric DoG and Gabor 
p5 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1));
p6 = signrank(meanR(DOIds_conewts,2),meanR(DOIds_conewts,1));

%% Additional figure 9: Model comparison of DOG vs Gabor (CV per accuracy, mean SSE & BIC

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
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

load Peraccuracy.mat
load SSE.mat
meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));
for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);

% BIC
N = 100;
numcrescentparams = 8;
numgaborparams = 8;
numDOGparams = 6;
load Crescenterror.mat
load Gaborerror.mat
load DOGerror.mat
CrescentBIC = N*log(Crescenterror(~Singleopponent & simplecells,:)/N) + numcrescentparams*log(N);
GaborBIC = N*log(Gaborerror(~Singleopponent & simplecells,:)/N) + numgaborparams*log(N);
DOGBIC = N*log(DOGerror(~Singleopponent & simplecells,:)/N) + numDOGparams*log(N);
BICs = [CrescentBIC GaborBIC DOGBIC];

bins2 = -10:1:10;
% Based CV peraccuracy, CV mean SSE, BICs 

figure(plot_counter); set(gcf,'Name','Model Comparison: DoG vs Gabor: Lum, L-M & S');
subplot(331),plot(meanperaccuracy(LumIds_conewts,2),meanperaccuracy(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(334),plot(meanSSE(LumIds_conewts,2),meanSSE(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 2],[0 2]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 2],'Ylim',[0 2],'XTick',0:1.0:2.0,'YTick',0:1.0:2.0); title('CV SSE');
subplot(337),plot(BICs(LumIds_conewts,2),BICs(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1300 -500],[-1300 -500]); xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1300 -500],'Xlim',[-1300 -500],'XTick',-1300:400:-500,'YTick',-1300:400:-500); title('BIC');
subplot(332),plot(meanperaccuracy(ColorOpponentIds_conewts,2),meanperaccuracy(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(335),plot(meanSSE(ColorOpponentIds_conewts,2),meanSSE(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 2],[0 2]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 2],'Ylim',[0 2],'XTick',0:1.0:2.0,'YTick',0:1.0:2.0); title('CV SSE');
subplot(338),plot(BICs(ColorOpponentIds_conewts,2),BICs(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1300 -500],[-1300 -500]); xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1300 -500],'Xlim',[-1300 -500],'XTick',-1300:400:-500,'YTick',-1300:400:-500); title('BIC');
subplot(333),plot(meanperaccuracy(Sconedominated_conewts,2),meanperaccuracy(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(336),plot(meanSSE(Sconedominated_conewts,2),meanSSE(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 2],[0 2]);  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 2],'Ylim',[0 2],'XTick',0:1.0:2.0,'YTick',0:1.0:2.0); title('CV SSE');
subplot(339),plot(BICs(Sconedominated_conewts,2),BICs(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1300 -500],[-1300 -500]); xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'TickDir','out','Ylim',[-1300 -500],'Xlim',[-1300 -500],'XTick',-1300:400:-500,'YTick',-1300:400:-500); title('BIC');
plot_counter = plot_counter + 1;

% Some more analyses for comparing the results between Lum, L-M and S cells: Gabor vs DOG
% sign rank sum: Mean per accuracy
p1 = signrank(meanperaccuracy(LumIds_conewts,2),meanperaccuracy(LumIds_conewts,3));
p2 = signrank(meanperaccuracy(ColorOpponentIds_conewts,2),meanperaccuracy(ColorOpponentIds_conewts,3));
p3 = signrank(meanperaccuracy(Sconedominated_conewts,2),meanperaccuracy(Sconedominated_conewts,3));

% Mann-Whitney U test
p4 = ranksum(meanperaccuracy(LumIds_conewts,2)-meanperaccuracy(LumIds_conewts,3),meanperaccuracy(ColorOpponentIds_conewts,2)-meanperaccuracy(ColorOpponentIds_conewts,3));
p5 = ranksum(meanperaccuracy(LumIds_conewts,2)-meanperaccuracy(LumIds_conewts,3),meanperaccuracy(Sconedominated_conewts,2)-meanperaccuracy(Sconedominated_conewts,3));
p6 = ranksum(meanperaccuracy(ColorOpponentIds_conewts,2)-meanperaccuracy(ColorOpponentIds_conewts,3),meanperaccuracy(Sconedominated_conewts,2)-meanperaccuracy(Sconedominated_conewts,3));

% sign rank sum: Mean SSE
p7 = signrank(meanSSE(LumIds_conewts,2),meanSSE(LumIds_conewts,3));
p8 = signrank(meanSSE(ColorOpponentIds_conewts,2),meanSSE(ColorOpponentIds_conewts,3));
p9 = signrank(meanSSE(Sconedominated_conewts,2),meanSSE(Sconedominated_conewts,3));

% Mann-Whitney U test
p10 = ranksum(meanSSE(LumIds_conewts,2)-meanSSE(LumIds_conewts,3),meanSSE(ColorOpponentIds_conewts,2)-meanSSE(ColorOpponentIds_conewts,3));
p11 = ranksum(meanSSE(LumIds_conewts,2)-meanSSE(LumIds_conewts,3),meanSSE(Sconedominated_conewts,2)-meanSSE(Sconedominated_conewts,3));
p12 = ranksum(meanSSE(ColorOpponentIds_conewts,2)-meanSSE(ColorOpponentIds_conewts,3),meanSSE(Sconedominated_conewts,2)-meanSSE(Sconedominated_conewts,3));

% sign rank sum: BIC
p13 = signrank(BICs(LumIds_conewts,2),BICs(LumIds_conewts,3));
p14 = signrank(BICs(ColorOpponentIds_conewts,2),BICs(ColorOpponentIds_conewts,3));
p15 = signrank(BICs(Sconedominated_conewts,2),BICs(Sconedominated_conewts,3));

% Mann-Whitney U test
p16 = ranksum(BICs(LumIds_conewts,2)-BICs(LumIds_conewts,3),BICs(ColorOpponentIds_conewts,2)-BICs(ColorOpponentIds_conewts,3));
p17 = ranksum(BICs(LumIds_conewts,2)-BICs(LumIds_conewts,3),BICs(Sconedominated_conewts,2)-BICs(Sconedominated_conewts,3));
p18 = ranksum(BICs(ColorOpponentIds_conewts,2)-BICs(ColorOpponentIds_conewts,3),BICs(Sconedominated_conewts,2)-BICs(Sconedominated_conewts,3));


%% Additional figure 10: Model comparison of Crescent vs Gabor (CV per accuracy, mean SSE & BIC)

if ~exist('plot_counter')
    plot_counter = 1;
end
% load Output_ListWN.mat
load Output_ListWN2.mat
load Singleopponent.mat
load modelfits.mat
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);

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

load Peraccuracy.mat
load SSE.mat

meanperaccuracy = zeros(size(Peraccuracy));
meanSSE = zeros(size(SSE));

for ii = 1:size(Peraccuracy,1)
    for jj = 1:size(Peraccuracy,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy{ii,jj});
        meanSSE(ii,jj) = mean(SSE{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);

% BIC
N = 100;
numcrescentparams = 8;
numgaborparams = 8;
numDOGparams = 6;
load Crescenterror.mat
load Gaborerror.mat
load DOGerror.mat
CrescentBIC = N*log(Crescenterror(~Singleopponent & simplecells,:)/N) + numcrescentparams*log(N);
GaborBIC = N*log(Gaborerror(~Singleopponent & simplecells,:)/N) + numgaborparams*log(N);
DOGBIC = N*log(DOGerror(~Singleopponent & simplecells,:)/N) + numDOGparams*log(N);
BICs = [CrescentBIC GaborBIC DOGBIC];

bins2 = -10:1:10;
% Based CV peraccuracy, CV mean SSE, BICs 

figure(plot_counter); set(gcf,'Name','Model Comparison: Crescent vs Gabor: Lum, L-M & S');
subplot(331),plot(meanperaccuracy(LumIds_conewts,2),meanperaccuracy(LumIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(334),plot(meanSSE(LumIds_conewts,2),meanSSE(LumIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 2],[0 2]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 2],'Ylim',[0 2],'XTick',0:1.0:2.0,'YTick',0:1.0:2.0); title('CV SSE');
subplot(337),plot(BICs(LumIds_conewts,2),BICs(LumIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1300 -500],[-1300 -500]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'TickDir','out','Ylim',[-1300 -500],'Xlim',[-1300 -500],'XTick',-1300:400:-500,'YTick',-1300:400:-500); title('BIC');
subplot(332),plot(meanperaccuracy(ColorOpponentIds_conewts,2),meanperaccuracy(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(335),plot(meanSSE(ColorOpponentIds_conewts,2),meanSSE(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0 2],[0 2]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 2],'Ylim',[0 2],'XTick',0:1.0:2.0,'YTick',0:1.0:2.0); title('CV SSE');
subplot(338),plot(BICs(ColorOpponentIds_conewts,2),BICs(ColorOpponentIds_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1300 -500],[-1300 -500]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'TickDir','out','Ylim',[-1300 -500],'Xlim',[-1300 -500],'XTick',-1300:400:-500,'YTick',-1300:400:-500); title('BIC');
subplot(333),plot(meanperaccuracy(Sconedominated_conewts,2),meanperaccuracy(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on; 
line([0.5 0.9],[0.5 0.9]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0.5 0.9],'Ylim',[0.5 0.9],'XTick',0.5:0.1:0.9,'YTick',0.5:0.1:0.9); title('CV per accuracy');
subplot(336),plot(meanSSE(Sconedominated_conewts,2),meanSSE(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([0 2],[0 2]);  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 2],'Ylim',[0 2],'XTick',0:1.0:2.0,'YTick',0:1.0:2.0); title('CV SSE');
subplot(339),plot(BICs(Sconedominated_conewts,2),BICs(Sconedominated_conewts,1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
line([-1300 -500],[-1300 -500]); xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'TickDir','out','Ylim',[-1300 -500],'Xlim',[-1300 -500],'XTick',-1300:400:-500,'YTick',-1300:400:-500); title('BIC');
plot_counter = plot_counter + 1;

% Some more analyses for comparing the results between Lum, L-M and S cells: Gabor vs Crescent
% sign rank sum: Mean per accuracy
p1 = signrank(meanperaccuracy(LumIds_conewts,2),meanperaccuracy(LumIds_conewts,1));
p2 = signrank(meanperaccuracy(ColorOpponentIds_conewts,2),meanperaccuracy(ColorOpponentIds_conewts,1));
p3 = signrank(meanperaccuracy(Sconedominated_conewts,2),meanperaccuracy(Sconedominated_conewts,1));

% Mann-Whitney U test
p4 = ranksum(meanperaccuracy(LumIds_conewts,2)-meanperaccuracy(LumIds_conewts,1),meanperaccuracy(ColorOpponentIds_conewts,2)-meanperaccuracy(ColorOpponentIds_conewts,1));
p5 = ranksum(meanperaccuracy(LumIds_conewts,2)-meanperaccuracy(LumIds_conewts,1),meanperaccuracy(Sconedominated_conewts,2)-meanperaccuracy(Sconedominated_conewts,1));
p6 = ranksum(meanperaccuracy(ColorOpponentIds_conewts,2)-meanperaccuracy(ColorOpponentIds_conewts,1),meanperaccuracy(Sconedominated_conewts,2)-meanperaccuracy(Sconedominated_conewts,1));

% sign rank sum: Mean SSE
p7 = signrank(meanSSE(LumIds_conewts,2),meanSSE(LumIds_conewts,1));
p8 = signrank(meanSSE(ColorOpponentIds_conewts,2),meanSSE(ColorOpponentIds_conewts,1));
p9 = signrank(meanSSE(Sconedominated_conewts,2),meanSSE(Sconedominated_conewts,1));

% Mann-Whitney U test
p10 = ranksum(meanSSE(LumIds_conewts,2)-meanSSE(LumIds_conewts,1),meanSSE(ColorOpponentIds_conewts,2)-meanSSE(ColorOpponentIds_conewts,1));
p11 = ranksum(meanSSE(LumIds_conewts,2)-meanSSE(LumIds_conewts,1),meanSSE(Sconedominated_conewts,2)-meanSSE(Sconedominated_conewts,1));
p12 = ranksum(meanSSE(ColorOpponentIds_conewts,2)-meanSSE(ColorOpponentIds_conewts,1),meanSSE(Sconedominated_conewts,2)-meanSSE(Sconedominated_conewts,1));

% sign rank sum: BIC
p13 = signrank(BICs(LumIds_conewts,2),BICs(LumIds_conewts,1));
p14 = signrank(BICs(ColorOpponentIds_conewts,2),BICs(ColorOpponentIds_conewts,1));
p15 = signrank(BICs(Sconedominated_conewts,2),BICs(Sconedominated_conewts,1));

% Mann-Whitney U test
p16 = ranksum(BICs(LumIds_conewts,2)-BICs(LumIds_conewts,1),BICs(ColorOpponentIds_conewts,2)-BICs(ColorOpponentIds_conewts,1));
p17 = ranksum(BICs(LumIds_conewts,2)-BICs(LumIds_conewts,1),BICs(Sconedominated_conewts,2)-BICs(Sconedominated_conewts,1));
p18 = ranksum(BICs(ColorOpponentIds_conewts,2)-BICs(ColorOpponentIds_conewts,1),BICs(Sconedominated_conewts,2)-BICs(Sconedominated_conewts,1));


%% STAs of L+M, L-M & S-cone sensitive cells
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
Lumind = ind(LumIds_conewts); Lumsubplot = ceil(sqrt(numel(Lumind)));
COind = ind(ColorOpponentIds_conewts); COsubplot = ceil(sqrt(numel(COind)));
Sconeind = ind(Sconedominated_conewts); Sconesubplot = ceil(sqrt(numel(Sconeind)));

% plotting STAs of lumninance simple cells
figure(plot_counter); set(gcf,'Name','Luminance cells');
for ii = 1:numel(Lumind)
    tmp_vec_gun = Output_List{Lumind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(10,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of L-M simple cells
figure(plot_counter); set(gcf,'Name','L-M: SVD');
for ii = 1:numel(COind)
    tmp_vec_gun = Output_List{COind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(10,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

% plotting STAs of S cone-dominated simple cells
figure(plot_counter); set(gcf,'Name','S cone: SVD');
for ii = 1:numel(Sconeind)
    tmp_vec_gun = Output_List{Sconeind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(10,10,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;

%% Single opponent cells, Cone weights and STAs

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

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;


% calculating the cone weights
conewts_svdSO = cell2mat(Output_List(Singleopponent & simplecells,23)');
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

figure(plot_counter); set(gcf,'Name','Cone wts:SO') 
subplot(211); plot(conewts_svdSO(1,:),conewts_svdSO(2,:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('L'), ylabel('M');
subplot(212); plot(conewts_svdSO(3,:),conewts_svdSO(2,:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out');  plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k');  xlabel('S'), ylabel('M');
plot_counter = plot_counter + 1;

SOind = find(Singleopponent) ;
% plotting STAs of Single-opponent cells
figure(plot_counter); set(gcf,'Name','SO: SVD');
for ii = 1:numel(SOind)
    tmp_vec_gun = Output_List{SOind(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
 
    subplot(8,8,ii); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
end
plot_counter = plot_counter + 1;