% A new script for preparing figure for the DO spatial structure paper
% Author - Abhishek De, 5/19
% Thesemdata were submitted to the JoV first submission

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

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
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

figure(plot_counter); set(gcf,'Name','Stim->STA');
for ii = 1:3
    im = rand(10,10,3);
    avg = mean(mean(im));
    im = im - repmat(avg,[10 10 1]);
    subplot(2,4,ii),image(im./(6*max(abs(im(:))))+.5); set(gca,'XTick',[],'YTick',[]); axis square; 
end
idx = COind(32);
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
wt = Mrgbtocc*Output_List{idx,5}; wt = sign(wt(2))*wt;
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
mean_expvariance_max = 100*mean(expvariance_max,1);
std_expvariance_max = 100*std(expvariance_max,1);
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

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% % calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];


figure(plot_counter); set(gcf,'Name','Cone wts') 
subplot(211);plot(conewts_svd(1,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(1,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
subplot(212); plot(conewts_svd(3,LumIds_conewts),conewts_svd(2,LumIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(3,ColorOpponentIds_conewts),conewts_svd(2,ColorOpponentIds_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(3,Sconedominated_conewts),conewts_svd(2,Sconedominated_conewts),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
plot(conewts_svd(3,Other_conewts),conewts_svd(2,Other_conewts),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('S'), ylabel('M');
plot_counter = plot_counter + 1;

% Further classification of S-cone sensitive cells
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
SvsLM = Sconesensitive(:,sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==-1);
SMvsL = Sconesensitive(:,sign(Sconesensitive(1,:))==-1 & sign(Sconesensitive(3,:))==1);
SLvsM = Sconesensitive(:,sign(Sconesensitive(1,:))==-1 & sign(Sconesensitive(3,:))==-1);

figure(plot_counter); set(gcf,'Name','Cone wts for S cells') 
plot(SvsLM(1,:),SvsLM(2,:),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.1],'MarkerEdgeColor',[1 1 1]); hold on;
plot(SMvsL(1,:),SMvsL(2,:),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0.7 0.3 0.1],'MarkerEdgeColor',[1 1 1]);
plot(SLvsM(1,:),SLvsM(2,:),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0.3 0.7 0.1],'MarkerEdgeColor',[1 1 1]); 
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  xlabel('L'), ylabel('M');
plot_counter = plot_counter + 1;

%% Figure 3: Model fits for each cell type
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
ind = find(~Singleopponent & simplecells);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
Lumind = ind(LumIds_conewts); 
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% idx = [Lumind(8); COind(30); Sconeind(20)];
idx = [Lumind(20); Lumind(2); COind(32); COind(27); Sconeind(15); Sconeind(8)];
figure(plot_counter);
L = numel(idx);
swt_tmp = [-1 -1 1 1 -1 1];
gf_tmp = [-1 1 1 1 -1 -1];
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
    % Cone weights
    wt = Mrgbtocc * Output_List{jj,5};
    wt = sign(wt(2))*wt;
    wt = wt/sum(abs(wt));
    subplot(5,L,ii+L); bar(wt','FaceColor','k'); set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square;
    % Spatial weighting function
    im = swt_tmp(ii)*Output_List{jj,4};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+2*L); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    % Gabor fit
    im = gf_tmp(ii)*modelfits{jj,2};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+3*L);image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
    % DoG fit
    im = -1*dogf_tmp(ii)*modelfits{jj,3};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+4*L); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));
end
plot_counter = plot_counter + 1;

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

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
SNR = Zmax(Z_cellsofinterest)./Z(Z_cellsofinterest,1);
SNR = SNR(~Singleopponent & simplecells);

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

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

idx = [20; 2; 32; 27; 18; 9];
bins2 = -10:1:10;
% Based on CV-R
figure(plot_counter); set(gcf,'Name','Model Comparison: DoG vs Gabor');
subplot(311);  hold on;
for ii = 1:length(LumIds_conewts)
    if ii ~=20 | ii ~=2
        h(ii) = plot(meanR(LumIds_conewts(ii),2),meanR(LumIds_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 20
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
    if ii ~=30 | ii ~=25
        h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),2),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 30
        h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),2),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 25
        h(ii) = plot(meanR(ColorOpponentIds_conewts(ii),2),meanR(ColorOpponentIds_conewts(ii),3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
end
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

subplot(313), hold on;
for ii = 1:length(Sconedominated_conewts)
    if ii ~=15 | ii ~=8
        h(ii) = plot(meanR(Sconedominated_conewts(ii),2),meanR(Sconedominated_conewts(ii),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 15
        h(ii) = plot(meanR(Sconedominated_conewts(ii),2),meanR(Sconedominated_conewts(ii),3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', num2str(ii),''')']);
    end
    if ii == 8
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

SNR = Zmax(Z_cellsofinterest);
% SNR = cell2mat(Output_List(:,20)); % See PopAnalysis_createOutputList.m for what is being stored in the 20th column
SNR = SNR(~Singleopponent & simplecells);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

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
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
set(gca,'Tickdir','out','Xlim',[2.5 3.7],'Ylim',[0 1],'XTick',2.5:0.4:3.7,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(122); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[2.5 3.7],'Ylim',[0 1],'XTick',2.5:0.4:3.7,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;
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

% Need to quantify the the differences in DoG as a function of SNR:
celltype_dv = zeros(numel(group),3);
for ii = 1:size(celltype_dv,1)
    celltype_dv(ii,group(ii)) = 1;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
colors = {'green','red','blue'};
 
% Trying linear regression
figure(plot_counter);
subplot(121); plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,2),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,2),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); lims = get(gca,'Xlim');
X = [repmat(data1,1,3).*celltype_dv celltype_dv]; % Individual model
[b, bint, r_Gabor,~, stats] = regress(data2,X);
[b_full, bint_full, r_full_Gabor] = regress(data2,[data1 ones(size(data1))]);
for i = 1:3
    tmp = [b(i) b(i+3)]*[lims; 1 1];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[lims; 1 1],'k'); hold on;
axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;

subplot(122); plot(log10(SNR(LumIds_conewts)),meanR(LumIds_conewts,3),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(log10(SNR(ColorOpponentIds_conewts)),meanR(ColorOpponentIds_conewts,3),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(log10(SNR(Sconedominated_conewts)),meanR(Sconedominated_conewts,3),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); 
X = [repmat(data1,1,3).*celltype_dv celltype_dv];
[b, bint, r_DoG,~, stats] = regress(data3,X);
[b_full, bint_full, r_full_DoG] = regress(data3,[data1 ones(size(data1))]);
for i = 1:3
    tmp = [b(i) b(i+3)]*[lims; 1 1];
    plot(lims,tmp,'color',colors{i}); hold on;
end
plot(lims,[b_full(1) b_full(2)]*[lims; 1 1],'k'); hold on;
axis square; xlabel('log10(SNR)'); ylabel('R'); title('DoG'); hold off;
plot_counter = plot_counter + 1;

RSSE_Gabor = sum(r_Gabor.^2);
RSSE_full_Gabor = sum(r_full_Gabor.^2);
Fstat_Gabor = ((RSSE_full_Gabor - RSSE_Gabor)/4)/(RSSE_Gabor/(numel(data1)-6));
p_Gabor = 1-fcdf(Fstat_Gabor,4,(numel(data1)-6));

RSSE_DoG = sum(r_DoG.^2);
RSSE_full_DoG = sum(r_full_DoG.^2);
Fstat_DoG = ((RSSE_full_DoG - RSSE_DoG)/4)/(RSSE_DoG/(numel(data1)-6));
p_DoG = 1-fcdf(Fstat_DoG,4,(numel(data1)-6));

% Trying 2D KS test 
alpha = 0.05/3; % Applying Bonferroni correction
% Looking at Gabor differences
[H1, pValue1, KSstatistic1] = kstest_2s_2d([log10(SNR(LumIds_conewts)) meanR(LumIds_conewts,2)], [log10(SNR(ColorOpponentIds_conewts)) meanR(ColorOpponentIds_conewts,2)], alpha);
[H2, pValue2, KSstatistic2] = kstest_2s_2d([log10(SNR(LumIds_conewts)) meanR(LumIds_conewts,2)], [log10(SNR(Sconedominated_conewts)) meanR(Sconedominated_conewts,2)], alpha);
[H3, pValue3, KSstatistic3] = kstest_2s_2d([log10(SNR(Sconedominated_conewts)) meanR(Sconedominated_conewts,2)], [log10(SNR(ColorOpponentIds_conewts)) meanR(ColorOpponentIds_conewts,2)], alpha);

% Looking at DoG differences
[H4, pValue4, KSstatistic4] = kstest_2s_2d([log10(SNR(LumIds_conewts)) meanR(LumIds_conewts,3)], [log10(SNR(ColorOpponentIds_conewts)) meanR(ColorOpponentIds_conewts,3)], alpha);
[H5, pValue5, KSstatistic5] = kstest_2s_2d([log10(SNR(LumIds_conewts)) meanR(LumIds_conewts,3)], [log10(SNR(Sconedominated_conewts)) meanR(Sconedominated_conewts,3)], alpha);
[H6, pValue6, KSstatistic6] = kstest_2s_2d([log10(SNR(Sconedominated_conewts)) meanR(Sconedominated_conewts,3)], [log10(SNR(ColorOpponentIds_conewts)) meanR(ColorOpponentIds_conewts,3)], alpha);

% Fitting a saturating exponential to the R_Gabor 
f = @(b,x) b(3)./(1+exp(b(1)*(x-b(2)))); % Objective Function
% B1 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),[-10 -5 eps],[-eps 5 1]);
% B2 = lsqcurvefit(f,[-1; 2; 0.5],log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),[-10 -5 eps],[-eps 5 1]);
B1 = fitsigmoidMLE(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),[meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2)*10000 10000-10000*meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2)]);
B2 = fitsigmoidMLE(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),[meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3)*10000 10000-10000*meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3)]);

figure(plot_counter); 
subplot(121); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],2),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(2.4,3.6,51),f(B1,linspace(2.4,3.6,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
subplot(122); plot(log10(SNR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts])),meanR([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts],3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
hold on; plot(linspace(2.4,3.6,51),f(B2,linspace(2.4,3.6,51)),'k','Linewidth',2);
set(gca,'Tickdir','out','Xlim',[2.4 3.6],'Ylim',[0 1],'XTick',2.4:0.4:3.6,'YTick',0:0.25:1.0); axis square; xlabel('log10(SNR)'); ylabel('R'); title('Gabor'); hold off;
plot_counter = plot_counter + 1;

Max_GaborR = B1(3);
Max_DoGR = B2(3);
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

numcells = size(Output_List,1);
% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

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
Lumid = logical(Lumid); COid = logical(COid); Sid = logical(Sid);

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
group = [ones(size(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)')); 2*ones(size(aspectratio_relevantcells(COid & meanR_GaborDOG>0)')); 3*ones(size(aspectratio_relevantcells(Sid & meanR_GaborDOG>0)'))];
p = kruskalwallis(data,group,'off');

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

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

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

idx = [Lumind(20); Lumind(2); COind(32); COind(27); Sconeind(15); Sconeind(8)];
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

% Analyzing spatial phases
% storing gabor phases of the cells
numcells = size(Output_List,1);
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

meanR_GaborCrescent = meanR(:,2)-meanR(:,1); % Gabor - Crescent
Lumid = zeros(size(aspectratio_relevantcells))'; Lumid(LumIds_conewts) = 1;
COid = zeros(size(aspectratio_relevantcells))'; COid(ColorOpponentIds_conewts) = 1;
Sid = zeros(size(aspectratio_relevantcells))'; Sid(Sconedominated_conewts) = 1;

% Mann-Whitney U test for aspect analysing ratios
p7 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborCrescent>0),aspectratio_relevantcells(Lumid & meanR_GaborCrescent<0));
p8 = ranksum(aspectratio_relevantcells(COid & meanR_GaborCrescent>0),aspectratio_relevantcells(COid & meanR_GaborCrescent<0));
p9 = ranksum(aspectratio_relevantcells(Sid & meanR_GaborCrescent>0),aspectratio_relevantcells(Sid & meanR_GaborCrescent<0));

p10 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborCrescent>0),aspectratio_relevantcells(COid & meanR_GaborCrescent>0));
p11 = ranksum(aspectratio_relevantcells(Sid & meanR_GaborCrescent>0),aspectratio_relevantcells(COid & meanR_GaborCrescent>0));
p12 = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborCrescent>0),aspectratio_relevantcells(Sid & meanR_GaborCrescent>0));

data = [aspectratio_relevantcells(Lumid & meanR_GaborCrescent>0)';aspectratio_relevantcells(COid & meanR_GaborCrescent>0)';aspectratio_relevantcells(Sid & meanR_GaborCrescent>0)'];
group = [ones(size(aspectratio_relevantcells(Lumid & meanR_GaborCrescent>0)')); 2*ones(size(aspectratio_relevantcells(COid & meanR_GaborCrescent>0)')); 3*ones(size(aspectratio_relevantcells(Sid & meanR_GaborCrescent>0)'))];
p13 = kruskalwallis(data,group);


% Mann-Whitney U test for phases 
p14 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborCrescent>0)),90-abs(90-gaborphases_relevantcells(COid & meanR_GaborCrescent>0)));
p15 = ranksum(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborCrescent>0)),90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborCrescent>0)));
p16 = ranksum(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborCrescent>0)),90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborCrescent>0)));

%% Additional figure 1: STAs of L+M, L-M & S-cone sensitive cells
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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];


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

%% Additional figure 2: Single opponent cells, Cone weights and STAs

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

RGB_svd = cell2mat(Output_List(Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svdSO = Mrgbtocc * RGB_svd;
conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

% calculating the cone weights
% conewts_svdSO = cell2mat(Output_List(Singleopponent & simplecells,23)');
% conewts_svdSO = conewts_svdSO./repmat(sum(abs(conewts_svdSO),1),[3 1]);
% conewts_svdSO = conewts_svdSO .* repmat(sign(conewts_svdSO(2,:)),[3 1]);

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

%% Additional figure 3: Model comparison of DOG vs Gabor (CV per accuracy, mean SSE & BIC

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

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


%% Additional figure 4: Model comparison of Crescent vs Gabor (CV per accuracy, mean SSE & BIC)

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

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
bkgndlms = cell2mat(Output_List(1,16)');
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');
conewts_svd = Mrgbtocc * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% calculating the cone weights
% conewts_svd = cell2mat(Output_List(~Singleopponent & simplecells,23)');
% conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
% conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

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


