% Figures used for Jneurophys rebuttal
% Author - Abhishek De, 

close all; clearvars;
plot_counter = 1;
% Section 1: Scatterplot showing the DoG and the crescent shapped model (reviewer's eyes only)

% Section 2: Scatterplot of RF size as a function with eccentricity (reviewer's eyes only)

% Section 3: Simulating Conway experiment, testing the Goldman-Lankow explanation of double-opponency. (Greg's code)

% Section 4: Gaussian RGB white noise in other color spaces (Greg's code)

% Section 5: Are the results consistent across the monkey genders?

% Section 6: Checking whether there are S-cone weights and BY DO cells within the simple cells

% Section 7: Trying a new criteria, simple cells -> cone non-opponent, DO cells -> cone opponent

% Section 8: Some more analysis of the spatial weighting function to address the reviewer's concern on BY DO cells 

% Section 9: Comparison of DOG gamma model with DOG, Gabor and non-concetric DOG model

% Section 10: Elliptical DOG (gamma) model fits for example cells (could be a reviewer figure)

% Section 11: Calculating a new RF area suggested by Greg (pi*sigma.^2/gamma: area of 1 STD ellipse)

% Section 12: A different representation for Figure 6, as per Reviewer 2's suggestion

%% Section 1: Scatterplot showing the DoG and the crescent shaped model (reviewer's eyes only)

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

%% Section 2: Scatterplot of RF size as a function with eccentricity (reviewer's eyes only)

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

idxs = [Lumind; COind; Sconeind];

ecc = [];
RFarea = [];
for aa = 1:numel(simplecells)
    
    ii = aa;
    % The signficant RFs are stored in column 18
    RFarea = [RFarea; sum(Output_List{ii,18}(:))*0.04];
    
    % RF eccentricity is stored in column 8
    tmpecc = cell2mat(Output_List(ii,8));
    ecc = [ecc; sqrt(sum(tmpecc.^2,2))/10];

end

% plotting the results 
figure(plot_counter);
subplot(131); plot(ecc(Lumind), sqrt(RFarea(Lumind)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Ecc'); ylabel('RF field size'); title('Simple'), set(gca,'Xlim',[0 8],'Ylim',[0 2]);
subplot(132); plot(ecc(COind), sqrt(RFarea(COind)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Ecc'); ylabel('RF field size'); title('DO LM'); set(gca,'Xlim',[0 8], 'Ylim',[0 2]);
subplot(133); plot(ecc(Sconeind), sqrt(RFarea(Sconeind)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('Ecc'); ylabel('RF field size'); title('DO S'); set(gca,'Xlim',[0 8],'Ylim', [0 2]);
plot_counter = plot_counter + 1;

% Some stats 
[r1,p1] = corr(ecc(Lumind), sqrt(RFarea(Lumind)),'type','Spearman','rows','complete');
[r2,p2] = corr(ecc(COind), sqrt(RFarea(COind)),'type','Spearman','rows','complete');
[r3,p3] = corr(ecc(Sconeind), sqrt(RFarea(Sconeind)),'type','Spearman','rows','complete');


%% Section 3: %% Simulating Conway experiment, testing the Goldman-Lankow explanation of double-opponency.
% Simulated L-cone map
sigma = .5; % Standard deviation of single-opponent (L-M) RF
bounds = 3;
boxwidth = .5;
niter = 2000;
data = [];
for i = 1:niter
    inc_pos = unifrnd(-bounds,bounds,2,1);
    dec_pos = unifrnd(-bounds,bounds,2,1);
    while abs(inc_pos(1)-dec_pos(1))<boxwidth & abs(inc_pos(2)-dec_pos(2))<boxwidth
        dec_pos = unifrnd(-bounds,bounds,2,1); % No overlapping increments with decrements
    end
    r_inc_x = normcdf(inc_pos(1)+[-.5 .5]*boxwidth, 0, sigma);
    r_inc_y = normcdf(inc_pos(2)+[-.5 .5]*boxwidth, 0, sigma);
    r_inc = diff(r_inc_x)*diff(r_inc_y);
    r_dec_x = normcdf(dec_pos(1)+[-.5 .5]*boxwidth, 0, sigma);
    r_dec_y = normcdf(dec_pos(2)+[-.5 .5]*boxwidth, 0, sigma);
    r_dec = diff(r_dec_x)*diff(r_dec_y);
    r = r_inc-r_dec;
    data = [data; r inc_pos' dec_pos'];
end
% Calculating STA
[x,y] = meshgrid(linspace(-bounds, bounds,20));
inc_map = zeros(size(x));
dec_map = zeros(size(x));
for i = 1:numel(x)
    L_inc = abs(data(:,2)-x(i)) < boxwidth & abs(data(:,3)-y(i)) < boxwidth;
    L_dec = abs(data(:,4)-x(i)) < boxwidth & abs(data(:,5)-y(i)) < boxwidth;
    inc_map(i) = (data(:,1)'*L_inc);
    dec_map(i) = (data(:,1)'*L_dec);
end
% For image rendering, putting min at 0 and max at 255
min_val = min([inc_map(:); dec_map(:)]);
max_val = max([inc_map(:); dec_map(:)]);
figure; 
colormap(jet(255))
subplot(2,2,1);
image(255*(inc_map-min_val)./(max_val-min_val)); axis square; set(gca,'XTick',[],'YTick',[]);
subplot(2,2,2);
image(255*(dec_map-min_val)./(max_val-min_val)); axis square; set(gca,'XTick',[],'YTick',[]);
subplot(2,2,3); hold on;
plot(x(1,:),inc_map(10,:),'r-');
plot(x(1,:),dec_map(10,:),'b-');
hold off;
set(gcf,'renderer','painters');

%% Section 4: % Gaussian RGB white noise in other color spaces

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
figure;
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

%% Section 5: Are the results consistent across the monkey genders?
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

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
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
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

% Need to classify files as male vs. female
Gender = []; % Male = 0, Female = 1;

for ii = 1:numel(simplecells)
    if strcmp(char(Output_List{ii,1}(1)),'S') | strcmp(char(Output_List{ii,1}(1)),'K')
        Gender = [Gender; 1]; % Female
    else
        Gender = [Gender; 0]; % Male
    end
end
Gender = logical(Gender(ind));

figure(plot_counter); set(gcf,'Name','Model Comparison: Gabor vs. DoG');
subplot(231);  hold on;
plot(meanR(LumIds_conewts(Gender(LumIds_conewts)),2),meanR(LumIds_conewts(Gender(LumIds_conewts)),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(meanR(LumIds_conewts(~Gender(LumIds_conewts)),2),meanR(LumIds_conewts(~Gender(LumIds_conewts)),3),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');

subplot(232), hold on;
plot(meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),3),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

subplot(233), hold on;
plot(meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),3),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('DOG'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('S');

subplot(234);  hold on;
plot(meanR(LumIds_conewts(Gender(LumIds_conewts)),2),meanR(LumIds_conewts(Gender(LumIds_conewts)),1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(meanR(LumIds_conewts(~Gender(LumIds_conewts)),2),meanR(LumIds_conewts(~Gender(LumIds_conewts)),1),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('Lum');

subplot(235), hold on;
plot(meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),1),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('L-M');

subplot(236), hold on;
plot(meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),1),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),1),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
line([0 1],[0 1],'Color','k');  xlabel('Gabor'); ylabel('Crescent'); axis square; set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1],'XTick',0:0.5:1.0,'YTick',0:0.5:1.0); title('S');

plot_counter = plot_counter + 1;

DO_conewts = [ColorOpponentIds_conewts Sconedominated_conewts];
% Some additional stats test 
% Gabor vs. DoG
p1 = signrank(meanR(LumIds_conewts(Gender(LumIds_conewts)),2),meanR(LumIds_conewts(Gender(LumIds_conewts)),3)); % Male
p2 = signrank(meanR(LumIds_conewts(~Gender(LumIds_conewts)),2),meanR(LumIds_conewts(~Gender(LumIds_conewts)),3)); % Female
p3 = signrank(meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),3)); % Male
p4 = signrank(meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),3)); % Female
p5 = signrank(meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),3)); % Male
p6 = signrank(meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),3)); % Female

% Gabor vs. crescent
p7 = signrank(meanR(LumIds_conewts(Gender(LumIds_conewts)),2),meanR(LumIds_conewts(Gender(LumIds_conewts)),1)); % Male
p8 = signrank(meanR(LumIds_conewts(~Gender(LumIds_conewts)),2),meanR(LumIds_conewts(~Gender(LumIds_conewts)),1)); % Female
p9 = signrank(meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)),1)); % Male
p10 = signrank(meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),2),meanR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts)),1)); % Female
p11 = signrank(meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(Gender(Sconedominated_conewts)),1)); % Male
p12 = signrank(meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),2),meanR(Sconedominated_conewts(~Gender(Sconedominated_conewts)),1)); % Female


% Combining the DO cells: LM opponent and S-cone dominated 
p13 = signrank(meanR(DO_conewts(Gender(DO_conewts)),2),meanR(DO_conewts(Gender(DO_conewts)),3)); % Female
p14 = signrank(meanR(DO_conewts(~Gender(DO_conewts)),2),meanR(DO_conewts(~Gender(DO_conewts)),3)); % Male
p15 = signrank(meanR(DO_conewts(Gender(DO_conewts)),2),meanR(DO_conewts(Gender(DO_conewts)),1)); % Female
p16 = signrank(meanR(DO_conewts(~Gender(DO_conewts)),2),meanR(DO_conewts(~Gender(DO_conewts)),1)); % Male

% Checking the differences in SNR  between the male and female macaques
p17 = ranksum(SNR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts))),SNR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts)))); % L-M
p18 = ranksum(SNR(LumIds_conewts(~Gender(LumIds_conewts))),SNR(LumIds_conewts(Gender(LumIds_conewts)))); % L+M
p19 = ranksum(SNR(Sconedominated_conewts(~Gender(Sconedominated_conewts))),SNR(Sconedominated_conewts(Gender(Sconedominated_conewts)))); % L+M

% A new figure for showing the SNR distribution 
figure(plot_counter);
bins = logspace(0,log10(16),8); 
subplot(131); histogram(SNR(LumIds_conewts(~Gender(LumIds_conewts))), bins, 'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(SNR(LumIds_conewts(Gender(LumIds_conewts))), bins, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); axis square; 
xlabel('SNR'); ylabel('# cells'), set(gca,'Tickdir','out','XScale','log','Xlim',[1 16],'XTick',[1 2 4 8 16]); title('LUM'); legend('Male','Female'); hold off;
subplot(132); histogram(SNR(ColorOpponentIds_conewts(~Gender(ColorOpponentIds_conewts))), bins, 'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(SNR(ColorOpponentIds_conewts(Gender(ColorOpponentIds_conewts))), bins, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); axis square; 
xlabel('SNR'); ylabel('# cells'), set(gca,'Tickdir','out','XScale','log','Xlim',[1 16],'XTick',[1 2 4 8 16]); title('L-M'); hold off;
subplot(133); histogram(SNR(Sconedominated_conewts(~Gender(Sconedominated_conewts))), bins, 'DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
histogram(SNR(Sconedominated_conewts(Gender(Sconedominated_conewts))), bins, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); axis square; 
xlabel('SNR'); ylabel('# cells'), set(gca,'Tickdir','out','XScale','log','Xlim',[1 16],'XTick',[1 2 4 8 16]); title('S'); hold off;
plot_counter = plot_counter + 1;

% A statistical comparison suggested by Greg
p20 = ranksum(meanR(DO_conewts(Gender(DO_conewts)),2)-meanR(DO_conewts(Gender(DO_conewts)),3),meanR(DO_conewts(~Gender(DO_conewts)),2)-meanR(DO_conewts(~Gender(DO_conewts)),3)); % Female vs. Male
p21 = ranksum(meanR(DO_conewts(Gender(DO_conewts)),2)-meanR(DO_conewts(Gender(DO_conewts)),1),meanR(DO_conewts(~Gender(DO_conewts)),2)-meanR(DO_conewts(~Gender(DO_conewts)),1)); % Female vs. Male

%% Section 6: Checking whether there are S-cone weights and BY DO cells within the simple cells
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

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
COind = ind(ColorOpponentIds_conewts); 
Sconeind = ind(Sconedominated_conewts);

% Weighted STA is stored in the second column
BRgun_corr = []; % For storing the Blue-yellow STA 
for ii = 1:numel(Lumind)
    STA = Output_List{Lumind(ii),2};
    [r,p] = corr(STA(:,1),STA(:,3),'type','Spearman');
    BRgun_corr = [BRgun_corr; r p conewts_svd(3,LumIds_conewts(ii))];
end

%% Section 7: Trying a new criteria, simple cells -> cone non-opponent, DO cells -> cone opponent 
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

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts);
DOind = ind(DOIds_conewts); 

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

% Doing some stats
% DoG vs. Gabor 
p1 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,3));
p2 = signrank(meanR(DOIds_conewts,2),meanR(DOIds_conewts,3));

% Crescent vs. Gabor 
p3 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,1));
p4 = signrank(meanR(DOIds_conewts,2),meanR(DOIds_conewts,1));

% Gabor phases and aspect ratios 
Lumid = zeros(size(gaborphases_relevantcells))'; Lumid(LumIds_conewts) = 1;
DOid = zeros(size(gaborphases_relevantcells))'; DOid(DOIds_conewts) = 1;
Lumid = logical(Lumid); DOid = logical(DOid); 

LUM_phase = mean(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)));
DO_phase = mean(90-abs(90-gaborphases_relevantcells(DOid & meanR_GaborDOG>0)));

LUM_aspectratio = median(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0));
DO_aspectratio = median(aspectratio_relevantcells(DOid & meanR_GaborDOG>0));

p = ranksum(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0), aspectratio_relevantcells(DOid & meanR_GaborDOG>0));

%% Section 8:  Some more analysis of the spatial weighting function to address the reviewer's concern on BY DO cells 

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


% Z-transforming the R values
data1 = [SNR(LumIds_conewts); SNR(ColorOpponentIds_conewts); SNR(Sconedominated_conewts)];
data1log = log10(data1);
data2 = atanh([meanR(LumIds_conewts,2); meanR(ColorOpponentIds_conewts,2); meanR(Sconedominated_conewts,2)]); % Gabor 
data3 = atanh([meanR(LumIds_conewts,3); meanR(ColorOpponentIds_conewts,3); meanR(Sconedominated_conewts,3)]); % DoG
group = [ones(size(SNR(LumIds_conewts))); 2*ones(size(SNR(ColorOpponentIds_conewts))) ; 3*ones(size(SNR(Sconedominated_conewts)))];

% Comparing the SNRs of LM opponent and S-cone sensitive cells 
p1 = ranksum(log10(SNR(Sconedominated_conewts)), log10(SNR(ColorOpponentIds_conewts)));


% A different way of calculating the SNR 
N = cell2mat(Output_List(:,21)); % Noise 
S = cell2mat(Output_List(:,17)); % Signal
SNR_new = S(:,1)./N(:,1);
p2 = ranksum(log10(SNR_new(Sconedominated_conewts)), log10(SNR_new(ColorOpponentIds_conewts)));

%% Section 9: Comparison of DOG gamma model with DOG, Gabor and non-concetric DOG model

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

load Peraccuracy_gamma.mat
load SSE_gamma.mat
load Deviation_gamma.mat
meanperaccuracy = zeros(size(Peraccuracy_gamma));
meanSSE = zeros(size(SSE_gamma));
meanR = zeros(size(SSE_gamma));
for ii = 1:size(Peraccuracy_gamma,1)
    for jj = 1:size(Peraccuracy_gamma,2)
        meanperaccuracy(ii,jj) = mean(Peraccuracy_gamma{ii,jj});
        meanSSE(ii,jj) = mean(SSE_gamma{ii,jj});
        meanR(ii,jj) = mean(cos(Deviation_gamma{ii,jj}*pi/180));
    end
end
meanperaccuracy = meanperaccuracy(~Singleopponent & simplecells,:);
meanSSE = meanSSE(~Singleopponent & simplecells,:);
meanR = meanR(~Singleopponent & simplecells,:);

DOIds_conewts = [ColorOpponentIds_conewts Sconedominated_conewts];
% Doing some stats
p1 = signrank(meanR(LumIds_conewts,2),meanR(LumIds_conewts,6)); % Gabor vs. DOGgamma
p2 = signrank(meanR(LumIds_conewts,3),meanR(LumIds_conewts,6)); % DOG vs. DOGgamma
p3 = signrank(meanR(LumIds_conewts,1),meanR(LumIds_conewts,6)); % Crescent vs. DOGgamma

p4 = signrank(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,6)); % Gabor vs. DOGgamma
p5 = signrank(meanR(ColorOpponentIds_conewts,3),meanR(ColorOpponentIds_conewts,6)); % DOG vs. DOGgamma
p6 = signrank(meanR(ColorOpponentIds_conewts,1),meanR(ColorOpponentIds_conewts,6)); % Crescent vs. DOGgamma

p7 = signrank(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,6)); % Gabor vs. DOGgamma
p8 = signrank(meanR(Sconedominated_conewts,3),meanR(Sconedominated_conewts,6)); % DOG vs. DOGgamma
p9 = signrank(meanR(Sconedominated_conewts,1),meanR(Sconedominated_conewts,6)); % Crescent vs. DOGgamma

p10 = signrank(meanR(DOIds_conewts,2),meanR(DOIds_conewts,6)); % Gabor vs. DOGgamma
p11 = signrank(meanR(DOIds_conewts,3),meanR(DOIds_conewts,6)); % DOG vs. DOGgamma
p12 = signrank(meanR(DOIds_conewts,1),meanR(DOIds_conewts,6)); % Crescent vs. DOGgamma
p13 = signrank(meanR(DOIds_conewts,1),meanR(DOIds_conewts,2)); % Crescent vs. Gabor

% Gabor vs. Crescent
p14 = signrank(meanR(LumIds_conewts,1),meanR(LumIds_conewts,2)); % Crescent vs. Gabor

% Plotting the results 
figure(plot_counter);
subplot(331); hold on; plot(meanR(LumIds_conewts,2),meanR(LumIds_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Gabor'); ylabel('Gamma'); hold off; 
subplot(332); hold on; plot(meanR(LumIds_conewts,3),meanR(LumIds_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('DOG'); ylabel('Gamma'); hold off;
subplot(333); hold on; plot(meanR(LumIds_conewts,1),meanR(LumIds_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Crescent'); ylabel('Gamma'); hold off;
subplot(334); hold on; plot(meanR(ColorOpponentIds_conewts,2),meanR(ColorOpponentIds_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Gabor'); ylabel('Gamma'); hold off; 
subplot(335); hold on; plot(meanR(ColorOpponentIds_conewts,3),meanR(ColorOpponentIds_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('DOG'); ylabel('Gamma'); hold off;
subplot(336); hold on; plot(meanR(ColorOpponentIds_conewts,1),meanR(ColorOpponentIds_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Crescent'); ylabel('Gamma'); hold off;
subplot(337); hold on; plot(meanR(Sconedominated_conewts,2),meanR(Sconedominated_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Gabor'); ylabel('Gamma'); hold off; 
subplot(338); hold on; plot(meanR(Sconedominated_conewts,3),meanR(Sconedominated_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('DOG'); ylabel('Gamma'); hold off;
subplot(339); hold on; plot(meanR(Sconedominated_conewts,1),meanR(Sconedominated_conewts,6),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
axis square; plot([0 1],[0 1],'k'); set(gca,'Tickdir','out','Xlim',[0 1],'Ylim',[0 1]); xlabel('Crescent'); ylabel('Gamma'); hold off;
plot_counter = plot_counter + 1;

%% Section 10:  Elliptical DOG (gamma) model fits for example cells (could be a reviewer figure)

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
    im = -1*dogf_tmp(ii)*modelfits{jj,3};5
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+3*L); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));


     % DoG gamma fit
    im = -1*dogf_tmp(ii)*modelfits{jj,6};
    im = sigmoid(im,lambda,0)-0.5;
    subplot(5,L,ii+4*L); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255));

end
set(gcf,'renderer','painters');
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

%% Section 11: Calculating a new RF area suggested by Greg (pi*sigma.^2/gamma: area of 1 STD ellipse)

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
subplot(121); plot(ecc([Lumind; COind; Sconeind]),RFarea_gaborfit([Lumind; COind; Sconeind]),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out', 'Xlim', [0 10],'YScale','log','Ylim',[0.005 10]); title('All cells'); ylabel('RF area'); xlabel('Eccentricity in degrees')
subplot(122); plot(ecc(Lumind),RFarea_gaborfit(Lumind),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(ecc(COind),RFarea_gaborfit(COind),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(ecc(Sconeind),RFarea_gaborfit(Sconeind),'o','MarkerFaceColor',[0 0.5 1.0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out', 'Xlim', [0 10],'YScale','log','Ylim',[0.005 10]);  ylabel('RF area'); 
legend('LUM','L-M','S'); xlabel('Eccentricity in degrees');
plot_counter = plot_counter + 1;

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

%% Section 12: A different representation for Figure 6, as per Reviewer 2's suggestion

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
Lumid = logical(Lumid); COid = logical(COid); Sid = logical(Sid);

% All cells and Best fitting Gabor cells
figure(plot_counter); set(gcf,'Name','Gabor Phases & aspect ratios');
subplot(321),histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG<0)),bins1,'DisplayStyle','stairs','EdgeColor',[0 0.5 1]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(Lumid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 30],'YTick',0:10:30); axis square; title('Lum'); 
subplot(322),histogram(aspectratio_relevantcells(Lumid & meanR_GaborDOG<0),bins2,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanR_GaborDOG<0)),30,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Lumid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Ylim',[0 30],'YTick',0:15:30,'Xlim',[0.3 10],'XTick',[0.3 1 3 10], 'Xscale','log'); axis square; title('Lum'); 
subplot(323),histogram(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG<0)),bins1,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(COid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 30],'YTick',0:10:30); axis square; title('L-M'); 
subplot(324),histogram(aspectratio_relevantcells(COid & meanR_GaborDOG<0),bins2,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0]);  hold on; plot(median(aspectratio_relevantcells(COid & meanR_GaborDOG<0)),30,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(COid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(COid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0.3 10],'XTick',[0.3 1 3 10],'Ylim',[0 30],'YTick',0:15:30,'Xscale','log'); axis square;title('L-M'); 
subplot(325),histogram(90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG<0)),bins1,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0]); hold on;
histogram(90-abs(90-gaborphases_relevantcells(Sid & meanR_GaborDOG>0)),bins1,'FaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0 90],'XTick',0:45:90,'Ylim',[0 10],'YTick',0:5:10); axis square; title('S');
subplot(326),histogram(aspectratio_relevantcells(Sid & meanR_GaborDOG<0),bins2,'DisplayStyle','stairs','EdgeColor',[0 0.5 1.0]); hold on; plot(median(aspectratio_relevantcells(Sid & meanR_GaborDOG<0)),15,'kv','MarkerFaceColor',[1 1 1]);
histogram(aspectratio_relevantcells(Sid & meanR_GaborDOG>0),bins2,'FaceColor',[0 0 0]); hold on; plot(median(aspectratio_relevantcells(Sid & meanR_GaborDOG>0)),0,'kv','MarkerFaceColor',[0 0 0]); set(gca,'Tickdir','out','Xlim',[0.3 10],'XTick',[0.3 1 3 10],'Ylim',[0 15],'YTick',0:5:15,'Xscale','log'); axis square; title('S');
set(gcf,'renderer','painters');
plot_counter = plot_counter + 1;

