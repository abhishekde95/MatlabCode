% Writing this script to come up with better ways of identifying single
% opponent cells
% Author - Abhishek De, 4/20

close all; clearvars;
plot_counter = 1;

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
ind = find(~Singleopponent & simplecells);

load Deviation.mat
meanR = zeros(size(Deviation));
Ratio_of_power = [];
pval = [];
for ii = 1:size(Output_List,1)
    im = Output_List{ii,4};
    fim = fft2(im);
    Ratio_of_power = [Ratio_of_power; max(abs(fim(:)))/abs(fim(1))]; 
    
    for jj = 1:size(Deviation,2)
        meanR(ii,jj) = mean(cos(Deviation{ii,jj}*pi/180));
    end 
    
    % Calculating p value 
    stddevSTA = unique(Output_List{ii,27});
    STA = Output_List{ii,2};
    normSTA = norm(STA(:));
    nspikes = Output_List{ii,15};
    stdfft = stddevSTA*5/(normSTA*sqrt(nspikes)*sqrt(3));
    pval = [pval; 1-normcdf(max(abs(fim(:))),abs(fim(1)),stdfft)];
end

x = logspace(0,4,51);
figure(plot_counter); subplot(223); histogram(Ratio_of_power(simplecells & pval>=0.01),x); hold on;
histogram(Ratio_of_power(simplecells & pval<0.01),x); axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[1 2500]);
ylabel('# cells'); xlabel('Highest response divided by response in lowest freq'); legend('p>0.01','p<0.01');
subplot(221); plot(Ratio_of_power(simplecells),meanR(simplecells,3),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Ratio_of_power(simplecells),meanR(simplecells,4),'o','MarkerSize',4,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Ratio_of_power(simplecells),meanR(simplecells,2),'o','MarkerSize',4,'MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Ratio_of_power(simplecells),meanR(simplecells,1),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Ratio_of_power(simplecells),meanR(simplecells,5),'o','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[1 2500]); ylabel('R[DOG,SG,Gabor]'); xlabel('Highest response divided by response in lowest freq'); 
legend('DoG','SG','Gabor','Crescent','S2G'); hold off;
subplot(222); histogram(meanR(simplecells,3),0:0.05:1.0,'Displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2); hold on;
histogram(meanR(simplecells,4),0:0.05:1.0,'Displaystyle','stairs','EdgeColor',[0 1 0],'Linewidth',2); hold on;
histogram(meanR(simplecells,2),0:0.05:1.0,'Displaystyle','stairs','EdgeColor',[0 0 1],'Linewidth',2); hold on;
histogram(meanR(simplecells,1),0:0.05:1.0,'Displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
histogram(meanR(simplecells,5),0:0.05:1.0,'Displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
axis square; set(gca,'Tickdir','out','Xlim',[0 1]); xlabel('R[DOG,SG,Gabor]'); ylabel('# cells');  
legend('DoG','SG','Gabor','Crescent','S2G'); hold off;


% Next I want to compare the in which frequency bin is the highest energy stored for single opponent and non-single opponent cells
conewts_svd = cell2mat(Output_List((~Singleopponent) & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

% Segregating cells into different categories
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); 
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Creating arrays for storing the spatial frequencies where the energy is highest 
Freq_SO = [];
Freq_nonSO = [];

% Creating arrays for storing ratio of powers
Ratio_of_power_SO = [];
Ratio_of_power_nonSO = [];

% Creating the spatial frequency grid 
f = (-5/10:1/10:4/10)*(2.5);
[Sx, Sy] = meshgrid(f,f);
Sf = sqrt(Sx.^2+Sy.^2);

AllIdxs = ind([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]);
for ii = 1:size(Output_List,1) 
    
    % Spatial RF
    im = Output_List{ii,4};
    fim = fftshift(fft2(im));
    
    [a,b] = find(abs(fim)==max(abs(fim(:))));
    tmp = unique(sqrt((a-6).^2 + (b-6).^2))*0.25;
    Rp = (max(abs(fim(:)))/abs(fim(6,6)));
    
    if numel(tmp)>1
        keyboard;
    end
    
    if Singleopponent(ii)
        Freq_SO = [Freq_SO; tmp];
        Ratio_of_power_SO = [Ratio_of_power_SO; Rp];
    end
    
    if ismember(ii,AllIdxs)
        Freq_nonSO = [Freq_nonSO; tmp];
        Ratio_of_power_nonSO = [Ratio_of_power_nonSO; Rp];
    end  
    
end

figure(plot_counter); subplot(224); plot(Freq_nonSO(pval(AllIdxs)<0.01),Ratio_of_power_nonSO(pval(AllIdxs)<0.01),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(Freq_nonSO(pval(AllIdxs)>=0.01),Ratio_of_power_nonSO(pval(AllIdxs)>=0.01),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
plot(Freq_SO,Ratio_of_power_SO,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([eps 1.5],[1 1],'k'); axis square; set(gca,'Tickdir','out','Xlim',[0 1.5],'XTick',[0:0.5:1.5],'YScale','log','Ylim',[0.5 10000]);
xlabel('Spatial frequency with highest response'); ylabel('Ratio of power'); legend('nonSO','SO pval','SO'); hold off;
plot_counter = plot_counter + 1;


% plotting the Ratio of power for LUM, DO and S-cone dominated cells 
figure(plot_counter); subplot(311); histogram(Ratio_of_power(ind(LumIds_conewts)),x,'DisplayStyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); hold on;
plot([1.2 1.2],[0 10],'k','Linewidth',2);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[1 2500]); ylabel('# cells'); xlabel('Highest response divided by response in lowest freq'); title('LUM');
subplot(312); histogram(Ratio_of_power(ind(ColorOpponentIds_conewts)),x,'DisplayStyle','stairs','EdgeColor',[1 0 0],'Linewidth',2); hold on;
plot([1.2 1.2],[0 20],'k','Linewidth',2);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[1 2500]); ylabel('# cells'); xlabel('Highest response divided by response in lowest freq'); title('DO(LM)');
subplot(313); histogram(Ratio_of_power(ind(Sconedominated_conewts)),x,'DisplayStyle','stairs','EdgeColor',[0 0.5 1],'Linewidth',2); hold on;
plot([1.2 1.2],[0 10],'k','Linewidth',2);
axis square; set(gca,'Tickdir','out','XScale','log','Xlim',[1 2500]); ylabel('# cells'); xlabel('Highest response divided by response in lowest freq'); title('DO(S)');
plot_counter = plot_counter + 1;
%% Plotting the spatial RF and their ffts for all single opponent cells:
% Data driven approach 
% Ratio_of_power<2 
if ~exist('plot_counter')
    plot_counter = 1;
end
lambda = 12;
CRITERIA = Ratio_of_power<2 & simplecells;
% N = find(Singleopponent==1)';
N = find(CRITERIA)';
[~,I] = sort(Ratio_of_power(N));
N = N(I);
count = 1;
numplots = ceil(sqrt(numel(N)));
for ii = N 
    
    % STA 
    STA = Output_List{ii,2};
    normfactor = 0.5/(max(abs(STA(:)))+0.01);
    STA = normfactor*STA + 0.5;
    STA = reshape(STA,[10 10 3]);
    
    % Spatial RF
    im = Output_List{ii,4};
    fim = fftshift(fft2(im));
    
    % Modifying im for plotting purposes
    im = sigmoid(im,lambda,0)-0.5;
    
   
    if CRITERIA(ii)
        figure(plot_counter); subplot(numplots,numplots,count); image(STA); set(gca,'XTick',[],'YTick',[]); axis square; title(num2str(Ratio_of_power(ii),3));
        if ismember(ii,AllIdxs)
            set(gca,'XColor','r','YColor','r');
        end
%         figure(plot_counter+1); subplot(numplots,numplots,count); imagesc(255*(im./(2*max(abs(im(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]);
        figure(plot_counter+2); subplot(numplots,numplots,count); imagesc(abs(fim)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title(num2str(Ratio_of_power(ii),3));
        if ismember(ii,AllIdxs)
            set(gca,'XColor','r','YColor','r');
        end
        count = count + 1;
    end
end
plot_counter = plot_counter + 3;

%% Checking the 2-D ffts of the cells where Gaussian fits are the superior
% Model based approach of segregating cells 
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
ind = find(~Singleopponent & simplecells);

conewts_svd = cell2mat(Output_List((~Singleopponent) & simplecells,23)');
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); 
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Creating the spatial frequency grid 
f = (-5/10:1/10:4/10)*2.5;
[Sx, Sy] = meshgrid(f,f);
Sf = sqrt(Sx.^2+Sy.^2);
Sf(Sf==0) = eps;

lambda = 12;
count = 1;
[~,I] = max(meanR,[],2);
CRITERIA = I>=4 & simplecells; 
numplots = ceil(sqrt(sum(CRITERIA)));
AllIdxs = ind([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]);
% for cell where DoG model performs better than the Gabor model
for jj = 1:numel(CRITERIA) 
    ii = jj;
    
    % STA 
    STA = Output_List{ii,2};
    normfactor = 0.5/(max(abs(STA(:)))+0.01);
    STA = normfactor*STA + 0.5;
    STA = reshape(STA,[10 10 3]);
    
    % Spatial RF
    im = Output_List{ii,4};
    fim = fftshift(fft2(im));
    
    % Modifying im for plotting purposes
    im = sigmoid(im,lambda,0)-0.5;
   
    if CRITERIA(ii)
        disp(ii);
        figure(plot_counter); subplot(numplots,numplots,count); image(STA); set(gca,'XTick',[],'YTick',[]); axis square; title(num2str(Ratio_of_power(ii),3));
        if ismember(ii,AllIdxs)
            set(gca,'XColor','r','YColor','r');
        end
%         figure(plot_counter+1); subplot(numplots,numplots,count); imagesc(255*(im./(2*max(abs(im(:))))+.5)); colormap(gray(255)); axis square; set(gca,'XTick',[],'YTick',[]);
        figure(plot_counter+2); subplot(numplots,numplots,count); imagesc(abs(fim)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title(num2str(Ratio_of_power(ii),3));
        if ismember(ii,AllIdxs)
            set(gca,'XColor','r','YColor','r');
        end
        count = count + 1;
    end
end
plot_counter = plot_counter + 3;

