% This script will be for analyzing and visualizing the results from the
% already fit data (Data is fit by the file PopAnalysis_GaborvsDOG_WN.m)
% Author - Abhishek De, 1/19
close all; clearvars;
load Output_ListWN.mat
load DOGAICWN.mat
load DOGBICWN.mat
load DOGerrorWN.mat
load GaborAICWN.mat
load GaborBICWN.mat
load GaborerrorWN.mat

[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_LvsM, spikeIdx_LvsM] = fnamesFromTxt2('LvsM.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');
Input_List = [filename_Lum; filename_LvsM; filename_ColorOpponent];
spikeIdx = [spikeIdx_Lum; spikeIdx_LvsM; spikeIdx_ColorOpponent];
cellIds = [repmat({'Lum'},size(filename_Lum)); repmat({'LvsM'},size(filename_LvsM)); repmat({'ColorOpponent'},size(filename_ColorOpponent))];
LumIds = strcmp(cellIds,'Lum');
LvsMIds = strcmp(cellIds,'LvsM');
ColorOpponentIds = strcmp(cellIds,'ColorOpponent');

N = 100;
numgaborparams = 8;
numDOGparams = 6;
bins = -100:10:300;
GaborBIC = N*log(Gaborerror/N) + numgaborparams*log(N);
DOGBIC = N*log(DOGerror/N) + numDOGparams*log(N);
deltaBIC = DOGBIC-GaborBIC;
GaborAIC = 2*numgaborparams + N*log(Gaborerror/N);
DOGAIC = 2*numDOGparams + N*log(DOGerror/N);
deltaAIC = DOGAIC-GaborAIC;
plot_counter = 1;
figure(plot_counter); set(gcf,'Name','model comparison');
subplot(331); plot(Gaborerror(LumIds),DOGerror(LumIds),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[0 10],'Ylim',[0 10],'TickDir','out'); line([0 10],[0 10]);xlabel('Gabor error'); ylabel('DOG Error'); title('Lum'); hold off;
subplot(332); plot(Gaborerror(LvsMIds),DOGerror(LvsMIds),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[0 10],'Ylim',[0 10],'TickDir','out'); line([0 10],[0 10]);xlabel('Gabor error'); ylabel('DOG Error'); title('LvsM'); hold off;
subplot(333); plot(Gaborerror(ColorOpponentIds),DOGerror(ColorOpponentIds),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[0 10],'Ylim',[0 10],'TickDir','out'); line([0 10],[0 10]);xlabel('Gabor error'); ylabel('DOG Error'); title('ColorOpponent'); hold off;
subplot(334); plot(GaborBIC(LumIds),DOGBIC(LumIds),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[-1000 -100],'Ylim',[-1000 -100],'TickDir','out'); line([-1000 -100],[-1000 -100]);xlabel('Gabor BIC'); ylabel('DOG BIC'); title('Lum'); hold off;
subplot(335); plot(GaborBIC(LvsMIds),DOGBIC(LvsMIds),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[-1000 -100],'Ylim',[-1000 -100],'TickDir','out'); line([-1000 -100],[-1000 -100]);xlabel('Gabor BIC'); ylabel('DOG BIC'); title('LvsM'); hold off;
subplot(336); plot(GaborBIC(ColorOpponentIds),DOGBIC(ColorOpponentIds),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[-1000 -100],'Ylim',[-1000 -100],'TickDir','out'); line([-1000 -100],[-1000 -100]);xlabel('Gabor BIC'); ylabel('DOG BIC'); title('Color Opponent'); hold off;
subplot(337); histogram(deltaBIC,bins); hold on; histogram(deltaBIC(LumIds),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(LumIds)),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('Lum'); set(gca,'TickDir','out');hold off;
subplot(338); histogram(deltaBIC,bins); hold on; histogram(deltaBIC(LvsMIds),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(LvsMIds)),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('LvsM'); set(gca,'TickDir','out');hold off;
subplot(339); histogram(deltaBIC,bins); hold on; histogram(deltaBIC(ColorOpponentIds),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(ColorOpponentIds)),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('ColorOpponent'); set(gca,'TickDir','out');hold off;
plot_counter = plot_counter+1;

x = [ones(numel(filename_Lum),1); 2*ones(numel(filename_LvsM),1); 3*ones(numel(filename_ColorOpponent),1)];
p = anova1(deltaBIC,x,'off');
figure(plot_counter); set(gcf,'Name','Boxplots for Lum, LvsM, ColorOpponent: ANOVA 1 way');
boxplot(deltaBIC,x,'Notch','on','Labels',{'Lum','LvsM','Color Opp'}); title(strcat('deltaBIC',{' '},num2str(p,3))); ylabel('delta BIC'); axis square;
plot_counter = plot_counter + 1;

% Now I will select 25 most negative and 25 most positive deltaBIC cells and see what the STA looks like: just to visually confirm if the fitting is working 
numsubplots = 6;
[~,id] = sort(deltaBIC);
mostneg = id(1:numsubplots^2);
mostpos = id(numel(deltaBIC):-1:numel(deltaBIC)-numsubplots^2+1); 
N = [mostneg; mostpos];
count = 1;
for ii = 1:numel(N)
    tmp_vec_gun = Output_List{N(ii),2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
     
    figure(plot_counter),subplot(numsubplots,numsubplots,count); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    count = count + 1;
    if count > numsubplots^2
        plot_counter = plot_counter + 1;
        count  = 1;
    end
end
plot_counter = plot_counter + 1;    
    
figure(plot_counter);bar([sum(LumIds(mostneg)) sum(LvsMIds(mostneg)) sum(ColorOpponentIds(mostneg)); sum(LumIds(mostpos)) sum(LvsMIds(mostpos)) sum(ColorOpponentIds(mostpos))]);
legend('Lum','LvsM','Color Opp');set(gca,'xticklabel',{'best DOG','best Gabor'}); ylabel('no. of cells');
plot_counter = plot_counter + 1;

%% Selecting all phases from Gabor fits where deltaBIC is positive
bins = -1000:10:1000;
Lumgaborphases = abs(Gaborphase(deltaBIC>0 & LumIds)*180/pi);
LvsMgaborphases = abs(Gaborphase(deltaBIC>0 & LvsMIds)*180/pi);
ColorOpponentgaborphases = abs(Gaborphase(deltaBIC>0 & ColorOpponentIds)*180/pi);
x = [ones(numel(Lumgaborphases),1); 2*ones(numel(LvsMgaborphases),1); 3*ones(numel(ColorOpponentgaborphases),1)];
p = anova1([Lumgaborphases;LvsMgaborphases;ColorOpponentgaborphases],x,'off');
figure(plot_counter); set(gcf,'Name','Phases of Gabor fits');
subplot(221); histogram(Lumgaborphases,bins); xlabel('phase'); title('Lum');
subplot(222); histogram(LvsMgaborphases,bins); xlabel('phase'); title('LvsM');
subplot(223); histogram(ColorOpponentgaborphases,bins); xlabel('phase'); title('ColorOpponent');
subplot(224); boxplot([Lumgaborphases;LvsMgaborphases;ColorOpponentgaborphases],x,'Notch','on','Labels',{'Lum','LvsM','Color Opp'}); title(strcat('phase',{' '},num2str(p,3)));axis square;
plot_counter = plot_counter + 1;
