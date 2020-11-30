%close all
%clear all


%% User-defined variables
%Specify cone weights
%a=-275.6;
a=275.6;
b=209.68;
c=32.64;

%Specify color directions (in cone contrast units)
LMSdir(1,:)=[.025 .025 0];
LMSdir(2,:)=[-.01 .01 0];
LMSdir(3,:)=[0 0 .01];

%Specify staircasing variables
numReversals=7;
thresh=20;
stepSize=0.5;
scale=0.5;

%Specify how many phases to be completed
phases=10;

%Specify bootstrapping variables
nbootiter=100;

%Specify number of grand loops
grandloops=1;


%% Non-user-defined variables
load rig1mon
global M 
M=monspd*fundamentals;
global bkgndLMS 
bkgndLMS=bkgndrgb'*M;
global ConeWeights 
ConeWeights=[a b c];


%% Start grand loop
errorMatrix=[];
planepredrespMatrix=[];
correlationMatrix=[];
SSEmatrix=nan(grandloops,1);
SSEmixmatrix=nan(grandloops,nbootiter);
planeparamsmatrix=nan(grandloops,3);
onecolumn=[];

for q=1:grandloops
clear trialspec


%% 3+ phases in 1 loop
[trialspec]=trialspecs(LMSdir,numReversals,thresh,stepSize,scale,phases);


%% Plot these stimuli in color space and fit plane
[xformmat,scaled,Loog,planeparams,SSE]=plotincolorspace(trialspec,thresh);
planeparamsmatrix(q,:)=((planeparams'*xformmat')*thresh);
SSEmatrix(q)=SSE;


 %% Shuffle residuals to find a distribution of planeSSE/quadSSE; Display in histogram
% SSEmixmatrix(:,q)=NaN(nbootiter,1);
% [SSEmixmatrix] = SSEmixemup(scaled,Loog,xformmat,planeparams,nbootiter,SSEmixmatrix,q);
% 
% onecolumn=[onecolumn;SSEmixmatrix(:,q)];
% figure(4);clf;hold on
% title('Combined SSE Bootstrapping loops')
% hist(onecolumn,100)
% drawnow


%% Plot SSEmatrix in histogram
% figure(5);clf;
% set(gca,'XLim',[0 10])
% hist(SSEmatrix,100)
% title('NeuroThresh Distribution of planeSSE/quadSSE')
% drawnow


%% End grand loop
end


%% Display results

% Correlation Analysis - THIS IS LOCATED IN FILE NINEPOINTS.M
% correlationmean=mean(correlationMatrix);
% correlationstd=std(correlationMatrix);
% fprintf('The mean of the correlation coefficient distribution composed of %u samples is %f\n', q, correlationmean)
% fprintf('The standard deviation of the correlation coefficient distribution composed of %u samples is %f\n', q, correlationstd)


% Estimated error of plane fitting
% errorAngle=acosd((ConeWeights/norm(ConeWeights))*((planeparams'*xformmat')/norm(planeparams'*xformmat'))');
% errorMatrix=[errorMatrix, errorAngle];
% fprintf('The angle of the error in degrees is %f\n', errorAngle)


% SSE Neurothresh Model
% meanSSE=mean(SSEmatrix);
% stdSSE=std(SSEmatrix);
% figure; hold on;
% title('SSE Neurothresh Model')
% hist(SSEmatrix,100)
% legend(['Mean SSE =',num2str(meanSSE),   'Mean std =',num2str(stdSSE)]); 


% SSE Bootstrapping
% meanSSEmix=mean(mean(SSEmixmatrix)); 
% stdSSEmix=std(meanSSEmix);
% % for i=1:mixloops %Displays results for individual loops
% %     figure; hold on;
% %     title(['SSE Bootstrapping loop',num2str(i)])
% %     hist(SSEmixmatrix(:,i),100)
% %     legend(['Mean SSE =',num2str(meanSSEmix(i)),'Mean std =',num2str(stdSSEmix(i))],1); 
% % end
% figure; hold on;
% title('Combined SSE Bootstrapping loops')
% text(4,400,['Mean SSE =',num2str(meanmeanSSEmix)])
% text(4,350,['Mean std =',num2str(meanstdSSEmix)])
% hist(SSEmixmatrix,100)


% Display the means and standard deviations
% fprintf('The mean SSE composed of %u samples is %f\n', q, meanSSE)
% fprintf('The standard deviation of the SSEs composed of %u samples is %f\n', q, stdSSE)
% fprintf('The mean bootstrapped SSE composed of %u samples (each with %u shuffles)  is %f\n', q, nbootiter, meanmeanSSEmix)
% fprintf('The standard deviation of the bootstrapped SSEs composed of %u samples (each with %u shuffles) is %f\n', q, nbootiter, meanstdSSEmix)

% Display how well Planeparams predicts Cone Weights

%% Make ConeWeightsDist movie without actual ConeWeights

% figure; axes; hold on;
% plot3(ConeWeightsDist(:,1),ConeWeightsDist(:,2),ConeWeightsDist(:,3),'bo','MarkerFaceColor','blue')
% hold on; %grid on;
% %plot3(ConeWeights(1),ConeWeights(2),ConeWeights(3),'oy','MarkerFaceColor','yellow')
% set(gcf,'Color',[0 0 0]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
% %set(gca,'CameraViewAngleMode','manual')
% set(gca,'Position',[.15 .15 .7 .7]);
% %set(gca,'Color','none');
% set(gcf,'Position',[ 10   321   672   604]);
% %set(gcf,'Position',[100   500   400   0]);

% 
% viewangles = [180:3:360]+90;
% viewangles(end) = [];
% 
% clear M;
% for i = 1:length(viewangles)
%     axis vis3d
%     set(gca,'View',[viewangles(i) 22])
%     M(i) = getframe(gcf);
% end
% 
% repeat = 1;     %default = 1
% pSearch = 1;    %default = 0
% bSearch = 1;    %default = 1
% reference = 1;  %default = 0
% pixRange = 10;  %default = 10
% iFrame = 8;     %default = 8
% pFrame = 10;    %default = 10
% bFrame = 25;    %default = 25
% 
% options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
% mpgwrite(M, gray, 'ConeWeightsPlot1.mpg', options);


%% Make ConeWeightsDist movie with actual ConeWeights
% figure; axes; hold on;
% plot3(ConeWeightsDist(:,1),ConeWeightsDist(:,2),ConeWeightsDist(:,3),'bo','MarkerFaceColor','blue')
% hold on; %grid on;
% plot3(ConeWeights(1),ConeWeights(2),ConeWeights(3),'oy','MarkerFaceColor','yellow')
% set(gcf,'Color',[0 0 0]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','off');
% %set(gca,'CameraViewAngleMode','manual')
% set(gca,'Position',[.15 .15 .7 .7]);
% %set(gca,'Color','none');
% set(gcf,'Position',[ 10   321   672   604]);
% %set(gcf,'Position',[100   500   400   0]);
% 
% 
% viewangles = [0:3:360]+90;
% viewangles(end) = [];
% 
% clear M;
% for i = 1:length(viewangles)
%     axis vis3d
%     set(gca,'View',[viewangles(i) 22])
%     M(i) = getframe(gcf);
% end
% 
% repeat = 1;     %default = 1
% pSearch = 1;    %default = 0
% bSearch = 1;    %default = 1
% reference = 1;  %default = 0
% pixRange = 10;  %default = 10
% iFrame = 8;     %default = 8
% pFrame = 10;    %default = 10
% bFrame = 25;    %default = 25
% 
% options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
% mpgwrite(M, gray, 'ConeWeightsPlot2.mpg', options);
% 
