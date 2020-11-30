function ConeCenters = ConeGrid(StimArea,Eccentricity)
% function [ConeCentersX,ConeCentersY,ConeType]=ConeGrid(StimArea,Eccentricity)
% StimArea and Eccentricity must be defined in visual angle
% Calculates Coordinates of ConeCenters is 2D plane expressed in um
% Created Apr_2011 Angueyra

%Area=VisAngle2Retina_Area(Eccentricity,StimArea); %from visual angle^2 to
%mm^2 (really from mm^2 to DVA)
%Area = StimArea; 
%Area = (12.75 * tand(sqrt(StimArea))).^2
Area = (-55.3947 + 13.1579 * sqrt(0.152 * sqrt(StimArea) + 17.7226)).^2;
Density = VisAngle2Retina_Ecc(Eccentricity); %from visual angle to cones/mm^2
nCones=Area*Density; %number of cones
SConeRatio=VisAngle2Retina_SCones(Eccentricity); %from visual angle to Scones ratio (S./(M+L+S))

ConeSize=2*1e-3; %in mm (this actually changes with eccentricity)
ConeSpacing=1/Density; %in mm
SConeSpacing=1/(Density*SConeRatio); %in mm
Rad3Over2 = sqrt(3) / 2;

%random positional noise
noiseVar = ConeSize/4; %arbitrary number for now. Has this been measured?
Rad3Over2 = sqrt(3) / 2;

%SCones
SgridSize = ceil(sqrt(nCones*SConeRatio));
if rem(SgridSize,2) == 1, SgridSize = SgridSize+1; end %easier to work with even numbers
[SX SY] = meshgrid(1:1:SgridSize); %grid is square
n = size(SX,1);
SX = Rad3Over2 * SX + randn(1,1)*noiseVar;
SY = SY + repmat([0 0.5],[n,n/2]);

% Normalized from 0 to 1 and little translation so that overlay with L-M
% cones is more centered
SX = SX ./ SgridSize;
SY = SY  ./ SgridSize;
SX = SX - min(min(SX))/2.2;
SY = SY - min(min(SY))/2.2;

% Rescale to um
SX= SX * sqrt(Area);
SY= SY * sqrt(Area);

% Jitter Cone Centers
rand('seed',1);
SL = length(SX);
SX = SX + noiseVar/Rad3Over2*randn(SL,SL);
SY = SY + noiseVar*randn(SL,SL);

gridSize = ceil(sqrt(nCones));
if rem(gridSize,2) == 1, gridSize = gridSize+1; end %easier to work with even numbers
[X Y] = meshgrid(1:1:gridSize); %grid is square
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

% Normalized from 0 to 1
X = X ./ gridSize;
Y = Y  ./ gridSize;

% Rescale to um
X= X * sqrt(Area);
Y= Y * sqrt(Area);

% Jitter Cone Centers
rand('seed',1);
L = length(X);
X = X + noiseVar/Rad3Over2*randn(L,L);
Y = Y + noiseVar*randn(L,L);

% Nearest neighbour plot
% [vX,vY]=voronoi(X,Y);
% XLimits=[min(min(X))-ConeSpacing*2 max(max(X))+ConeSpacing*2];
% YLimits=[min(min(Y))-ConeSpacing*2 max(max(Y))+ConeSpacing*2];
% vX(vX>XLimits(2))=XLimits(2);
% vX(vX<XLimits(1))=XLimits(1);
% vY(vY>YLimits(2))=YLimits(2);
% vY(vY<YLimits(1))=YLimits(1);
% [v,c]=voronoin([X(:) Y(:)]);
% plot(X,Y,'k.')
% hold on
% plot(vX,vY,'b.-');
% hold off
% xlim([min(min(X))-ConeSpacing*2 max(max(X))+ConeSpacing*2])
% ylim([min(min(Y))-ConeSpacing*2 max(max(Y))+ConeSpacing*2])

%% Check for collisions
% Given Point Spread Function of the eye (especially for short wavelength)
% light, this might not be worth doing given the time it takes for large
% cone grids. 
% I should check if InterConeDistance > ConeSize or maybe have a cone size
% vs eccentricity subfunction
% If ICD = ConeSize, noiseVar should be 0 (hexagonal grid) (e.g. fovea)

% clear xCrash
% clear yCrash tCrash ttCrash Crash
% for i=1:SL
%     for j=1:SL
%         xCrash{i}{j}=find(X >= SX(i,j)-2*ConeSize & X <= SX(i,j)+2*ConeSize)';
%         yCrash{i}{j}=find(Y >= SY(i,j)-2*ConeSize & Y <= SY(i,j)+2*ConeSize)';
%         tCrash{i}{j}=intersect(xCrash{i}{j},yCrash{i}{j});
%     end
%     ttCrash{i}=cell2mat(tCrash{i});
% end
% Crash=cell2mat(ttCrash);
% clear tCrash tCrash
% 
% figure(1)
% clf
% sp1=subplot(1,2,1);
% set(sp1,'Position',[0.02 0.02 0.46 0.98])
% plot(X,Y,'k.')
% hold on
% plot(SX,SY,'bo')
% plot(reshape(SX,SL^2,1),reshape(SY,SL^2,1),'b.')
% for i=1:SL
%     for j=1:SL
%         circle([SX(i,j),SY(i,j)],ConeSize,20,'b-');
%     end
% end
% plot(X(Crash),Y(Crash),'c.','MarkerSize',18)
% for i=1:size(Crash,2)
%     circle([X(Crash(i)),Y(Crash(i))],ConeSize,20,'c-');
% end
% hold off
% xlim([min(min(X))-ConeSpacing/gridSize*2 max(max(X))+ConeSpacing/gridSize*2])
% ylim([min(min(Y))-ConeSpacing/gridSize*2 max(max(Y))+ConeSpacing/gridSize*2])
% axis square
% drawnow
% 
% % Move randomly the crashed cones until they're out of the way of S-cones
% % (are L and M cones going to crash too with realistic numbers?)
% % I could just give it a number of tries to find a position without
% % crashing and if it doesn't find it just remove cone (?)
% cr=ones(size(Crash));
% if noiseVar==0
%     noiseVar=ConeSize/100;
% end
% for i=1:size(Crash,2)
%     while cr(i)
%         X(Crash(i))=X(Crash(i))+(noiseVar/2/Rad3Over2*randn(1,1)*sqrt(Area));
%         Y(Crash(i))=Y(Crash(i))+(noiseVar/2*randn(1,1)*sqrt(Area));
%         % factor of 2 because X is center and not boundary
%         % Crash with S Cones
%         xReCrashS=find(SX >= X(Crash(i))-(2*ConeSize) & SX <= X(Crash(i))+(2*ConeSize));
%         yReCrashS=find(SY >= Y(Crash(i))-(2*ConeSize) & SY <= Y(Crash(i))+(2*ConeSize));
%         % Crash with L-M Cones
%         xReCrash=find(X >= X(Crash(i))-(ConeSize) & X <= X(Crash(i))+(ConeSize));
%         yReCrash=find(Y >= Y(Crash(i))-(ConeSize) & Y <= Y(Crash(i))+(ConeSize));
% %         hold on
% %         circle([X(Crash(i)),Y(Crash(i))],ConeSize,20,'m-');
% %         drawnow  
%         if isempty(intersect(xReCrashS,yReCrashS)) %&& isempty(intersect(xReCrash,yReCrash))
%             cr(i)=0;
%             hold on
%             circle([X(Crash(i)),Y(Crash(i))],ConeSize,20,'k-');
%             hold off
%             drawnow
%         end
%     end 
% end
%%
% Assign typing to L and M Cones
randType=randn(size(X));
MX=X(randType>=0);
MY=Y(randType>=0);
LX=X(randType<0);
LY=Y(randType<0);

% 
% figure(1)
% % sp2=subplot(1,2,2);
% % set(sp2,'Position',[0.53 0.02 0.46 0.98])
% plot(X,Y,'k.')
% hold on
% plot(MX,MY,'g.')
% for i=1:numel(MX)
%     circle([MX(i),MY(i)],ConeSize,20,'g-');
% end
% plot(LX,LY,'r.')
% for i=1:numel(LX)
%     circle([LX(i),LY(i)],ConeSize,20,'r-');
% end
% plot(SX,SY,'bo')
% plot(reshape(SX,SL^2,1),reshape(SY,SL^2,1),'b.')
% for i=1:SL
%     for j=1:SL
%         circle([SX(i,j),SY(i,j)],ConeSize,20,'b-');
%     end
% end
% % plot(X(Crash),Y(Crash),'c.','MarkerSize',18)
% % for i=1:size(Crash,2)
% %     circle([X(Crash(i)),Y(Crash(i))],ConeSize,20,'c-');
% % end
% hold off
% xlim([min(min(X))-ConeSpacing/gridSize*2 max(max(X))+ConeSpacing/gridSize*2])
% ylim([min(min(Y))-ConeSpacing/gridSize*2 max(max(Y))+ConeSpacing/gridSize*2])
% axis square

%% Output
ConeCenters.MCones.X=MX;
ConeCenters.MCones.Y=MY;
ConeCenters.LCones.X=LX;
ConeCenters.LCones.Y=LY;
ConeCenters.SCones.X=reshape(SX,[],1);
ConeCenters.SCones.Y=reshape(SY,[],1);


end

%% Subfunctions
function cs_area=VisAngle2Retina_Area(Eccentricity,Area)
% 2nd order polynomial for macaque (Goodchild et al., 1996)
%Edges of stimulus

% JPW: This function converts distance from the fovea in mm to DVA... Juan
% has this in reverse, I think.  Also, Juan is using variable Area as a 
% distance, not an area.  distance = sqrt(Area), correct?
%Outer=(0.038*(Eccentricity+sqrt(Area)/2)^2)+(4.21*(Eccentricity+sqrt(Area)/2))+0.1; %JPW
%Inner=(0.038*(Eccentricity-sqrt(Area)/2)^2)+(4.21*(Eccentricity-sqrt(Area)/2))+0.1; %JPW
Outer=(0.038*(Eccentricity+Area/2)^2)+(4.21*(Eccentricity+Area/2))+0.1; %JA
Inner=(0.038*(Eccentricity-Area/2)^2)+(4.21*(Eccentricity-Area/2))+0.1; %JA

StimDiameter=(Outer-Inner); %in mm
cs_area=pi*(StimDiameter/2)^2; %in mm^2
end

function cs_density=VisAngle2Retina_Ecc(Eccentricity)
coeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; %for macaque (Goodchild et al., 1996)
cs_density = (coeffs(1)*(exp(coeffs(2)*Eccentricity)))+...
             (coeffs(3)*(exp(coeffs(4)*Eccentricity)))+...
             (coeffs(5)*(exp(coeffs(6)*Eccentricity)));
cs_density = cs_density*1e3; %??? in cones/mm^2 Left out in Goodchild et al., 1996?
end

function SConeRatio=VisAngle2Retina_SCones(Eccentricity)
SConeRatio=0.06;
end

function H=circle(center,radius,NOP,style)
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
H=plot(X,Y,style);
end