function cs_area=VisAngle2Retina_Area(Eccentricity,Area)
% 2nd order polynomial for macaque (Goodchild et al., 1996)
%Edges of stimulus
Outer=(0.038*(Eccentricity+Area/2)^2)+(4.21*(Eccentricity+Area/2))+0.1;
Inner=(0.038*(Eccentricity-Area/2)^2)+(4.21*(Eccentricity-Area/2))+0.1;
StimDiameter=(Outer-Inner); %in mm
cs_area=pi*(StimDiameter/2)^2; %in mm^2
end