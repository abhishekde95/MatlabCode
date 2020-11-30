% This code was initially created to work with FlashTutorial by JAngueyra.  
%Isolated as a function to use with Parent_Flash_SigNoise_Script (not 
%recommended due to time of calculation) -JPW 
% Tweaked May_2011 by JPW

function plot_cone_grid(Cones)

figure

% MCones
plot(Cones.Centers.MCones.X,Cones.Centers.MCones.Y,'g.')
hold on
for i=1:numel(Cones.Centers.MCones.X)
%     circle([Cones.Centers.MCones.X(i),Cones.Centers.MCones.Y(i)],sqrt(Cones.CollectingArea/pi),20,'g-');
    circle([Cones.Centers.MCones.X(i),Cones.Centers.MCones.Y(i)],Cones.Size/2,20,'g-');
end
% LCones
plot(Cones.Centers.LCones.X,Cones.Centers.LCones.Y,'r.')
for i=1:numel(Cones.Centers.LCones.X)
%     circle([Cones.Centers.LCones.X(i),Cones.Centers.LCones.Y(i)],sqrt(Cones.CollectingArea/pi),20,'r-');
    circle([Cones.Centers.LCones.X(i),Cones.Centers.LCones.Y(i)],Cones.Size/2,20,'r-');
end
% SCones
plot(Cones.Centers.SCones.X,Cones.Centers.SCones.Y,'b.')
for i=1:numel(Cones.Centers.SCones.X)
%     circle([Cones.Centers.SCones.X(i),Cones.Centers.SCones.Y(i)],sqrt(Cones.CollectingArea/pi),20,'b-');
    circle([Cones.Centers.SCones.X(i),Cones.Centers.SCones.Y(i)],Cones.Size/2,20,'b-');
end
hold off
axis tight
axis square
xlabel('Retina (mm)')
ylabel('Retina (mm)')