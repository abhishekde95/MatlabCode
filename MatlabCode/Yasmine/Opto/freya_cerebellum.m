%% Yasmine El-Shamayleh (7/2015)

% Cerebellum extent in Freya's SC chamber

%% DRIVE COORDS AND LOGICAL [X, Y, L)
% 0: no cerebellum
% 1: cerebellum found


coords =[
    0 -6 1;
    0 -3 1;
    0 -0 0;
    0 -1.5 0;
    -4 -4 1;
    +4 -4 0; 
    +2 -4 1
    -5.5 -5 1;
    -5.5 -2 0;
    -5.5 -3 0;
    -5.5 -4 1;
    -5.5 -2.5 0;
    -3.5 -3.5 1;
    -3.5 -5 1;
    -1 -5 1;
    +1 -5 1;
    0  -4.5 1;
    -2.5 -4.5 1;
    -4.5 -4.5 1];

%% PLOT

figure; hold on; box off; axis square;
plot([0 0],[-6 6],'k--');
plot([-6 6],[0 0],'k--');

for i = 1:length(coords)
    if coords(i,3)==1
        plot(coords(i,1), coords(i,2), 'b.', 'Markersize', 30);
    elseif  coords(i,3)==0
        plot(coords(i,1), coords(i,2), 'rx', 'Markersize', 30);
    end
end

set(gca,'xlim',[-6 6], 'ylim', [-6 6]);
set(gca,'xtick',-6:1:6, 'ytick', -6:1:6, 'tickDir', 'out');




title('Cerebellum in Freya''s SC chamber')











