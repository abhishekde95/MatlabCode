% Takes an RF as an input and gives back few isoeccentricity points
% Author - Abhishek De, 5/18
close all;
clearvars;
plot_counter = 1;
RF_loc = [47 -60]; % X & Y in tens of degrees
angle = 30; % in degrees
RFhorizontal_meridian = [norm(RF_loc) 0];
RFvertical_meridian = [0 -norm(RF_loc)];

[THETA,R] = cart2pol(RF_loc(1),RF_loc(2));
[Rx,Ry] = pol2cart([THETA-angle*pi/180 THETA+angle*pi/180],[norm(RF_loc) norm(RF_loc)]);
all_locs = [RF_loc; RFhorizontal_meridian; RFvertical_meridian;[Rx' Ry']];
disp([all_locs sqrt(all_locs(:,1).^2+all_locs(:,2).^2)]);
