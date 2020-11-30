% Figuring out the the brain regions in Apollo's chamber
% Author - Abhishek De, 12/18
close all; clearvars;
% left brain from Atlas
% [A1,slices1,scaling1] = getRegionMask('PHT00','V1');
% [A2,slices2,scaling2] = getRegionMask('PHT00','V2');
% hP1 = viewRegionMask(A1,slices1,scaling1,[],1,1,'blue'); % visualizing V1
% hP2 = viewRegionMask(A2,slices2,scaling2,[],1,2,'red'); % visualizing V2

% right brain from another atlas
[A1,slices1,scaling1] = getRegionMask('LVE00_on_F99','V1');
[A21,slices21,scaling21] = getRegionMask('LVE00_on_F99','V2d'); % dorsal V2
[A22,slices22,scaling22] = getRegionMask('LVE00_on_F99','V2v'); % ventral V2
[A3,slices3,scaling3] = getRegionMask('LVE00_on_F99','V3'); % V3
[A4,slices4,scaling4] = getRegionMask('LVE00_on_F99','V4'); % dorsal V4
hP1 = viewRegionMask(A1,slices1,scaling1,[],1,1,'blue',0.2); % visualizing V1
hP21 = viewRegionMask(A21,slices21,scaling21,[],1,1,'red',0.4); % visualizing V2d
hP22 = viewRegionMask(A22,slices22,scaling22,[],1,1,'red',0.4); % visualizing V2v
hP3 = viewRegionMask(A3,slices3,scaling3,[],1,1,[0 1 1],0.6); % visualizing V3
hP4 = viewRegionMask(A4,slices4,scaling4,[],1,1,'green',0.75); % visualizing V4

% idx = -10; idy = -40; idz = 0;
% [x,y,z] = sph2cart(-pi/6,pi/2,10);
% figure(1); hold on; plot3([x x+idx],[y y+idy],[z z+idz],'Linewidth',2)
% hold off;

