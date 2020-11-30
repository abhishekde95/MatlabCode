
% Mask - subunit movie
close all; 
nframes = 100;
BandW = 0;
figure;
axes;
clear M;

writerObj = VideoWriter('mask_brd.mp4');
open(writerObj);

for i = 1:nframes
    a = 128*ones(10,10,3);
    rnd_num = normrnd(128,30,2,3);
%     keyboard
    for j = 1:3
        a(2:5,5:6,j) = rnd_num(1,j)*ones(4,2); % subunit 1
        a(2:5,7:8,j) = rnd_num(2,j)*ones(4,2); % subunit 2
    end
    a = uint8(a);
    image(a);
    colormap(gray)
    axis square;
    clear a;
 
    set(gca,'XTick',[],'YTick',[],'Visible','off');
    axis image;
    frame = getframe(gca);
    writeVideo(writerObj,frame);
end

close(writerObj);
