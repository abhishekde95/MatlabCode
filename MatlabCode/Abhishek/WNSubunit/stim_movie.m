%% 
% Section 24
% Stimulus movie
close all; 
nframes = 30;
BandW = 0;
figure;
axes;
clear M;

writerObj = VideoWriter('check_brd.mp4');
open(writerObj);

for i = 1:nframes
    a = uint8(normrnd(128,30,10,10,3));

    if (BandW == 1)
        a(:,:,2) = a(:,:,1); 
        a(:,:,3) = a(:,:,1); 
    end
    image(a);
    colormap(gray)
    axis square;

 
    set(gca,'XTick',[],'YTick',[],'Visible','off');
    axis image;
    frame = getframe(gca);
    writeVideo(writerObj,frame);
end

close(writerObj);
