% Implementing the zero order reverse filtering - Tao et al. 2017
% Author - Abhishek De, 11/19

close all; clearvars;

% Creating a bandpass butterworth filter
fs = 40000; % Sampling rate 
samplingrate  = 1/fs;
f1 = 150; % Low frequency cut off
f2 = 8000; % HIgh frequency cut off

[bl,al] = butter(3,f2/(fs/2),'low'); % 3 poles for 8000 Hz, steeper fall off;
[bh,ah] = butter(1,f1/(fs/2),'high'); % 1 pole for 150 Hz

Hd = dfilt.cascade(dfilt.df1(bh,ah),dfilt.df1(bl,al)); % Cascaded filter 

% Visualizing the bandpass filter 
fvtool(Hd);

% Need to create a dummy filtered signal - start with a sine wave 
t = 0:samplingrate:0.001;
Input = -sin(t*5500);

% Zero order reverse filtering algorithm 
Niter = 1;
X = Input;

figure(1); plot(t,Input,'k','Linewidth',2); hold on;
for ii = 1:Niter
    X = X + (Input - filter(Hd,X)); 
    plot(t,X./max(X),'color',[1-(ii/10) ii/10 0],'Linewidth',2); drawnow;
end
axis square; set(gca,'Tickdir','out','Xlim',[0 0.001]); xlabel('time (s)'); legend('Input','Output'); hold off;
Output = X;

 

