function Output = reversefilter(Input,Niter)

fs = 40000; % Sampling rate
samplingrate  = 1/fs;
f1 = 150; % Low frequency cut off
f2 = 8000; % HIgh frequency cut off

[bl,al] = butter(3,f2/(fs/2),'low'); % 3 poles for 8000 Hz, steeper fall off;
[bh,ah] = butter(1,f1/(fs/2),'high'); % 1 pole for 150 Hz

Hd = dfilt.cascade(dfilt.df1(bh,ah),dfilt.df1(bl,al)); % Cascaded filter
X = Input;
for ii = 1:Niter
    X = X + (Input - filter(Hd,X));
end
Output = X;
end

