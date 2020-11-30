%fitting size tuning curves

function [xval, fitval, peak, peakval] = fitPeak(peakmethod, radii, dataset) %calculate radius that gave peak response
if peakmethod == 1 %Radius that gave the peak response
    xval = radii;
    fitval = dataset;
    peakval = max(fitval);
    peak = xval(fitval == peakval);
elseif peakmethod == 2 % Peak of spline fit
    xval = linspace(min(radii),max(radii),100);
    pp = csaps(radii,dataset,.01);
    fitval = fnval(pp,xval);
    peakval = max(fitval);
    peak = xval(fitval == peakval);
elseif peakmethod == 3 % Peak of Gaussian fit
    x = radii;
    y = dataset + 10; %positive values only (will subtrast 10 below)
    [~, ~, ~, ~, xx, yyfit] = fitGaussian(x,y);
    xval = xx;
    fitval = yyfit - 10; %subtracting 10 (see comment above)
    peakval = max(fitval);
    peak = xval(fitval == peakval);
end
end
