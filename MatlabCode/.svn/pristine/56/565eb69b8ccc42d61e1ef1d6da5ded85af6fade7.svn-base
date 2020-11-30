function im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg) 
% im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg)    
%
% Returns an image of a Gabor with an (optional) luminance edge running
% through it. No out-of-gamut checking is done so be careful or just
% normalize your image after it's returned by this function.
%
%   INPUTS
%         bkgndrgb: Background RGB values (three-element vector with arguments 0:1)
%         gaborrgb: RGB values of the Gabor (delta from backgroundRGB)
%         edgergb: RGB values of the edge (delta from backgroundRGB)
%         theta: angle of the Gabor in radians
%         lambda: wavelength of the Gabor
%         sigma: standard deviation of the Gabor (in DVA)
%         gamma: aspect ratio of the Gabor (1 = symmetric)
%         phi: phase of the Gabor (in rad)
%         xoff: center of the Gabor in X
%         yoff: center of the Gabor in Y
%         etheta: angle of the edge (relative to the angle of the Gabor)
%         edisp: edge displacement (relative to the center of the Gabor)
%         gausslim: plotinf limits (in Gaussian probability: 0:1)
%         pixperdeg: pixels per degree

    stimsizeindeg = norminv([1-gausslim gausslim],0,1)*abs(sigma)/min([1 gamma]);
    %stimsizeindeg = [min(stimsizeindeg(1),-.65) max(stimsizeindeg(2),.65)];
    stimsizeinpix = round(2*stimsizeindeg(2)*pixperdeg);
    interval = linspace(stimsizeindeg(1), stimsizeindeg(2), stimsizeinpix);
  %  [X, Y] = meshgrid(interval-xoff/2,interval+yoff/2);
    %size(interval)
    [X, Y] = meshgrid(interval-xoff,interval+yoff); 
    xprime = X.*cos(-theta)+Y.*sin(-theta);
    yprime = -X.*sin(-theta)+Y.*cos(-theta);
    fittedgabor = exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos(2.*pi.*yprime./lambda-phi);
    orientation = theta+etheta;
    if (mod(orientation, pi/2) < 1e-5 || mod (orientation, pi/2) > pi/2-1e-5)
        orientation = pi/2*round(orientation/(pi/2));
    end
    edgeim = cosd(orientation*180/pi)*(Y+edisp*sigma*cosd(orientation*180/pi))+...
            sind(orientation*180/pi)*(X+edisp*sigma*sind(orientation*180/pi));
    edgeim(edgeim == 0) = eps;
    edgeim = sign(edgeim).*exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2));
    im = zeros(stimsizeinpix,stimsizeinpix,3);
    for plane = 1:3
        im(:,:,plane) = (fittedgabor.*gaborrgb(plane))+(edgeim*edgergb(plane))+bkgndrgb(plane);
    end
end
