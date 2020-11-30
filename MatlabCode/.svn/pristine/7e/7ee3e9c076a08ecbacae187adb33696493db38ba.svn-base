function im = returnfittedimage(out,var)

interval = [1:10];
centerinterval = ceil(median(interval));
interval = interval-centerinterval;
[X, Y] = meshgrid(interval(1:1:end),interval(1:1:end));

if strcmp(var,'Crescent')
    X1 = X-out.muxc; Y1 = Y+out.muyc; % For center, So negative numbers mean down
    X2 = X-out.muxs; Y2 = Y+out.muys; % For surround, So negative numbers mean down
    im = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(X1.^2+Y1.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(X2.^2+Y2.^2)/out.surroundsigma^2);
elseif strcmp(var,'Gabor')
    X = X-out.xoffset; Y = Y+out.yoffset; % So negative numbers mean down
    xprime = X.*cos(-out.theta)+Y.*sin(-out.theta);
    yprime = -X.*sin(-out.theta)+Y.*cos(-out.theta);
    im = out.amplitude*exp(-(xprime.^2+out.gamma.^2.*yprime.^2)./(2.*out.sigma.^2)).*cos((2.*pi.*yprime./out.lambda)-out.phi);   
elseif strcmp(var,'DOG')
    X = X-out.mux; Y = Y+out.muy; 
    im = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(X.^2+Y.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(X.^2+Y.^2)/out.surroundsigma^2);
end

end

