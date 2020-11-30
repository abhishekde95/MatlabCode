function [out,fittedDOG,R_square,Pearson_R] = DOGfit_wgamma(im,DOGfitparams)

% This function also includes a gamma parameter that allows the standard
% deviations to be different along the X and Y dimension
if (size(im,1) ~= size(im,2))
        error('input matrix must be square');
    end
    nstixperside = size(im,1);
    im = im./norm(im(:));
    
    % Getting initguess parameters
    imsq = im.^2;
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    initguess.xoffset = (interval*sum(imsq)')./sum(imsq(:));
    initguess.yoffset = -(interval*sum(imsq')')./sum(imsq(:));
    

    % Looking at the power spectrum
    powerspectra = abs(fftshift(fft2(im))).^2;
    powerspectra = powerspectra.*(powerspectra > .1*max(powerspectra(:)));
    powerspectra = powerspectra./sum(powerspectra(:));

    [i,j] = meshgrid(interval,interval);
    angles = atan2(i,j)-pi/2;
    amp = sqrt(i.^2+j.^2);
    [peaksfy,peaksfx] = ind2sub(size(powerspectra),find(powerspectra == max(powerspectra(:)),1));
    cosines = powerspectra.*cos(angles);
    sines = powerspectra.*sin(angles);
    if (peaksfy<=ceil(median(interval)))
        meansin = mean(mean(sines(1:centerinterval,:)));
    else
        meansin = mean(mean(sines(centerinterval:nstixperside,:)));
    end
    if (peaksfx<=ceil(median(interval)))
        meancos = mean(mean(cosines(:,1:centerinterval)));
    else
        meancos = mean(mean(cosines(:,centerinterval:nstixperside)));
    end
    theta = atan2(meansin,meancos);
    theta = mod(theta + pi/2,2*pi);
    
    if isempty(DOGfitparams)
        
        % Preparing for fitting
        gaincenter = 1;
        gainsurround = 0.5;
        mux = initguess.xoffset;
        muy = initguess.yoffset;
        centersigma = 2;
        surroundsigma = 5;
        gamma = 1;
        
    else
        gaincenter = DOGfitparams.gaincenter;
        gainsurround = DOGfitparams.gaincenter/2;
        mux = DOGfitparams.mux;
        muy = DOGfitparams.muy;
        centersigma = DOGfitparams.centersigma;
        surroundsigma = 2*DOGfitparams.centersigma;
        gamma = 1;

    end
   
    % Doing the fitting
    options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
    A = [0 0 1 0 0 0 0 0;... % Upper bound on mux
        0 0 -1 0 0 0 0 0;... % Lower bound on mux
        0 0 0 1 0 0 0 0;... % Upper bound on muy
        0 0 0 -1 0 0 0 0;... % Lower bound on muy
        0 0 0 0 -1 0 0 0;... % center sigma > 0
        0 0 0 0 1 0 0 0;... % center sigma < max_c
        0 0 0 0 0 -1 0 0;... % surround sigma > 0
        0 0 0 0 0 1 0 0;... % surround sigma < max_s
        0 0 0 0 0 0 -1 0;... % gamma > min_c
        0 0 0 0 0 0 1 0]; % gamma < max_s

    b = [centerinterval centerinterval centerinterval centerinterval -eps centerinterval -eps centerinterval -0.1 10];
    [x1,fval1,exitflag1] = fmincon(@(x)DOGfiterrgamma_AD(x,im),[gaincenter,gainsurround,mux,muy,centersigma,surroundsigma,gamma,theta],A,b,[],[],[],[],@(x)nonlincondition(x),options);
    [x2,fval2,exitflag2] = fmincon(@(x)DOGfiterrgamma_AD(x,-1*im),[gaincenter,gainsurround,mux,muy,centersigma,surroundsigma,gamma,theta],A,b,[],[],[],[],@(x)nonlincondition(x),options);
    if fval1<fval2
        x = x1;
        fval = fval1;
        exitflag = exitflag1;
    else
        x = x2;
        fval = fval2;
        exitflag = exitflag2;
    end
    % Start of parameter
    out.gaincenter = x(1);
    out.gainsurround = x(2);
    out.mux = x(3);
    out.muy = x(4);
    out.centersigma = x(5);
    out.surroundsigma = x(6);
    out.gamma = x(7);
    out.theta = x(8);
    % End of parameter
    
    out.exitflag = exitflag;
    out.fval = fval;
    % Calculating the fitted DOG
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    [X, Y] = meshgrid(interval(1:1:end),interval);
    X = X-out.mux; Y = Y+out.muy; % So negative numbers mean down
    xprime = X.*cos(-out.theta)+Y.*sin(-out.theta);
    yprime = -X.*sin(-out.theta)+Y.*cos(-out.theta);
    fittedDOG = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(xprime.^2+out.gamma*yprime.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(xprime.^2+out.gamma*yprime.^2)/out.surroundsigma^2);
    fittedDOG = fittedDOG/norm(fittedDOG(:));
    
    % Calculation of coefficient of determination: R-square
    numerator = fval; %min([sum((fittedDOG(:)+im(:)).^2) sum((fittedDOG(:)-im(:)).^2)]);
    denominator = sum(im(:).^2.*(im(:)-mean(im(:))).^2);
    R_square = 1-(numerator/denominator);
    Pearson_R = acos(abs(corr(fittedDOG(:),im(:))))*180/pi;
    
    function error = DOGfiterrgamma_AD(params,im)
        Gcenter = params(1);
        Gsurround = params(2);
        x_offset = params(3);
        y_offset = params(4);
        sigmaC = params(5);
        sigmaS = params(6);
        Gamma = params(7);
        Theta = params(8);
        
        interval = [1:size(im,1)];
        interval = interval-ceil(median(interval));
        [X, Y] = meshgrid(interval,interval);
        X = X-x_offset; Y = Y+y_offset; % So negative numbers mean down
        xprime = X.*cos(-Theta)+Y.*sin(-Theta);
        yprime = -X.*sin(-Theta)+Y.*cos(-Theta);
        DOG = Gcenter*sqrt(1/(2*pi*sigmaC^2))*exp(-(xprime.^2+Gamma*yprime.^2)/sigmaC^2) - Gsurround*sqrt(1/(2*pi*sigmaS^2))*exp(-(xprime.^2+Gamma*yprime.^2)/sigmaS^2);
        error = sum(im(:).^2.*(DOG(:)-im(:)).^2);
    end

    function [c,ceq] = nonlincondition(params)
        A_center = params(1);
        A_surround = params(2);
        sigmaC = params(5);
        sigmaS = params(6);
        c(1) = A_surround - A_center ;
        c(2) = -1*A_center;
        c(3) = -1*A_surround;
        c(4) = A_surround - 200;
        c(5) = A_center - 200;
        c(6) = sigmaC - sigmaS;
        ceq = [];
    end
    
end
