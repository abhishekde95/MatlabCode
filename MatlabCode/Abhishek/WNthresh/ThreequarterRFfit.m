function [out,fittedRF] = ThreequarterRFfit(im)
if (size(im,1) ~= size(im,2))
        error('input matrix must be square');
    end
    nstixperside = size(im,1);
    im = im./max(abs(im(:)));

    % Getting initguess parameters
    imsq = im.^2;
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    initguess.xoffset = (interval*sum(imsq)')./sum(imsq(:));
    initguess.yoffset = -(interval*sum(imsq')')./sum(imsq(:));

    
    % Preparing for fitting
    gaincenter = 1;
    gainsurround = 0.5;
    mux = initguess.xoffset;
    muy = initguess.yoffset;
    centersigma = 2;
    surroundsigma = 5;
   
    
    % Doing the fitting
    options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
%     keyboard;
    A = [0 0 0 0 -1 0;...
        0 0 0 0 0 -1;...
        0 0 0 0 1 -1];
    b = [-eps -eps eps];
    [x,fval,exitflag] = fmincon(@(x)ThreequarterRFfitter_AD(x,im),[gaincenter,gainsurround,mux,muy,centersigma,surroundsigma],A,b,[],[],[],[],[],options);
    out.gaincenter = x(1);
    out.gainsurround = x(2);
    out.mux = x(3);
    out.muy = x(4);
    out.centersigma = x(5);
    out.surroundsigma = x(6);
    out.exitflag = exitflag;
    out.fval = fval;
    % Calculating the fitted DOG
    [X, Y] = meshgrid(interval(1:1:end),interval);
    X = X-out.mux; Y = Y+out.muy; % So negative numbers mean down
    fittedRF = out.gaincenter*exp(-(X.^2+Y.^2)/out.centersigma^2) - out.gainsurround*exp(-(X.^2+Y.^2)/out.surroundsigma^2);

    
    function error = ThreequarterRFfitter_AD(params,im)
        Gcenter = params(1);
        Gsurround = params(2);
        x_offset = params(3);
        y_offset = params(4);
        sigmaC = params(5);
        sigmaS = params(6);
        
        interval = [1:size(im,1)];
        interval = interval-ceil(median(interval));
        [X, Y] = meshgrid(interval,interval);
        X = X-x_offset; Y = Y+y_offset; % So negative numbers mean down
        DOG = Gcenter*sqrt(1/(2*pi*sigmaC^2))*exp(-(X.^2+Y.^2)/sigmaC^2) - Gsurround*sqrt(1/(2*pi*sigmaS^2))*exp(-(X.^2+Y.^2)/sigmaS^2);
        error = sum(im(:).^2.*(DOG(:)-im(:)).^2);
    end
    
end


