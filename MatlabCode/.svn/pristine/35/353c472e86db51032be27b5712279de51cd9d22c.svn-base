function [out,fittedDOG,R_square] = DOGfit(im)
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

    % Preparing for fitting
    gaincenter = 1;
    gainsurround = 0.5;
    mux = initguess.xoffset;
    muy = initguess.yoffset;
    centersigma = 2;
    surroundsigma = 5;
   
    % Doing the fitting
    options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
    A = [0 0 1 0 0 0;... % Upper bound on mux
        0 0 -1 0 0 0;... % Lower bound on mux
        0 0 1 0 0 0;... % Upper bound on muy
        0 0 -1 0 0 0;... % Lower bound on muy
        0 0 0 0 -1 0;... % center sigma > 0
        0 0 0 0 1 0;... % center sigma < max_c
        0 0 0 0 0 -1;... % surround sigma > 0
        0 0 0 0 0 1]; % surround sigma < max_s

    b = [centerinterval centerinterval centerinterval centerinterval -eps centerinterval -eps centerinterval];
    [x1,fval1,exitflag1] = fmincon(@(x)DOGfiterr_AD(x,im),[gaincenter,gainsurround,mux,muy,centersigma,surroundsigma],A,b,[],[],[],[],@(x)nonlincondition(x),options);
    [x2,fval2,exitflag2] = fmincon(@(x)DOGfiterr_AD(x,-1*im),[gaincenter,gainsurround,mux,muy,centersigma,surroundsigma],A,b,[],[],[],[],@(x)nonlincondition(x),options);
    if fval1<fval2
        x = x1;
        fval = fval1;
        exitflag = exitflag1;
    else
        x = x2;
        fval = fval2;
        exitflag = exitflag2;
    end
    out.gaincenter = x(1);
    out.gainsurround = x(2);
    out.mux = x(3);
    out.muy = x(4);
    out.centersigma = x(5);
    out.surroundsigma = x(6);
    out.exitflag = exitflag;
    out.fval = fval;
    % Calculating the fitted DOG
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    [X, Y] = meshgrid(interval(1:1:end),interval);
    X = X-out.mux; Y = Y+out.muy; % So negative numbers mean down
    fittedDOG = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(X.^2+Y.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(X.^2+Y.^2)/out.surroundsigma^2);
    R_square = 1-fval/sum(im(:).^2);
    
    function error = DOGfiterr_AD(params,im)
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
