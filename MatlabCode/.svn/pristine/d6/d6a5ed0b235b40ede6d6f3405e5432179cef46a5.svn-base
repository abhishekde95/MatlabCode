
function [out,fittedCrescent,R_square] = Crescentfit_AD(im,DOGfitparams)
if (size(im,1) ~= size(im,2))
        error('input matrix must be square');
    end
    nstixperside = size(im,1);
    im = im./norm(im(:));
    imsq = im.^2;
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    
    if isempty(DOGfitparams)
        % Getting initguess parameters
        initguess.xoffset1 = (interval*sum(imsq)')./sum(imsq(:));
        initguess.yoffset1 = -(interval*sum(imsq')')./sum(imsq(:));
        initguess.xoffset2 = (interval*sum(imsq)')./sum(imsq(:));
        initguess.yoffset2 = -(interval*sum(imsq')')./sum(imsq(:));
        
        % Preparing for fitting
        gaincenter = 1;
        gainsurround = 0.5;
        muxc = initguess.xoffset1; % Center coordinates
        muyc = initguess.yoffset1;
        muxs = initguess.xoffset2; % Surround coordinates
        muys = initguess.yoffset2;
        centersigma = 2;
        surroundsigma = 5;
    else
        gaincenter = DOGfitparams.gaincenter;
        gainsurround = DOGfitparams.gainsurround;
        muxc = DOGfitparams.mux; % Center coordinates
        muyc = DOGfitparams.muy;
        muxs = DOGfitparams.mux; % Surround coordinates
        muys = DOGfitparams.muy;
        centersigma = DOGfitparams.centersigma;
        surroundsigma = DOGfitparams.surroundsigma;
    end
   
    % Doing the fitting
    options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
%     keyboard;
    A = [-1 0 0 0 0 0 0 0;... % Gain center > 0
        0 -1 0 0 0 0 0 0;... % Gain surround > 0
        -1 1 0 0 0 0 0 0;... % Gain center > Gain surround
        1 0 0 0 0 0 0 0;... % Upper bound gain center
        0 0 1 0 0 0 0 0;... % Upper bound on muxc
        0 0 -1 0 0 0 0 0;... % Lower bound on muxc
        0 0 0 1 0 0 0 0;... % Upper bound on muyc
        0 0 0 -1 0 0 0 0;... % Lower bound on muyc
        0 0 0 0 1 0 0 0;... % Upper bound on muxs
        0 0 0 0 -1 0 0 0;... % Lower bound on muxs
        0 0 0 0 0 1 0 0;... % Upper bound on muys
        0 0 0 0 0 -1 0 0;... % Lower bound on muys
        0 0 0 0 0 0 -1 0;... % center sigma > 0
        0 0 0 0 0 0 1 0;... % center sigma < max_c
        0 0 0 0 0 0 0 -1;... % surround sigma > 0
        0 0 0 0 0 0 0 1;... % surround sigma < max_s
        0 0 0 0 0 0 1 -1]; % center sigma < surround sigma
    b = [-eps -eps -eps 200 centerinterval centerinterval centerinterval centerinterval centerinterval centerinterval centerinterval centerinterval -eps centerinterval -eps centerinterval eps];
    [x1,fval1,exitflag1] = fmincon(@(x)Crescentfiterr_AD(x,im),[gaincenter,gainsurround,muxc,muyc,muxs,muys,centersigma,surroundsigma],A,b,[],[],[],[],@(x)nonlincondition(x),options);
    [x2,fval2,exitflag2] = fmincon(@(x)Crescentfiterr_AD(x,-1*im),[gaincenter,gainsurround,muxc,muyc,muxs,muys,centersigma,surroundsigma],A,b,[],[],[],[],@(x)nonlincondition(x),options);
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
    out.muxc = x(3);
    out.muyc = x(4);
    out.muxs = x(5);
    out.muys = x(6);
    out.centersigma = x(7);
    out.surroundsigma = x(8);
    out.exitflag = exitflag;
    out.fval = fval;
    % Calculating the fitted DOG
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    [X, Y] = meshgrid(interval(1:1:end),interval);
    X1 = X-out.muxc; Y1 = Y+out.muyc; % For center, So negative numbers mean down
    X2 = X-out.muxs; Y2 = Y+out.muys; % For surround, So negative numbers mean down
    fittedCrescent = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(X1.^2+Y1.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(X2.^2+Y2.^2)/out.surroundsigma^2);
    R_square = 1-fval/sum(im(:).^2);
     
    function error = Crescentfiterr_AD(params,im)
        Gcenter = params(1);
        Gsurround = params(2);
        x_offset1 = params(3);
        y_offset1 = params(4);
        x_offset2 = params(5);
        y_offset2 = params(6);
        sigmaC = params(7);
        sigmaS = params(8);
        
        interval = [1:size(im,1)];
        interval = interval-ceil(median(interval));
        [X, Y] = meshgrid(interval,interval);
        X1 = X-x_offset1; Y1 = Y+y_offset1;
        X2 = X-x_offset2; Y2 = Y+y_offset2;
        Crescent = Gcenter*sqrt(1/(2*pi*sigmaC^2))*exp(-(X1.^2+Y1.^2)/sigmaC^2) - Gsurround*sqrt(1/(2*pi*sigmaS^2))*exp(-(X2.^2+Y2.^2)/sigmaS^2);
        error = sum(im(:).^2.*(Crescent(:)-im(:)).^2);
    end

    function [c,ceq] = nonlincondition(params)
        x_offset1 = params(3);
        y_offset1 = params(4);
        x_offset2 = params(5);
        y_offset2 = params(6);
        sigmaC = params(7);
        sigmaS = params(8);
        c = sqrt((x_offset1-x_offset2)^2 + (y_offset1-y_offset2)^2) - (sigmaC  + sigmaS);
        ceq = [];
        
    end
    
end
