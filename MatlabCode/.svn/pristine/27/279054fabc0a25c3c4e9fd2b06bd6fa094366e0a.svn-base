function [out,fittedSGaussian2D,R_square,Pearson_R] = SingleGaussian2Dfit(im)
% A derivative of SingleGaussianfit
% Author - Abhishek De, 4/20
% This function fits Gaussian with variable 2D variances: 6 parameter model

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
    mux = initguess.xoffset;
    muy = initguess.yoffset;
    sigma_x = 2;
    sigma_y = 2;
    theta = 0;
   
    % Doing the fitting
    options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
    A1 = [0 1 0 0 0 0;... % Upper bound on mux
        0 -1 0 0 0 0;... % Lower bound on mux
        0 0 1 0 0 0;... % Upper bound on muy
        0 0 -1 0 0 0;... % Lower bound on muy
        0 0 0 -1 0 0;... % Sigma_x > 0
        0 0 0 1 0 0;... % Sigma_y < max_c
        0 0 0 0 -1 0;... % Sigma_y > 0
        0 0 0 0 1 0]; % Sigma_y < max_c
        
    b = [centerinterval centerinterval centerinterval centerinterval -eps centerinterval -eps centerinterval];
    [x1,fval1,exitflag1] = fmincon(@(x)SG2Dfiterr_AD(x,im),[gaincenter,mux,muy,sigma_x,sigma_y,theta],A1,b,[],[],[],[],@(x)nonlincondition(x),options);
    [x2,fval2,exitflag2] = fmincon(@(x)SG2Dfiterr_AD(x,-1*im),[gaincenter,mux,muy,sigma_x,sigma_y,theta],A1,b,[],[],[],[],@(x)nonlincondition(x),options);
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
    out.mux = x(2);
    out.muy = x(3);
    out.sigma_x = x(4);
    out.sigma_y = x(5);
    out.theta = x(6);
    out.exitflag = exitflag;
    out.fval = fval;
    
    % Calculating the fitted Single-Gaussian
    interval = [1:nstixperside];
    centerinterval = ceil(median(interval));
    interval = interval-centerinterval;
    [X, Y] = meshgrid(interval(1:1:end),interval);
    X = X-out.mux; Y = Y+out.muy; % So negative numbers mean down
    a = cos(out.theta)^2/(2*out.sigma_x^2) + sin(out.theta)^2/(2*out.sigma_y^2);
    b = -sin(2*out.theta)^2/(4*out.sigma_x^2) + sin(2*out.theta)^2/(4*out.sigma_y^2);
    c = sin(out.theta)^2/(2*out.sigma_x^2) + cos(out.theta)^2/(2*out.sigma_y^2);
    fittedSGaussian2D = out.gaincenter*exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
    fittedSGaussian2D = fittedSGaussian2D/norm(fittedSGaussian2D(:));
    
    % Calculation of coefficient of determination: R-square
    numerator = fval; %min([sum((fittedSGaussian2D(:)+im(:)).^2) sum((fittedSGaussian2D(:)-im(:)).^2)]);
    denominator = sum(im(:).^2.*(im(:)-mean(im(:))).^2);
    R_square = 1-(numerator/denominator);
    Pearson_R = acos(abs(corr(fittedSGaussian2D(:),im(:))))*180/pi;
    
    function error = SG2Dfiterr_AD(params,im)
        Gcenter = params(1);
        x_offset = params(2);
        y_offset = params(3);
        sigma_X = params(4);
        sigma_Y = params(5);
        THETA = params(6);
        interval = [1:size(im,1)];
        interval = interval-ceil(median(interval));
        [X, Y] = meshgrid(interval,interval);
        X = X-x_offset; Y = Y+y_offset; % So negative numbers mean down
        A = cos(THETA)^2/(2*sigma_X^2) + sin(THETA)^2/(2*sigma_Y^2);
        B = -sin(2*THETA)^2/(4*sigma_X^2) + sin(2*THETA)^2/(4*sigma_Y^2);
        C = sin(THETA)^2/(2*sigma_X^2) + cos(THETA)^2/(2*sigma_Y^2);
        
        S2G = Gcenter*exp(-(A*X.^2 + 2*B*X.*Y + C*Y.^2));
        error = sum(im(:).^2.*(S2G(:)-im(:)).^2);
    end

    function [c,ceq] = nonlincondition(params)
        A_center = params(1);
        sigmaX = params(4);  
        sigmaY = params(5); 
        c(1) = -1*A_center ;
        c(2) = A_center - 200;
        c(3) = -1*sigmaX;
        c(4) = -1*sigmaY;
        ceq = [];
    end
    
end


