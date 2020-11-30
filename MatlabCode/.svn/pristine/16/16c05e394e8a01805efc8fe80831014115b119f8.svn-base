function out = surfacefiterr3(data, params, Loog)

    % Function for computing the error of a parametric function fit to 
    % 3-D data assuming radial error.  Fitting an ellipsoid:
    %    ax^2+by^2+cz^2+dxy+exz+fyz-1 = 0
    % Fit is accomplished by minimizing the squared difference between 
    % the log of the radii and the log of the predicted radii (assuming
    % multiplicative error, which makes the analysis color space
    % independent).
    %
    % It's a little ugly that we have separate fitting routines (and
    % models) for planar and non-planar surfaces.  Eventually it would be
    % nice to have a single model and fitting routine.
    %
    % GDLH 9/26/10
    
    nparams = length(params);
    if (nparams ~=6)
        error('Wrong number of parameters passed to surfacefiterr3');
    end

    if (nargin < 3)
        Loog = false(size(data,1),1);
    else
        Loog = logical(Loog);
    end
    
    a = params(1);
    b = params(2);
    c = params(3);
    d = params(4);
    e = params(5);
    f = params(6);
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);

    % First converting data into polar coordinates
    [th,ph,r] = cart2sph(-x,-y,-z);
    
    predr2 = 1./(a.*(cos(ph).*cos(th)).^2 +...
                b.*(cos(ph).*sin(th)).^2 +...
                c.*sin(ph).^2 + ...
                2*d.*cos(ph).*cos(th).*cos(ph).*sin(th) +...
                2*e.*cos(ph).*cos(th).*sin(ph) +...
                2*f.*cos(ph).*sin(th).*sin(ph));
            
    if (any(predr2<0))
        disp('negative predicted r');
    %    keyboard
    end
            
    predr = sqrt(abs(predr2));
    %   figure; hold on;
    %   plot(log(predr),log(r),'k.');
            
    resid = log(abs(predr))-log(r);
    if (any(eig([params(1) params(4) params(5); params(4) params(2) params(6); params(5) params(6) params(3)]) < 0))
         out = 10^10;
         disp('non-positive definite matrix');
       %  keyboard
         return
    end
        
    out = sum(resid.^2);
end