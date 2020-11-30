function [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog, opt)
    
% [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog, opt)
%
% Function for fitting surfaces to NeuroThresh data.  This will have to be
% updated if I ever come up with a unifying framework for fitting various
% types of quadratic (and other?) surfaces.  Right now, fitting ellipsoids
% is handled completely differently from fitting "simple" quadratic
% surfaces.
% 
% INPUTS
% Scaled is an n x 3 matrix of staircase termination points.
% Loog is a binary vector indicating which color directions went out of gamut.
%   (1 = OOG, 0 = not OOG) 
% Set opt to 'ellipsoid' to fit an ellipsoid, otherwise a simple quadratic
% is fit.
%
% OUTPUTS
% planeparams = [a b c] where ax+by+cz+1=0;
% quadparams = [a b c d e f] where 
%       ax+by+cz+dx^2+ey^2+fxy+1 = 0  (if fitting a "simple quadratic")
%       ax^2+by^2+cz^2+dxy+exz+fyz-1 = 0  (if fitting an ellipsoid)
% planeSSE and quadSSE are the sum of squared errors associated with the
% fits.
% planeparams and quadparams are passed back in a coordinate frame that
% differs from the original coordinate frame by a transformation captured 
% by the xformmat matrix.  We switch to the new coordinate frame to assist
% the fitting algorithm (by spreading out the points) and creating a
% situation in which z is nearly a function of x and y (justifying no
% quadratic z terms for the "simple quadratic").
%
% GDLH 10/8/10

    if (nargin < 3)
        opt = '';
    end
    if strcmp(opt,'ellipsoid')
        ELLIPSOID = 1;
    else
        ELLIPSOID = 0;
    end 
   
    [v,d] = eig(cov([scaled(~Loog,:); -scaled(~Loog,:)]));
    d = diag(d);
    whtmat = v*diag(sqrt(1./d));
    scaled = scaled*whtmat;
    [th,ph,r] = cart2sph(scaled(:,1),scaled(:,2),scaled(:,3));
 
    errs = zeros(50,50);
    tmp = linspace(0,pi,size(errs,1));

    for i = 1:length(tmp)
        for j = 1:length(tmp)
            [tmpa, tmpb, tmpc] = sph2cart(tmp(i),tmp(j),1);
            predr = -1./(tmpa.*cos(ph).*cos(th)+tmpb.*cos(ph).*tmpc.*sin(th)+sin(ph));
            predr = predr.*geomean(r./abs(predr));
            resid = log(abs(predr))-log(r);
            resid(Loog&(abs(predr) > r)) = 0;
            errs(i,j) = sum(resid.^2);
        end
    end
    [cand_i, cand_j] = ind2sub(size(errs),find(errs(:) < prctile(errs(:),30)));
  %  figure;
  %  imagesc(errs < prctile(errs(:),30))
    % Need to trim down the number of candidates in an intelligent way
    planeparams = [];
    planeSSE = sum(r.^2);
    for i = 1:size(cand_i,1)
        [a, b, c] = sph2cart(tmp(cand_i(i)), tmp(cand_j(i)),1);
        predr = -1./(a.*cos(ph).*cos(th)+b.*cos(ph).*sin(th)+c.*sin(ph));
        shiftfactor = mean(log(abs(predr))-log(r));
        a = a*exp(shiftfactor);
        b = b*exp(shiftfactor);
        c = c*exp(shiftfactor);
        
        orthoabc = scaled -(scaled*([a b c]'./norm([a b c])))*([a b c]./norm([a b c]));
        [~,~,v] = svd(orthoabc);
        tmprotmat = MakeOrtho([[a; b; c;] v(:,[1 2])]);
        tmprotmat = [tmprotmat(:,[2 3]) tmprotmat(:,1)];
        xyz = scaled(:,[1 2 3])*tmprotmat;
        rotabc = tmprotmat\[a;b;c];
        
        % Doing gradient descent for planar fit.
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
        [tmpplaneparams, tmpplaneSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr2(xyz, x, Loog),rotabc,options);
        if (tmpplaneSSE < planeSSE)
            planeparams = tmpplaneparams;
            planeSSE = tmpplaneSSE;
            rotmat = tmprotmat;
        end
    end
 
    % Changing coordinate frames for gradient descent fits of planar and
    % quadratic surfaces.
    % Convention: rotating data points so that z is in the abc
    % direction and x and y are uncorrelated.  HOpe this helps the
    % fitting.
    
    % Rotating orthogonal to planeparams
    % Not sure if this help yet.
    orthoabc = scaled-(scaled*(planeparams./norm(planeparams)))*(planeparams./norm(planeparams))';
    [junk1,junk2,v] = svd(orthoabc);
    tmprotmat = MakeOrtho([planeparams v(:,[1 2])]);
    tmprotmat = [tmprotmat(:,[2 3]) tmprotmat(:,1)];
    rotmat = rotmat*tmprotmat;
    xyz = scaled(:,[1 2 3])*rotmat;
    planeparams = tmprotmat\planeparams;

    if (ELLIPSOID)
        D = [xyz(:,1) .* xyz(:,1),...
            xyz(:,2) .* xyz(:,2),...
            xyz(:,3) .* xyz(:,3),...
            2*xyz(:,1) .* xyz(:,2),...
            2*xyz(:,1) .* xyz(:,3),...
            2*xyz(:,2) .* xyz(:,3)];
        lssoln = (D' * D) \(D' * ones(size(xyz,1),1));
        [quadparams, quadSSE, exitflag] = fminsearch(@(x) surfacefiterr3(xyz,x),lssoln, options);
        disp('Fitting an ellipsoid');
    else % non-ellipsoid quadratic
        initguess = [planeparams;0;0;0];
        [quadparams, quadSSE, exitflag, out4] = fminsearch(@(x) surfacefiterr2(xyz, x, Loog),initguess,options);
    end
    xformmat = whtmat*rotmat;
end
