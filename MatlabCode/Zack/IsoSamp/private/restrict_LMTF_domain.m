% When the user enters in values to restrict either r, theta, or z, those ranges
% are passed as `restrict` to this function. The goal is to take those
% restrictions and return the proper wedge of the LMTF surface. Right now,
% restrict(3:4) are the lower and upper theta bounds. Unfortunately, this
% function doesn't return a proper wedge given a restriction on theta.
%
% This function is _not_ relevant when special_case > 0. In those cases r and
% theta are fixed, and so any changes this function makes are ignored.

% Destructively modifies "domain" to only include stimuli within "restrict". 
% If there is a restriction on theta, makes a meshgrid of 50x50 (hardcoded) 
% points in theta, TF, computes the corresponding radii (distance to the surface
% fit), and transforms (theta, TF, r) to (Lcc, Mcc) to find limits on the domain.
% Comment says that it doesn't work for theta. The reason this function exists 
% is that the user specifies stimulus limitations ("restrict") in cylindrical
%coordinates, but the "domain" is in Cartesian coordinates (and depends on 
% the thresholds, which come from the model fit). The reason we need the
% "domain" variable at all is because we are generating points on the surface 
% using "isosurface" which requires a meshgrid that is tight around the 
% surface.

function domain = restrict_LMTF_domain(domain, fpar, restrict)
if ~any(restrict(3:6)), return; end

if any(restrict(5:6)) % cut down TF
    domain(3,1) = restrict(5);
    domain(3,2) = restrict(6);
end

% We need to adjust the bounds for L and M but we don't know what the maximal
% threshold is in the restricted theta, TF wedge. So we calculate those
% thresholds and update the bounds on L and M in the domain matrix.
if any(restrict(3:4))
    npoints = 50;
    r_th = linspace(restrict(3), restrict(4), npoints);
    r_tf = linspace(domain(3,1), domain(3,2), npoints);
    [r_TH, r_TF] = meshgrid(r_th, r_tf);
    r_R = LMTF_thresh_from_model(r_TH, r_TF, fpar);
    [Lp,Mp] = pol2cart(r_TH, r_R);
    domain(1,1) = min(Lp(:));
    domain(1,2) = max(Lp(:));
    domain(2,1) = min(Mp(:));
    domain(2,2) = max(Mp(:));
end
