function err = ellipsefiterr(params, xy, Loog)
% err = ellipsefiterr(params, xy, Loog)
%
% fminsearch only searches over a vector of 3 parameters for an
% ellipsoid fit. 
% params(1) = length of semimajor axis
% params(2) = length of semiminor axis
% params(3) = rotation of ellipse

[angs,r] = cart2pol(xy(:,1),xy(:,2));
predr = zeros(length(r),1);

a = params(1);
b = params(2);
th = params(3);

% R = (b^2-a^2)*cos(2*angs-2*th)+a^2+b^2;
% Q = sqrt(2)*a*b*sqrt(R);
% predr = Q./R; % From Wikipedia
predr = (a*b)./sqrt((b*cos(angs-th)).^2 + (a*sin(angs-th)).^2); %also wikipedia


resid = log(predr)-log(r);
resid(Loog&(predr > r)) = 0;
err = sum(abs(resid));
%err = sum(resid.^2);