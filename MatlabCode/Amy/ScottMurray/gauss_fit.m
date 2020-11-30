function [y, x] = gauss_fit(E,data);

%E = "eccentricity" (what is now radius)
% returns the fit y and variable x
% x(1) = width
% x(2) = horizontal shift
% x(3) = amplitude
% x(4) = vertical shift

% these initial conditions work for "normalized" data (scaled to max value)
[x val] = fminsearch('gauss_fun',[.5 2 .5 .5],optimset('MaxFunEvals', 100000),E,data);
%these initial conditions work for non-normalized data (these were fairly
%randomly chosen. could probably do better.
%[x val] = fminsearch('gauss_fun',[2.5 4 22 -2],optimset('MaxFunEvals', 100000),E,data);

S = x(1); %width (std)
D = x(2);   %shift
A = x(3);   %peak
V = x(4); %vertical shift

y = (A*exp(-(E-D).^2/(2*S.^2)) / (S*sqrt(2*pi)))+V;
