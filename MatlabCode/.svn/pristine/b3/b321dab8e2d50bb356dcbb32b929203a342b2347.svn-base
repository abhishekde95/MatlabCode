% Based on the Mathematica code found in this answer:
% http://mathematica.stackexchange.com/a/103557
%
% 2016/03/10 - ZALB

function [x, y, z] = ellipsorand(a, b, c, n)
R = mvnrnd(zeros(1, 3), diag([a^2 b^2 c^3]), n);
XYZ = bsxfun(@rdivide, a*b*c*R, sqrt(R.^2 * [b^2*c^2; a^2*c^2; a^2*b^2]));
x = XYZ(:,1);
y = XYZ(:,2);
z = XYZ(:,3);
