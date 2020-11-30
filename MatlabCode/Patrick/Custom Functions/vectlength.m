function [length] = vectlength(pt1,pt2)

% Returns the norm of a vector.  For Matricies, each row is a unique
% vector.

length = sqrt(sum(diff([pt1;pt2]).^2,2));