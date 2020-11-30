function [c,ceq] = GLMSSurfNonLinConst(params,LMvals,nsp,surftype,errortype,extra)
%function [c,ceq] = GLMSSurfNonLinConst(params)

% The condition here is that the first sigma value should be smaller than
% the others, ensuring that the tuning is in the principle direction.

% 
c(1) = (params(2) - params(3));
c(2) = (params(2) - params(4));
c(3) = (params(2) - params(5));

% c should be less than 0. Used here to make sure first sigma value is the smallest of the 4.
badL = sign(c) == 1;
c = sum(c(badL).^2).^2;

% ceq should be 0. Just setting it to 0.
ceq = 0;






