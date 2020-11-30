function [differece_in_f, fmap, a, W, sqrtL, obj] = objFun_NewtonsMethod(f0, datastruct)

K = datastruct.K;
r = datastruct.r;

[g, dg,ddg] = logexp1(f0);

% g(g==0)= 1e-6;
z = (r.*(g.*ddg - dg.^2) - g.^2.*ddg)./g.^2;
% z(z>=0) = -1e-6; % to make diag(L) always positive 
z(z>-1e-6) = -1e-3; % to make diag(L) always positive 
% i don't know why, but some z's are very very small like 1e-28.
L = - diag(z);

muf = datastruct.muf;

sqrtL = sqrt(L);

B = eye(size(L)) + sqrtL*K*sqrtL;

W = chol(B, 'lower');

b = L*f0 + (r./g.*dg - dg);
sqrtLinvBsqrtL = sqrtL*(W'\(W\sqrtL));
a = b - sqrtLinvBsqrtL*(K*b + muf);

fmap = K*a + muf;

% check if obj is increasing
gfmap = logexp1(fmap);
obj = -.5*a'*fmap + r'*log(gfmap) - sum(gfmap);

differece_in_f = sum((fmap - f0).^2)/length(fmap);
