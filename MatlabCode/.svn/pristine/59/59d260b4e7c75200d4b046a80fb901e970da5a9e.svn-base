function [differece_in_h, hmap, a, W, sqrtL, L] = objFun_NewtonsMethod_nse(h0, datastruct)

H = datastruct.H;
L = diag(exp(h0));
r = datastruct.r;
muf = datastruct.muf;

sqrtL = sqrt(L);

B = eye(size(L)) + sqrtL*H*sqrtL;

W = chol(B, 'lower');

b = L*h0 + r - exp(h0);

M = W\sqrtL;
sqrtLinvBsqrtL = M'*M;
a = b - sqrtLinvBsqrtL*(H*b + muf);

hmap = H*a + muf;

differece_in_h = norm(hmap - h0);


