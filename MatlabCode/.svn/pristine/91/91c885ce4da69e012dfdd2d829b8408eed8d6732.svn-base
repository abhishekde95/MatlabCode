function neglogev = computeLogevidence_nse(hmap, W, muf, a, r)

logev = r'*hmap - sum(exp(hmap)) - logdetns(W) - 0.5*(hmap - muf)'*a;

neglogev = - logev;