function f = computeLogevidence(fmap, W, muf, a, datastruct)

r = datastruct.r;
g = logexp1(fmap);    
logev = r'*log(g) - sum(g) - logdetns(W) - 0.5*(fmap - muf)'*a;
f = - logev;

% if nargout>1
%     theta1 = param(1);
%     theta2 = param(2);
%     
%     dfdmuf =  sum(a); % der. wrt muf
%     
%     dKdtheta1 = 1./theta1*K; 
%     dfdtheta1 = 0.5*a'*dKdtheta1*a -0.5*trace(sqrtLinvBsqrtL*dKdtheta1);
%     
%     dKdtheta2 = 0.5*1./(theta2).^2.*norm_mat.*K;
%     dfdtheta2 =  0.5*a'*dKdtheta2*a -0.5*trace(sqrtLinvBsqrtL*dKdtheta2);
%     
%     df = [dfdmuf; dfdtheta1; dfdtheta2];
%     df = - df; 
%     
% end
