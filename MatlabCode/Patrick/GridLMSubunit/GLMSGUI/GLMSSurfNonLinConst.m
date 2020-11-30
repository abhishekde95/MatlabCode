function [c,ceq] = GLMSSurfNonLinConst(params,LMvals,nsp,surftype,errortype,extra)

% The condition here is that the first sigma value should be smaller than
% the others, ensuring that the tuning is in the principle direction.

%if strcmp(errortype,'NegativeBinomial')
kappa = params(end);
mu = ComputeNakaRushtonJPW(params,LMvals,surftype);
sig = sqrt(mu + (kappa * mu.^2));
ceq = min(sig.^2-mu) - abs(min(sig.^2-mu))-.0000001;
%c(1) = max(mu-sig.^2-.00001);


%if strcmp(surftype,'surface7')
    c = params(2) - params(3:end-4);
%elseif strcmp(surftype,'surface8')
%    c = params(2) - min(params(3:5));
%end

%ceq = [];

% if nargout > 2
%     gradc = zeros(numel(params),numel(c));
%     if any(c > 0)
%         L = c > 0;
%         gradc(end,L) = -c(L).^2;
%         gradc(end,~L) = 0;
%     end
%     if c(1) > 0
%         gradc(1,3:numel(params)-4) = 1;
%     end
%     if c(2) > 0
%         gradc(2,end) = 1;
%     end
%     gradceq = [];
% end

%GC = [];
%GCeq = [];



