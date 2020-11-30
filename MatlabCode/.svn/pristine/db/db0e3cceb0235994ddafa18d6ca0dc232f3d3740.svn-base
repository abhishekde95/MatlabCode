function [negloglike,dL,H] = negloglike_glm_basis_AD(prs,NL,xconvki,yconvhi,y,dt,refreshRate,L2pen,N_spatialparams)
% Adapted from Ali Weber
% function [negloglike,dL,H] = negloglike_glm_basis(prs,NL,xconvki,yconvhi,y,dt,refreshRate)
%   prs:        vector of parameters, coefficients for basis functions in order
%                     [spatial filter; kprs (stim filter); hprs (post-spike filter); dc]
%   NL:         function handle for nonlinearity
%   xconvki:    stimulus convolved with each filter basis vector,
%                     upsampled to match response sampling rate
%   yconvhi:    response vector convolved with each filter basis vector
%   y:          response vector (zeros and ones)
%   dt:         time scale of y (in frames/stimulus frame)
%   refreshRate: refresh rate of stimulus (frames/sec)
%   L2pen:      penalty on L2 norm of prs vector
%
% gives same output as negloglike_glm_basis_old, neglogli_GLM, and Loss_GLM_logli

%% calculate negative log likelihood


nkbasis = size(xconvki,2); % number of basis functions for k
sprs = prs(1:N_spatialparams); % Spatial params 
kprs = prs(N_spatialparams+1:N_spatialparams+nkbasis); % temporal basis functions weighted by given parameters
hprs = prs(N_spatialparams+nkbasis+1:end-1); %  basis functions weighted by given parameters
dc = prs(end); % dc current (accounts for mean spike rate)

% Need to modify the xconvki
xconvki_mod = [repmat(xconvki*kprs,[1 N_spatialparams]) xconvki*sum(sprs)];
xconvk_dc = xconvki_mod*[sprs; kprs] + dc;
yconvh = yconvhi*hprs; % same as output of spikeconv_mex (line 56, negloglig_GLM)

g = xconvk_dc+yconvh; % g = loglambda for NL = @exp, same as Iinj (line 56, neglogli_GLM)
lambda = feval(NL,g);

negloglike = -y'*g + dt*sum(lambda)/refreshRate + L2pen*prs'*prs;  % negative log likelihood

%% calculate negative gradient
if nargout>1
    dL = zeros(size(prs));
    prsMat = [xconvki_mod yconvhi ones(size(xconvki,1),1)];
    for pr = 1:length(prs)
        dL(pr) = -sum(prsMat(logical(y),pr)) + dt/refreshRate*sum(prsMat(:,pr).*lambda) + L2pen*2*prs(pr);
    end
end

%% calculate negative Hessian
if nargout>2
    H = zeros(length(prs));
    prsMat = [xconvki_mod yconvhi ones(size(xconvki,1),1)];
    for pr1 = 1:length(prs)
        for pr2 = pr1:length(prs)
            H(pr1,pr2) = dt/refreshRate*sum(prsMat(:,pr1).*prsMat(:,pr2).*lambda) + L2pen*2*(pr1==pr2);
            H(pr2,pr1) = H(pr1,pr2);
        end
    end
end
