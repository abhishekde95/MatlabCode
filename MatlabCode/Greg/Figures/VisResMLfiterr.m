function negloglik = VisResMLfiterr(varargin)
% Support function for ML fitting of weights in a whitenoise experiment
% err = VisResMLfiterr(params,stro,whichframes,spikename,maxT)
% The first four parameters are parameters of the Naka-Rushton function:
% (A, exponent, C50, baseline)
% The other parameters are the weights
% For the moment, fitting every weight in two frames irrespective of whether
% pixels are inside or outside the RF.
%
% This is a support function that is called by VisResFigs.m
persistent stro;
persistent whichframes;
persistent spikename;
persistent maxT;
persistent nframes;
persistent counter;
persistent M;

if nargin == 0 % should never get here
    keyboard 
    negloglik = 1
    return
end

if strcmp(varargin{1},'init')
   stro = varargin{2};
   whichframes = varargin{3};
   spikename = varargin{4};
   maxT = varargin{5};
   if nargin > 5
       M = varargin{6};
   else
       M = [];
   end
   nframes = sum(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'num_frames')));
   negloglik = nan;
   counter = 0;
else
    params = varargin{1};
    A = params(1);
    n = params(2);
    C50 = params(3);
    baseline = params(4);
    weights = params(5:end);
    if length(weights)<stro.sum.exptParams.nstixperside.^2*3 % Using color/space separable filters
       weights = weights(1:end-3)*weights(end-2:end)'; % blowing weights back up to full size
    end
    weights = weights(:)./norm(weights(:));
    
    if ~isempty(M) % line below: everything in braces are the "initargs". The first is what we're projecting onto.
       out = getWhtnsStats(stro,maxT,'STPROJmod', {weights, min(whichframes)-1, nframes, [stro.sum.exptParams.nstixperside.^2 3 length(whichframes)] }, spikename,M);
    else
       out = getWhtnsStats(stro,maxT,'STPROJmod', {weights, min(whichframes)-1, nframes, [stro.sum.exptParams.nstixperside.^2 3 length(whichframes)] }, spikename);        
    end
    contrast = out{1};
    responses = out{2};
    
    lambda = A.*(max(contrast,0).^n./(C50.^n+max(contrast,0).^n))+baseline;
    negloglik = -1*sum(responses.*log(lambda)-lambda)
    counter = counter+1;
    counter
end