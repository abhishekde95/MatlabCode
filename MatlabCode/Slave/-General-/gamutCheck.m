function [in_gamut,gamut_scalars] = gamutCheck(cc, bkgndrgb, M, tail)
% function [in_gamut,gamut_scalars] = gamutCheck(cc, bkgndrgb, M, tail)
%
% Will determine whether a requested cone contrast is physically obtainable on the monitor.
%
%   OUTPUT
%       in_gamut: 1 = in gamut, 0 = not in gamut;
%       gamut_scalars: the largest scalar by which cc can be multiplied by and
%       still fit inside the gamut.
%
%   INPUT
%       bkgndrgb: RGBs of the background in normalized intensity units; cc: the
%       cone contrasts desired as columns (direction and amplitude); M: the 'M'
%       matrix obtained by taking the dot products of the phosphor emission
%       spectra and the cone fundamentals; tail: 'both' means we're checking for
%       a biphasic stimulus (like a grating)and 'single' means we're checking
%       for a monophasic stimulus.
%
% NOTE: Each *column* of cc should be a different stimulus

if ~strcmp(tail,'both') && ~strcmp(tail,'single')
    error('Third argument must be either ''both'' or ''single''');
end

cc_size = size(cc);
if length(cc_size) ~= 2 || ~any(cc_size == 3)
    error('cc must be 3-by-N');
elseif cc_size(1) ~= 3
    cc = cc';
    cc_size = cc_size([2 1]);
end
n_dirs = cc_size(2);

bkgndrgb = bkgndrgb(:);
bkgndlms = M*bkgndrgb;
rgbs = M\(bkgndlms(:,ones(1,n_dirs)).*cc); % rgb is is delta units from bkgnd

scalefactors = bsxfun(@rdivide, [1-bkgndrgb; -bkgndrgb], [rgbs; rgbs]);

% No flipping about the origin if we're just considering monopolar stimuli
if strncmp(tail,'single',1)
    scalefactors(scalefactors < 0) = nan;
end

signs = sign(scalefactors);
scalefactors = abs(scalefactors*(1-2*eps)); % A hack to avoid round off errors
% abs to avoid flipping the polarity

% Do each outer product s_i*r_i', where s_i is the i-th column of the signed scalefactors matrix and
% r_i' is the i-th row from the rgbs' matrix, and stack each resulting matrix into the third
% dimension of collisionpts. Then add bkgndrgb to every row of collisionpts.
collisionpts = bsxfun(@plus, ...
    bsxfun(@times, reshape(signs.*scalefactors,[],1,n_dirs), reshape(rgbs,1,[],n_dirs)), ...
    bkgndrgb');

Lvalid = all(collisionpts>=0 & collisionpts<=1, 2);
gamut_scalars = zeros(1, n_dirs);
for ii = 1:n_dirs
    gamut_scalars(ii) = min(scalefactors(Lvalid(:,ii),ii));
end
in_gamut = gamut_scalars > 1;
