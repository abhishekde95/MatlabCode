function [stim,stringified,cc_at_thresh] = ChromCat_generate_stimuli(subjectID, tf, thresh_mult)
if nargin < 3, thresh_mult = 1.75; end

[M,bkgndrgb] = M_from_calibration('Dell4BitsCal.mat');
LMTFdata = load(fullfile(fileparts(which('IsoSampOnline')), 'private', 'data', 'LMTF.mat'));
fpar = LMTFdata.(subjectID).model;
theta = LMTFdata.(subjectID).theta;

f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));

b = abs(f1(tf)).^-1; % rg
a = abs(f2(tf)).^-1; % lum
tmptheta = linspace(0,pi/2,12)';
tmplumrg = [cos(tmptheta) sin(tmptheta)]; % first column is L+M, second is L-M
rotmat = [cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
tmplm = tmplumrg*rotmat; % in [L, M]
thtmp = atan2(tmplm(:,2),tmplm(:,1))-theta; % clockwise rotation from [L,M] to [a,b]
rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
stim = tmplm.*rtmp(:,[1 1])*thresh_mult;

[in_gamut,mult_to_edge] = gamutCheck([stim zeros(size(stim,1),1)]', bkgndrgb, M, 'both');

if ~all(in_gamut)
    warning('Not all requested stimuli are in gamut. Adjust the threshold multiplier to about %.3f', ...
        floor(thresh_mult*min(mult_to_edge)*1e3)/1e3);
end

stringified = deblank(char({
    sprintf('    float <subject>%d[][2] = {', tf)
    sprintf('        {%.9f,%.9f},\n', stim(2:end-1,:)')
    '    };'
    }));

cc_at_thresh = char({
    sprintf('        gl_cc_at_LMthresh = %.9f;', stim(end,1)/thresh_mult)
    sprintf('        gl_cc_at_LvMthresh = %.9f;', stim(1,1)/thresh_mult)
    });
