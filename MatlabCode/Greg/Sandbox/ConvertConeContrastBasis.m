function out = ConvertConeContrastBasis(Min, Mout, bkgndrgb, cc)
% lmsout = ConvertConeContrastBasis(Min, Mout, bkgndrgb, cc)
%
% Converts contrasts in one basis to contrasts in another basis.
% Useful for looking at NeuroThresh data in cone-contrast spaces 
% defined by different cone fundamentals.
%
% cc should be nx3.

if (size(cc,2) ~= 3)
    cc = cc';
end
if (size(cc,2) ~= 3)
    error('cc should be n x 3');
end
if (size(bkgndrgb,1) ~=3 & size(bkgndrgb,2) ~=3)
    error('bkgndrgb should have 3 elements');
end
if (size(bkgndrgb,1) ~=3)
    bkgndrgb = bkgndrgb';
end

Mconv = Mout/Min;  % matrix that converts cone excitations in stro file
n = size(cc,1);
% irrespective of which fundamentals were used) to 10 deg cone excitations
bkgndlms = Min*bkgndrgb;
lms = cc.*repmat(bkgndlms',n,1)+repmat(bkgndlms',n,1);
lms_out = lms*Mconv';
bkgndlms_out = Mconv*bkgndlms;

out = (lms_out-repmat(bkgndlms_out',n,1))./repmat(bkgndlms_out',n,1);
end


% % Debugging stuff
% load('Dell4BitsCal');
% cal = cals{end};
% load('T_cones_smj10');
% load('T_cones_smj');
% mon_spd = SplineRaw([380:2:780]', cal.P_device, [380:5:780]');
% M1 = T_cones_smj*mon_spd;
% M2 = T_cones_smj10*mon_spd;
% 
% 
% bkgndrgb = [.5 .5 .5];
% cc = [.1 0 -.05]';
% 
% % Making a stimulus in 2 deg space
% lmsbkgnd = M1*bkgndrgb';
% lmsstim = cc.*lmsbkgnd+lmsbkgnd
% rgbstim = inv(M1)*lmsstim
% 
% % now figuring out the cone contrast (in 10 deg) from rgbstim
% lmsbkgnd = M2*bkgndrgb';
% lmsstim = M2*rgbstim;
% (lmsstim-lmsbkgnd)./lmsbkgnd
% 
% ConvertConeContrastBasis(M1, M2, bkgndrgb, cc)