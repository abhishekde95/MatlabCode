function [RGB,LMS,XYZ] = calcRGB(reflectances,illidx,lo,hi,R,C)
global illuminants T_xyz fundamentals
radiance = [];

for ii = 1:size(reflectances,2)
    radiance = [radiance reflectances(:,ii)*illuminants(illidx,ii)]; % refletance x illumination 
end
LMS = radiance*fundamentals;
XYZ = radiance*T_xyz(:,lo:2:hi)';
XYZ1 = max(XYZ, 0);
XYZ1 = XYZ1/max(XYZ1(:));
XYZ1 = XW2RGBFormat(XYZ1,R,C);
RGB = XYZ2sRGB_exgamma(XYZ1); % Converting XYZ to sRGB format
RGB = max(RGB, 0);
RGB = min(RGB, 1);


end

