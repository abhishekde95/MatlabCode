function oiParams = getoiParams_abhi(fov)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
oiParams.pupilArea = 12.6; % in mm, data from CHass et al; 2016
oiParams.pupilDiamMm = 2*sqrt(oiParams.pupilArea/pi);
oiParams.fieldOfViewDegs = fov;
oiParams.offAxis = false;
oiParams.blur = false; % no blurring due to optics
oiParams.lens = false; % optical lens inactive

end

