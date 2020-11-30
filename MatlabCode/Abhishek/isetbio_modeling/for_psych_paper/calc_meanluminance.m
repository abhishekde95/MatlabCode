% calculate the mean luminance of the background color of the CRT monitor 
stro = nex2stro(findfile('G052912008','/Volumes/NO BACKUP/NexFiles/Greg/Psychophysics'));
% stro = nex2stro(findfile('M081417003.nex'));
tmp_spd = reshape(stro.sum.exptParams.mon_spd,101,3);
load('T_xyz1931')
vlambda = T_xyz1931(2,:);
[mon_spd] = SplineSpd([380:4:780]', tmp_spd, [380:5:780]');
lum = 680*vlambda*mon_spd;